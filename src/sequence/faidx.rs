use log::info;
use rayon::prelude::*;
use rust_htslib::faidx::Reader as FastaReader;
use std::collections::HashMap;
use std::ffi::CString;
use std::fmt;
use std::fs;
use std::os::raw::{c_char, c_int, c_void};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

// --- htslib BGZF FFI (htslib is linked into the binary via rust_htslib) ---
// Used for bulk region reads (bgzf_useek + bgzf_read) which are ~14x faster than
// htslib's faidx_fetch_seq64 byte-by-byte path (bgzf_getc + per-char ctype).
extern "C" {
    fn bgzf_open(path: *const c_char, mode: *const c_char) -> *mut c_void;
    fn bgzf_close(fp: *mut c_void) -> c_int;
    fn bgzf_index_load(fp: *mut c_void, bname: *const c_char, suffix: *const c_char) -> c_int;
    fn bgzf_useek(fp: *mut c_void, uoffset: i64, whence: c_int) -> i64;
    fn bgzf_read(fp: *mut c_void, data: *mut c_void, length: usize) -> isize;
    fn bgzf_set_cache_size(fp: *mut c_void, size: c_int);
}

// Plain (uncompressed) FASTA is served with thread-safe `pread` (File::read_at)
// of only the requested byte range into a per-thread reusable scratch buffer.
// We deliberately avoid mmap here: a shared mmap of a multi-GB FASTA faults the
// whole file into RSS and serializes concurrent page faults on the kernel
// mmap_lock, which regresses both memory and high-thread throughput. pread
// copies just the needed bytes, keeps RSS bounded, and has no shared lock.

/// Strip FASTA line breaks (fixed geometry) and uppercase, from a raw byte slice.
fn fill_strip_upper(raw: &[u8], lb: usize, nl: usize, start_col: usize, out: &mut Vec<u8>) {
    let base = out.len();
    let mut i = 0usize;
    let mut col = start_col;
    while i < raw.len() {
        let take = (lb - col).min(raw.len() - i);
        out.extend_from_slice(&raw[i..i + take]);
        i += take;
        col += take;
        if col == lb {
            i += nl;
            col = 0;
        }
    }
    out[base..].make_ascii_uppercase();
}

/// Strip FASTA line breaks and uppercase IN PLACE within `buf`, which holds the raw
/// region bytes (newlines included) starting at column `start_col`. Compacts with a
/// read/write cursor and truncates. Avoids a separate scratch buffer on the plain path.
fn strip_newlines_in_place(buf: &mut Vec<u8>, lb: usize, nl: usize, start_col: usize) {
    if nl == 0 {
        buf.make_ascii_uppercase();
        return;
    }
    let len = buf.len();
    let (mut r, mut w, mut col) = (0usize, 0usize, start_col);
    while r < len {
        let take = (lb - col).min(len - r);
        if w != r {
            buf.copy_within(r..r + take, w); // memmove; w <= r always
        }
        for b in &mut buf[w..w + take] {
            b.make_ascii_uppercase();
        }
        w += take;
        r += take;
        col += take;
        if col == lb {
            r += nl;
            col = 0;
        }
    }
    buf.truncate(w);
}

/// Per-contig layout from the .fai (byte offset of the sequence + line geometry),
/// enough to compute the uncompressed byte range of any [start,end) region.
#[derive(Clone, Copy)]
struct FaiEntry {
    file_idx: usize,
    length: u64,
    offset: u64,
    linebases: u64,
    linewidth: u64,
}

/// LRU cache for FASTA readers to bound file descriptor usage.
/// Each open reader holds ~1 FD (the BGZF handle); with many FASTA files
/// and many threads the total can exceed the process fd limit.
struct ReaderCache {
    entries: HashMap<usize, (FastaReader, u64)>,
    counter: u64,
    capacity: usize,
}

impl ReaderCache {
    fn new(capacity: usize) -> Self {
        Self {
            entries: HashMap::with_capacity(capacity),
            counter: 0,
            capacity,
        }
    }

    fn get_or_open(&mut self, fasta_idx: usize, fasta_path: &Path) -> &mut FastaReader {
        self.counter += 1;
        let counter = self.counter;

        if self.entries.contains_key(&fasta_idx) {
            let entry = self.entries.get_mut(&fasta_idx).unwrap();
            entry.1 = counter;
            return &mut entry.0;
        }

        // Evict least-recently-used reader if at capacity
        if self.entries.len() >= self.capacity {
            let lru_key = *self
                .entries
                .iter()
                .min_by_key(|(_, (_, ts))| *ts)
                .unwrap()
                .0;
            self.entries.remove(&lru_key);
        }

        let reader = FastaReader::from_path(fasta_path)
            .unwrap_or_else(|e| panic!("Failed to open FASTA '{}': {}", fasta_path.display(), e));
        self.entries.insert(fasta_idx, (reader, counter));
        &mut self.entries.get_mut(&fasta_idx).unwrap().0
    }
}

/// LRU cache of raw BGZF handles (with the .gzi loaded) for the bulk read path.
/// Each handle is one FD, same budget as ReaderCache.
struct BgzfCache {
    entries: HashMap<usize, (*mut c_void, u64)>,
    counter: u64,
    capacity: usize,
    /// Reusable raw-read scratch buffer (avoids re-faulting fresh pages per fetch).
    scratch: Vec<u8>,
}
// Only ever touched by the owning rayon thread (one slot per thread).
unsafe impl Send for BgzfCache {}

impl BgzfCache {
    fn new(capacity: usize) -> Self {
        Self {
            entries: HashMap::with_capacity(capacity),
            counter: 0,
            capacity,
            scratch: Vec::new(),
        }
    }

    /// Returns an open BGZF* with the .gzi index loaded, or null on failure.
    fn get_or_open(&mut self, fasta_idx: usize, fasta_path: &Path) -> *mut c_void {
        self.counter += 1;
        let counter = self.counter;

        if let Some(entry) = self.entries.get_mut(&fasta_idx) {
            entry.1 = counter;
            return entry.0;
        }

        if self.entries.len() >= self.capacity {
            let lru_key = *self
                .entries
                .iter()
                .min_by_key(|(_, (_, ts))| *ts)
                .unwrap()
                .0;
            if let Some((fp, _)) = self.entries.remove(&lru_key) {
                unsafe { bgzf_close(fp) };
            }
        }

        let path_c = CString::new(fasta_path.to_string_lossy().as_bytes()).unwrap();
        let mode_c = CString::new("r").unwrap();
        let gzi_c = CString::new(".gzi").unwrap();
        let fp = unsafe { bgzf_open(path_c.as_ptr(), mode_c.as_ptr()) };
        if fp.is_null() {
            return std::ptr::null_mut();
        }
        // Load the .gzi so bgzf_useek can map uncompressed offsets to blocks.
        if unsafe { bgzf_index_load(fp, path_c.as_ptr(), gzi_c.as_ptr()) } < 0 {
            unsafe { bgzf_close(fp) };
            return std::ptr::null_mut();
        }
        // Cache recently-decompressed blocks so nearby fetches don't re-inflate.
        unsafe { bgzf_set_cache_size(fp, 64 * 1024 * 1024) };
        self.entries.insert(fasta_idx, (fp, counter));
        fp
    }
}

impl Drop for BgzfCache {
    fn drop(&mut self) {
        for (_, (fp, _)) in self.entries.drain() {
            unsafe { bgzf_close(fp) };
        }
    }
}

/// Read the soft file-descriptor limit for this process.
fn get_soft_fd_limit() -> usize {
    fs::read_to_string("/proc/self/limits")
        .ok()
        .and_then(|contents| {
            contents
                .lines()
                .find(|l| l.starts_with("Max open files"))
                .and_then(|l| {
                    let token = l.split_whitespace().nth(3)?;
                    if token.eq_ignore_ascii_case("unlimited") {
                        Some(usize::MAX)
                    } else {
                        token.parse().ok()
                    }
                })
        })
        .unwrap_or(1024)
}

pub struct FastaIndex {
    fasta_paths: Vec<PathBuf>,
    sequence_to_file: HashMap<String, usize>,
    /// Per-contig .fai geometry for the bulk read path.
    fai_entries: HashMap<String, FaiEntry>,
    /// Whether each FASTA file has a .gzi (BGZF index) enabling the bulk path.
    has_gzi: Vec<bool>,
    /// Open File per plain (uncompressed) FASTA: regions served via pread (read_at).
    plain_files: Vec<Option<std::fs::File>>,
    /// Per-thread LRU-cached FASTA readers (fallback path). Indexed by rayon thread.
    thread_readers: Vec<Mutex<Option<ReaderCache>>>,
    /// Per-thread LRU-cached BGZF handles (bulk path). Indexed by rayon thread.
    thread_bgzf: Vec<Mutex<Option<BgzfCache>>>,
    readers_per_thread: usize,
}

impl fmt::Debug for FastaIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastaIndex")
            .field("fasta_paths", &self.fasta_paths)
            .field("num_thread_slots", &self.thread_readers.len())
            .field("readers_per_thread", &self.readers_per_thread)
            .finish_non_exhaustive()
    }
}

impl FastaIndex {
    pub fn build(fasta_files: &[String]) -> Result<Self, String> {
        if fasta_files.is_empty() {
            return Ok(Self {
                fasta_paths: Vec::new(),
                sequence_to_file: HashMap::new(),
                fai_entries: HashMap::new(),
                has_gzi: Vec::new(),
                plain_files: Vec::new(),
                thread_readers: Vec::new(),
                thread_bgzf: Vec::new(),
                readers_per_thread: 0,
            });
        }

        // Phase 1: Parallel - check/create .fai files and read their contents
        let fai_results: Vec<_> = fasta_files
            .par_iter()
            .enumerate()
            .map(|(idx, entry)| -> Result<(usize, PathBuf, String, bool), String> {
                let path = PathBuf::from(entry);
                if !path.exists() {
                    return Err(format!("FASTA file '{entry}' not found"));
                }

                let path_str = path
                    .to_str()
                    .ok_or_else(|| format!("FASTA path contains invalid UTF-8: {entry}"))?;
                let fai_path = PathBuf::from(format!("{path_str}.fai"));

                // Create .fai if missing
                if !fai_path.exists() {
                    info!("Creating FASTA index for '{}'...", entry);
                    if let Err(e) = FastaReader::from_path(&path) {
                        return Err(format!("Failed to create FASTA index for '{entry}': {e}"));
                    }
                }

                // Read .fai content
                let fai_content = fs::read_to_string(&fai_path).map_err(|e| {
                    format!("Failed to read FASTA index '{}': {}", fai_path.display(), e)
                })?;

                let has_gzi = PathBuf::from(format!("{path_str}.gzi")).exists();
                Ok((idx, path, fai_content, has_gzi))
            })
            .collect::<Result<Vec<_>, String>>()?;

        // Phase 2: Sequential - build the maps (avoids shared mutable state)
        let mut fasta_paths = vec![PathBuf::new(); fasta_files.len()];
        let mut sequence_to_file = HashMap::new();
        let mut fai_entries = HashMap::new();
        let mut has_gzi = vec![false; fasta_files.len()];

        for (idx, path, fai_content, gzi) in fai_results {
            has_gzi[idx] = gzi;
            for line in fai_content.lines() {
                let fields: Vec<&str> = line.split('\t').collect();
                let name = match fields.first().map(|s| s.trim()).filter(|s| !s.is_empty()) {
                    Some(n) => n,
                    None => continue,
                };
                sequence_to_file.entry(name.to_string()).or_insert(idx);
                // .fai: name, length, offset, linebases, linewidth
                if fields.len() >= 5 {
                    if let (Ok(length), Ok(offset), Ok(linebases), Ok(linewidth)) = (
                        fields[1].parse::<u64>(),
                        fields[2].parse::<u64>(),
                        fields[3].parse::<u64>(),
                        fields[4].parse::<u64>(),
                    ) {
                        if linebases > 0 {
                            fai_entries.entry(name.to_string()).or_insert(FaiEntry {
                                file_idx: idx,
                                length,
                                offset,
                                linebases,
                                linewidth,
                            });
                        }
                    }
                }
            }
            fasta_paths[idx] = path;
        }

        // Open a File per plain (uncompressed) FASTA so regions are served via
        // pread (read_at) of only the requested bytes. No mmap: keeps RSS bounded
        // and avoids the shared mmap_lock under many threads.
        let plain_files: Vec<Option<std::fs::File>> = fasta_paths
            .iter()
            .enumerate()
            .map(|(idx, p)| {
                if has_gzi[idx] {
                    return None; // BGZF -> bulk bgzf path
                }
                let mut magic = [0u8; 2];
                let is_plain = std::fs::File::open(p)
                    .and_then(|mut f| {
                        use std::io::Read;
                        f.read_exact(&mut magic)
                    })
                    .map(|_| magic != [0x1f, 0x8b]) // not gzip/BGZF magic
                    .unwrap_or(false);
                if is_plain {
                    match std::fs::File::open(p) {
                        Ok(f) => {
                            info!("opened plain FASTA '{}' for pread", p.display());
                            Some(f)
                        }
                        Err(_) => None,
                    }
                } else {
                    None
                }
            })
            .collect();

        let num_slots = rayon::current_num_threads() + 1;

        // Cap open readers per thread to stay within the fd limit.
        // Each open reader/handle holds ~1 FD (BGZF handle).
        let fd_limit = get_soft_fd_limit();
        let reserved_fds = 64; // stdin/stdout/stderr, PAF, logs, etc.
        let available = fd_limit.saturating_sub(reserved_fds);
        let readers_per_thread = (available / num_slots).max(1).min(fasta_files.len());

        if readers_per_thread < fasta_files.len() {
            info!(
                "FD limit={}, thread slots={}, max FASTA readers/thread={} (LRU eviction active, {} FASTA files)",
                if fd_limit == usize::MAX { "unlimited".to_string() } else { fd_limit.to_string() },
                num_slots,
                readers_per_thread,
                fasta_files.len(),
            );
        }

        Ok(Self {
            fasta_paths,
            sequence_to_file,
            fai_entries,
            has_gzi,
            plain_files,
            thread_readers: (0..num_slots).map(|_| Mutex::new(None)).collect(),
            thread_bgzf: (0..num_slots).map(|_| Mutex::new(None)).collect(),
            readers_per_thread,
        })
    }

    /// Uncompressed byte offset of base `pos` within a contig.
    #[inline]
    fn byte_offset(e: &FaiEntry, pos: u64) -> u64 {
        e.offset + (pos / e.linebases) * e.linewidth + (pos % e.linebases)
    }

    pub fn fetch_sequence(
        &self,
        seq_name: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, String> {
        let mut out = Vec::new();
        self.fetch_sequence_into(seq_name, start, end, &mut out)?;
        Ok(out)
    }

    /// Like `fetch_sequence` but writes into a caller-provided buffer (reused across
    /// records to avoid per-fetch allocation / page faults). Clears `out` first.
    pub fn fetch_sequence_into(
        &self,
        seq_name: &str,
        start: usize,
        end: usize,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        out.clear();
        if start >= end {
            return Ok(());
        }

        if let Some(e) = self.fai_entries.get(seq_name) {
            let end_clamped = (end as u64).min(e.length);
            if (start as u64) < end_clamped {
                // Plain path (uncompressed FASTA): pread the region, no decompression.
                if matches!(self.plain_files.get(e.file_idx), Some(Some(_)))
                    && self.fetch_plain_into(e, start as u64, end_clamped, out)?
                {
                    return Ok(());
                }
                // BGZF bulk path: .gzi + .fai geometry -> single bgzf_read of the region.
                if self.has_gzi[e.file_idx]
                    && self.fetch_bulk_into(seq_name, e, start as u64, end_clamped, out)?
                {
                    return Ok(());
                }
                out.clear();
            }
        }

        // Fallback: rust_htslib byte-by-byte faidx.
        self.fetch_fallback_into(seq_name, start, end, out)
    }

    /// Plain (uncompressed) FASTA region read via thread-safe pread (File::read_at)
    /// directly into `out` (newlines included), then stripped in place. One buffer,
    /// no scratch: RSS stays at the reused `out`, which matters for multi-MB regions.
    /// Returns true if it filled `out`, false if the caller should fall back.
    fn fetch_plain_into(
        &self,
        e: &FaiEntry,
        start: u64,
        end: u64,
        out: &mut Vec<u8>,
    ) -> Result<bool, String> {
        use std::os::unix::fs::FileExt;
        let f = match self.plain_files.get(e.file_idx) {
            Some(Some(f)) => f,
            _ => return Ok(false),
        };
        let byte_start = Self::byte_offset(e, start);
        let byte_end = Self::byte_offset(e, end);
        let nbytes = (byte_end - byte_start) as usize;

        // pread the raw region (newlines included) straight into `out`. pread is
        // thread-safe (no shared offset, no mmap_lock), so no per-thread locking.
        out.clear();
        out.reserve(nbytes);
        unsafe {
            out.set_len(nbytes);
        }
        let mut g = 0usize;
        while g < nbytes {
            match f.read_at(&mut out[g..], byte_start + g as u64) {
                Ok(0) => break, // EOF (e.g. last line has no trailing newline)
                Ok(r) => g += r,
                Err(ref er) if er.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(er) => {
                    unsafe {
                        out.set_len(0);
                    }
                    return Err(format!(
                        "pread failed for contig at byte {byte_start} (+{g}) in '{}': {er}",
                        self.fasta_paths[e.file_idx].display()
                    ));
                }
            }
        }
        unsafe {
            out.set_len(g);
        }

        let lb = e.linebases as usize;
        let nl = (e.linewidth - e.linebases) as usize;
        let col = (start % e.linebases) as usize;
        strip_newlines_in_place(out, lb, nl, col);
        Ok(true)
    }

    /// Returns true if it filled `out`, false if the caller should fall back.
    fn fetch_bulk_into(
        &self,
        seq_name: &str,
        e: &FaiEntry,
        start: u64,
        end: u64,
        out: &mut Vec<u8>,
    ) -> Result<bool, String> {
        let fasta_path = &self.fasta_paths[e.file_idx];
        let byte_start = Self::byte_offset(e, start);
        let byte_end = Self::byte_offset(e, end);
        let nbytes = (byte_end - byte_start) as usize;

        let thread_idx = rayon::current_thread_index().unwrap_or(self.thread_bgzf.len() - 1);
        let mut slot = self.thread_bgzf[thread_idx].lock().unwrap();
        let cache = slot.get_or_insert_with(|| BgzfCache::new(self.readers_per_thread));
        let fp = cache.get_or_open(e.file_idx, fasta_path);
        if fp.is_null() {
            return Ok(false); // fall back
        }

        // Reuse the per-thread scratch buffer for the raw region bytes (no per-fetch
        // allocation/page-faults). bgzf_read fills it; we only read the [..got] prefix.
        let raw = &mut cache.scratch;
        raw.clear();
        raw.reserve(nbytes);
        let got;
        unsafe {
            raw.set_len(nbytes);
            if bgzf_useek(fp, byte_start as i64, 0) < 0 {
                return Err(format!(
                    "bgzf_useek failed for {seq_name}:{start}-{end} in '{}'",
                    fasta_path.display()
                ));
            }
            let mut g = 0usize;
            while g < nbytes {
                let r = bgzf_read(fp, raw[g..].as_mut_ptr() as *mut c_void, nbytes - g);
                if r <= 0 {
                    break;
                }
                g += r as usize;
            }
            got = g;
        }

        let lb = e.linebases as usize;
        let nl = (e.linewidth - e.linebases) as usize;
        let col = (start % e.linebases) as usize;
        out.reserve((end - start) as usize);
        fill_strip_upper(&raw[..got], lb, nl, col, out);
        Ok(true)
    }

    fn fetch_fallback_into(
        &self,
        seq_name: &str,
        start: usize,
        end: usize,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        let fasta_idx = self
            .sequence_to_file
            .get(seq_name)
            .ok_or_else(|| format!("Sequence '{seq_name}' not found in supplied FASTA files"))?;

        let fasta_path = &self.fasta_paths[*fasta_idx];

        let thread_idx = rayon::current_thread_index().unwrap_or(self.thread_readers.len() - 1);
        let mut slot = self.thread_readers[thread_idx].lock().unwrap();
        let cache = slot.get_or_insert_with(|| ReaderCache::new(self.readers_per_thread));
        let reader = cache.get_or_open(*fasta_idx, fasta_path);

        let seq_vec = reader.fetch_seq(seq_name, start, end - 1).map_err(|e| {
            format!(
                "Failed to fetch {seq_name}:{start}-{end} from '{}': {e}",
                fasta_path.display()
            )
        })?;
        out.extend_from_slice(&seq_vec);
        out.make_ascii_uppercase();
        Ok(())
    }
}
