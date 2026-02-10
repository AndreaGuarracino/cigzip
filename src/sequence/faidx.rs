use log::info;
use rayon::prelude::*;
use rust_htslib::faidx::Reader as FastaReader;
use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

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

        let reader = FastaReader::from_path(fasta_path).unwrap_or_else(|e| {
            panic!("Failed to open FASTA '{}': {}", fasta_path.display(), e)
        });
        self.entries.insert(fasta_idx, (reader, counter));
        &mut self.entries.get_mut(&fasta_idx).unwrap().0
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
    /// Per-thread LRU-cached FASTA readers. Indexed by rayon::current_thread_index().
    thread_readers: Vec<Mutex<Option<ReaderCache>>>,
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
                thread_readers: Vec::new(),
                readers_per_thread: 0,
            });
        }

        // Phase 1: Parallel - check/create .fai files and read their contents
        let fai_results: Vec<_> = fasta_files
            .par_iter()
            .enumerate()
            .map(|(idx, entry)| -> Result<(usize, PathBuf, String), String> {
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

                Ok((idx, path, fai_content))
            })
            .collect::<Result<Vec<_>, String>>()?;

        // Phase 2: Sequential - build the HashMap (avoids shared mutable state)
        let mut fasta_paths = vec![PathBuf::new(); fasta_files.len()];
        let mut sequence_to_file = HashMap::new();

        for (idx, path, fai_content) in fai_results {
            for line in fai_content.lines() {
                let fields: Vec<&str> = line.split('\t').collect();
                if let Some(seq_name) = fields.first().map(|s| s.trim()).filter(|s| !s.is_empty()) {
                    sequence_to_file.entry(seq_name.to_string()).or_insert(idx);
                }
            }
            fasta_paths[idx] = path;
        }

        let num_slots = rayon::current_num_threads() + 1;

        // Cap open readers per thread to stay within the fd limit.
        // Each open FastaReader holds ~1 FD (BGZF handle).
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
            thread_readers: (0..num_slots).map(|_| Mutex::new(None)).collect(),
            readers_per_thread,
        })
    }

    pub fn fetch_sequence(
        &self,
        seq_name: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, String> {
        if start >= end {
            return Ok(Vec::new());
        }

        let fasta_idx = self
            .sequence_to_file
            .get(seq_name)
            .ok_or_else(|| format!("Sequence '{seq_name}' not found in supplied FASTA files"))?;

        let fasta_path = &self.fasta_paths[*fasta_idx];

        let thread_idx = rayon::current_thread_index()
            .unwrap_or(self.thread_readers.len() - 1);
        let mut slot = self.thread_readers[thread_idx].lock().unwrap();
        let cache = slot.get_or_insert_with(|| ReaderCache::new(self.readers_per_thread));
        let reader = cache.get_or_open(*fasta_idx, fasta_path);

        let mut seq_vec = reader.fetch_seq(seq_name, start, end - 1).map_err(|e| {
            format!(
                "Failed to fetch {seq_name}:{start}-{end} from '{}': {e}",
                fasta_path.display()
            )
        })?;
        seq_vec.iter_mut().for_each(|b| *b = b.to_ascii_uppercase());
        Ok(seq_vec)
    }
}
