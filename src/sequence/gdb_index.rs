// Read sequence from a FASTGA GDB (2-bit packed). The `.1gdb` (ONEcode) holds the scaffold/contig/gap
// skeleton; the hidden `.<name>.bps` holds all contig bases packed 2 bits/base (A=0,C=1,G=2,T=3,
// little-end within a byte, each contig byte-aligned at its `boff`). N-gaps are implicit: the coordinate
// holes between consecutive contigs of a scaffold. Fetch = unpack the covered contig spans and emit 'N'
// for gap spans.
//
// The `.bps` is served with thread-safe `pread` (File::read_at), like the plain-FASTA path: the bytes
// stay resident in the OS page cache (the `.bps` is ~1/4 the ASCII FASTA, so a whole human pangenome
// fits in RAM), but they are NOT mmap'd into our address space — so peak RSS stays at the small working
// set (aligners + buffers), matching ALNtoPAF, instead of counting the whole 2-bit file.
use onecode::OneFile;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::os::unix::fs::FileExt;
use std::path::Path;

const ACGT: [u8; 4] = [b'A', b'C', b'G', b'T'];

// Per-thread scratch for the raw 2-bit bytes of one contig span (reused across fetches).
thread_local! {
    static SCRATCH: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

/// One contig: scaffold-relative start (gaps included), base length, byte offset into the `.bps`.
struct Contig {
    sbeg: u64,
    clen: u64,
    boff: u64,
}

struct Scaffold {
    contigs: Vec<Contig>, // sorted by sbeg
}

pub struct GdbIndex {
    bps: File,
    bps_len: u64,
    scaffolds: HashMap<String, Scaffold>,
}

impl fmt::Debug for GdbIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("GdbIndex")
            .field("scaffolds", &self.scaffolds.len())
            .field("bps_bytes", &self.bps_len)
            .finish()
    }
}

impl GdbIndex {
    pub fn build(files: &[String]) -> Result<Self, String> {
        if files.len() != 1 {
            return Err(format!("GDB decode expects exactly one .1gdb file, got {}", files.len()));
        }
        let gdb_path = &files[0];

        // Hidden bases file: dir/.<stem>.bps , where <stem> = filename without the .1gdb/.gdb suffix.
        let p = Path::new(gdb_path);
        let stem = p
            .file_name()
            .and_then(|n| n.to_str())
            .map(|n| n.trim_end_matches(".1gdb").trim_end_matches(".gdb"))
            .ok_or_else(|| format!("bad GDB path: {gdb_path}"))?;
        let bps_path = p.with_file_name(format!(".{stem}.bps"));

        // Parse the ONEcode skeleton once. S=scaffold header, C=contig len, G=gap len.
        let mut one = OneFile::open_read(gdb_path, None, None, 1)
            .map_err(|e| format!("failed to open GDB '{gdb_path}': {e}"))?;
        let mut scaffolds: HashMap<String, Scaffold> = HashMap::new();
        let mut cur_name: Option<String> = None;
        let mut cur_contigs: Vec<Contig> = Vec::new();
        let mut spos: u64 = 0; // position within the current scaffold (gaps included)
        let mut boff: u64 = 0; // running byte offset in .bps (accumulates across ALL contigs)
        loop {
            let lt = one.read_line();
            if lt == '\0' {
                break;
            }
            match lt {
                'S' => {
                    if let Some(name) = cur_name.take() {
                        scaffolds.insert(name, Scaffold { contigs: std::mem::take(&mut cur_contigs) });
                    }
                    // The S string is the full FASTA header; the PAF uses only the first token (PanSN name).
                    let hdr = one.string().unwrap_or("");
                    let name = hdr.split_whitespace().next().unwrap_or("").to_string();
                    if name.is_empty() {
                        return Err(format!("GDB '{gdb_path}' has a scaffold with an empty name"));
                    }
                    cur_name = Some(name);
                    spos = 0;
                }
                'G' => spos += one.int(0) as u64,
                'C' => {
                    let clen = one.int(0) as u64;
                    cur_contigs.push(Contig { sbeg: spos, clen, boff });
                    spos += clen;
                    boff += (clen + 3) / 4; // COMPRESSED_LEN: 4 bases/byte, contig byte-aligned
                }
                _ => {} // 'f' base-frequency, 'M' soft-mask, provenance, etc.
            }
        }
        if let Some(name) = cur_name.take() {
            scaffolds.insert(name, Scaffold { contigs: cur_contigs });
        }
        if scaffolds.is_empty() {
            return Err(format!("GDB '{gdb_path}' contains no scaffolds"));
        }

        let bps = File::open(&bps_path)
            .map_err(|e| format!("failed to open GDB bases '{}': {e}", bps_path.display()))?;
        let bps_len = bps
            .metadata()
            .map_err(|e| format!("stat '{}': {e}", bps_path.display()))?
            .len();
        if boff > bps_len {
            return Err(format!(
                "GDB bases '{}' too short: skeleton needs {} bytes, file has {}",
                bps_path.display(), boff, bps_len
            ));
        }
        Ok(GdbIndex { bps, bps_len, scaffolds })
    }

    pub fn fetch_sequence(&self, seq_name: &str, start: usize, end: usize) -> Result<Vec<u8>, String> {
        let mut out = Vec::new();
        self.fetch_sequence_into(seq_name, start, end, &mut out)?;
        Ok(out)
    }

    pub fn fetch_sequence_into(
        &self,
        seq_name: &str,
        start: usize,
        end: usize,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        out.clear();
        let scaf = self
            .scaffolds
            .get(seq_name)
            .ok_or_else(|| format!("sequence '{seq_name}' not found in GDB"))?;
        // Return exactly end-start bases: contig bases where covered, 'N' everywhere else (gaps, and any
        // span past the last contig). The caller passes PAF (true-length) coordinates; FASTA has literal
        // N's in those spans, so we do too — never clamp to a reconstructed scaffold length.
        let start = start as u64;
        let end = end as u64;
        if start >= end {
            return Ok(());
        }
        out.reserve((end - start) as usize);

        // First contig that can overlap `start` (contigs sorted by sbeg).
        let first = scaf.contigs.partition_point(|c| c.sbeg + c.clen <= start);
        let mut pos = start;
        SCRATCH.with(|cell| -> Result<(), String> {
            let mut buf = cell.borrow_mut();
            for c in &scaf.contigs[first..] {
                if c.sbeg >= end {
                    break;
                }
                // Gap before this contig -> N's.
                while pos < c.sbeg && pos < end {
                    out.push(b'N');
                    pos += 1;
                }
                if pos >= end {
                    break;
                }
                let stop = end.min(c.sbeg + c.clen);
                // Contig-relative offsets [o0, o1); read the covering .bps bytes with one pread.
                let o0 = pos - c.sbeg;
                let o1 = stop - c.sbeg;
                let b0 = o0 / 4;
                let b1 = (o1 - 1) / 4 + 1;
                let nbytes = (b1 - b0) as usize;
                buf.resize(nbytes, 0);
                self.bps
                    .read_exact_at(&mut buf[..nbytes], c.boff + b0)
                    .map_err(|e| format!("GDB .bps read at {}: {e}", c.boff + b0))?;
                for p in o0..o1 {
                    let byte = buf[(p / 4 - b0) as usize];
                    out.push(ACGT[((byte >> ((p & 3) << 1)) & 3) as usize]);
                }
                pos = stop;
            }
            Ok(())
        })?;
        // Trailing gap up to `end`.
        while pos < end {
            out.push(b'N');
            pos += 1;
        }
        Ok(())
    }
}
