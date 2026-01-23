use log::info;
use rayon::prelude::*;
use rust_htslib::faidx::Reader as FastaReader;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

#[derive(Debug)]
pub struct FastaIndex {
    fasta_paths: Vec<PathBuf>,
    sequence_to_file: HashMap<String, usize>,
}

impl FastaIndex {
    pub fn build(fasta_files: &[String]) -> Result<Self, String> {
        if fasta_files.is_empty() {
            return Ok(Self {
                fasta_paths: Vec::new(),
                sequence_to_file: HashMap::new(),
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

        Ok(Self {
            fasta_paths,
            sequence_to_file,
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
        let reader = FastaReader::from_path(fasta_path)
            .map_err(|e| format!("Failed to open FASTA '{}': {}", fasta_path.display(), e))?;

        let raw_seq = reader.fetch_seq(seq_name, start, end - 1).map_err(|e| {
            format!(
                "Failed to fetch {seq_name}:{start}-{end} from '{}': {e}",
                fasta_path.display()
            )
        })?;
        let mut seq_vec = raw_seq.to_vec();
        unsafe { libc::free(raw_seq.as_ptr() as *mut std::ffi::c_void) };
        seq_vec.iter_mut().for_each(|b| *b = b.to_ascii_uppercase());
        Ok(seq_vec)
    }
}
