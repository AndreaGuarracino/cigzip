use rust_htslib::faidx::Reader as FastaReader;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::PathBuf;

pub fn collect_sequence_paths(
    mut files: Vec<String>,
    sequence_list: Option<String>,
) -> Result<Vec<String>, String> {
    if let Some(list_path) = sequence_list {
        let contents = fs::read_to_string(&list_path)
            .map_err(|e| format!("Failed to read --sequence-list '{list_path}': {e}"))?;
        for line in contents.lines() {
            let entry = line.trim();
            if entry.is_empty() || entry.starts_with('#') {
                continue;
            }
            files.push(entry.to_string());
        }
    }

    let mut seen = HashSet::new();
    files.retain(|path| seen.insert(path.clone()));
    Ok(files)
}

#[derive(Debug)]
pub struct SequenceIndex {
    fasta_paths: Vec<PathBuf>,
    sequence_to_file: HashMap<String, usize>,
}

impl SequenceIndex {
    pub fn build(fasta_files: &[String]) -> Result<Self, String> {
        let mut fasta_paths = Vec::new();
        let mut sequence_to_file = HashMap::new();

        for (idx, entry) in fasta_files.iter().enumerate() {
            let path = PathBuf::from(entry);
            if !path.exists() {
                return Err(format!("FASTA file '{entry}' not found"));
            }

            let path_str = path
                .to_str()
                .ok_or_else(|| format!("FASTA path contains invalid UTF-8: {entry}"))?
                .to_string();
            let fai_path = PathBuf::from(format!("{path_str}.fai"));

            if !fai_path.exists() {
                if let Err(e) = FastaReader::from_path(&path) {
                    return Err(format!("Failed to prepare FASTA index for '{entry}': {e}"));
                }
            }

            let fai_content = fs::read_to_string(&fai_path).map_err(|e| {
                format!("Failed to read FASTA index '{}': {}", fai_path.display(), e)
            })?;

            for line in fai_content.lines() {
                let fields: Vec<&str> = line.split('\t').collect();
                if let Some(seq_name) = fields.get(0).map(|s| s.trim()).filter(|s| !s.is_empty()) {
                    sequence_to_file.entry(seq_name.to_string()).or_insert(idx);
                }
            }

            fasta_paths.push(path);
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
