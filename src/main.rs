use clap::{Parser, ValueEnum};
use flate2::read::MultiGzDecoder;
#[cfg(debug_assertions)]
use indicatif::ProgressBar;
#[cfg(debug_assertions)]
use indicatif::ProgressStyle;
#[cfg(debug_assertions)]
use lib_tracepoints::{align_sequences_wfa, cigar_ops_to_cigar_string};
use lib_tracepoints::{
    cigar_to_mixed_tracepoints, cigar_to_tracepoints, cigar_to_variable_tracepoints, cigar_to_tracepoints_fastga,
    cigar_to_tracepoints_raw, cigar_to_variable_tracepoints_raw,
    cigar_to_tracepoints_diagonal, cigar_to_mixed_tracepoints_diagonal, cigar_to_variable_tracepoints_diagonal,
    mixed_tracepoints_to_cigar, tracepoints_to_cigar, variable_tracepoints_to_cigar, tracepoints_to_cigar_fastga,
    DistanceMode, MixedRepresentation,
};
#[cfg(debug_assertions)]
use lib_wfa2::affine_wavefront::AffineWavefronts;
use log::{error, info, debug};
use rayon::prelude::*;
use rust_htslib::faidx::Reader as FastaReader;
use std::fs::File;
use std::fmt;
use std::io::{self, BufRead, BufReader};

/// Tracepoint representation type
#[derive(Debug, Clone, ValueEnum)]
enum TracepointType {
    /// Standard tracepoints
    Standard,
    /// Mixed representation (preserves S/H/P/N CIGAR operations)
    Mixed,
    /// Variable tracepoints representation
    Variable,
    /// FastGA tracepoints representation
    Fastga,
}

impl fmt::Display for TracepointType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TracepointType::Standard => write!(f, "standard"),
            TracepointType::Mixed => write!(f, "mixed"),
            TracepointType::Variable => write!(f, "variable"),
            TracepointType::Fastga => write!(f, "fastga"),
        }
    }
}

/// Common options shared between all commands
#[derive(Parser, Debug)]
struct CommonOpts {
    /// PAF file for alignments (use "-" to read from standard input)
    #[arg(short = 'p', long = "paf")]
    paf: String,

    /// Number of threads to use (default: 4)
    #[arg(long = "threads", default_value_t = 4)]
    threads: usize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[arg(short, long, default_value = "0")]
    verbose: u8,
}

#[derive(Parser, Debug)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Compression of alignments
    Compress {
        #[clap(flatten)]
        common: CommonOpts,

        /// Tracepoint type
        #[arg(long = "type", default_value_t = TracepointType::Standard)]
        tp_type: TracepointType,

        /// Max-diff value for tracepoints (default: 32; 100 if type is fastga)
        #[arg(long)]
        max_diff: Option<usize>,
    },
    /// Decompression of alignments
    Decompress {
        #[clap(flatten)]
        common: CommonOpts,

        /// Tracepoint type
        #[arg(long = "type", default_value_t = TracepointType::Standard)]
        tp_type: TracepointType,

        /// FASTA file for query sequences
        #[arg(short = 'q', long = "query-fasta")]
        query_fasta: String,

        /// FASTA file for target sequences
        #[arg(short = 't', long = "target-fasta")]
        target_fasta: String,

        /// Gap penalties in the format mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2 (only used without fastga type, default: 5,8,2,24,1)
        #[arg(long)]
        penalties: Option<String>,

        /// Trace spacing for fastga (only used with fastga type, default: 100)
        #[arg(long)]
        trace_spacing: Option<usize>,
    },
    /// Run debugging mode (only available in debug builds)
    #[cfg(debug_assertions)]
    Debug {
        /// PAF file for alignments (use "-" to read from standard input)
        #[arg(short = 'p', long = "paf")]
        paf: Option<String>,

        /// Number of threads to use (default: 4)
        #[arg(long = "threads", default_value_t = 4)]
        threads: usize,

        /// FASTA file for query sequences
        #[arg(short = 'q', long = "query-fasta")]
        query_fasta: Option<String>,

        /// FASTA file for target sequences
        #[arg(short = 't', long = "target-fasta")]
        target_fasta: Option<String>,

        /// Gap penalties in the format mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2
        #[arg(long, default_value = "5,8,2,24,1")]
        penalties: String,

        /// Max-diff value for tracepoints
        #[arg(long, default_value = "32")]
        max_diff: usize,

        /// Verbosity level (0 = error, 1 = info, 2 = debug)
        #[arg(short, long, default_value = "0")]
        verbose: u8,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments.
    let args = Args::parse();

    match args {
        Args::Compress {
            common,
            tp_type,
            max_diff,
        } => {
            setup_logger(common.verbose);

            // Determine max_diff: use provided value, or default to 100 if fastga, else 32
            let is_fastga = matches!(tp_type, TracepointType::Fastga);
            let max_diff = max_diff.unwrap_or(if is_fastga { 100 } else { 32 });

            info!("Converting CIGAR to {} tracepoints ({}={})", tp_type,  if is_fastga {"trace_spacing"} else {"max_diff"}, max_diff);

            // Set the thread pool size
            rayon::ThreadPoolBuilder::new()
                .num_threads(common.threads)
                .build_global()?;

            // Open the PAF file (or use stdin if "-" is provided).
            let paf_reader = get_paf_reader(&common.paf)?;

            // Process in chunks
            let chunk_size = std::cmp::max(common.threads * 100, 1000);
            let mut lines = Vec::with_capacity(chunk_size);
            for line_result in paf_reader.lines() {
                match line_result {
                    Ok(line) => {
                        if line.trim().is_empty() || line.starts_with('#') {
                            continue;
                        }

                        lines.push(line);

                        if lines.len() >= chunk_size {
                            // Process current chunk in parallel
                            process_compress_chunk(&lines, &tp_type, max_diff);
                            lines.clear();
                        }
                    }
                    Err(e) => return Err(e.into()),
                }
            }

            // Process remaining lines
            if !lines.is_empty() {
                process_compress_chunk(&lines, &tp_type, max_diff);
            }
        }
        Args::Decompress {
            common,
            tp_type,
            query_fasta,
            target_fasta,
            penalties,
            trace_spacing,
        } => {
            setup_logger(common.verbose);
            info!("Converting {} tracepoints to CIGAR", tp_type);

            // Set the thread pool size
            rayon::ThreadPoolBuilder::new()
                .num_threads(common.threads)
                .build_global()?;

            // Determine if we're using fastga based on tp_type
            let is_fastga = matches!(tp_type, TracepointType::Fastga);

            // Validate and apply conditional defaults
            let trace_spacing = if is_fastga {
                trace_spacing.unwrap_or(100)
            } else {
                if trace_spacing.is_some() {
                    error!("--trace-spacing should only be used with --type fastga");
                    std::process::exit(1);
                }
                0 // Not used when not fastga
            };

            let penalties_str = if is_fastga {
                if penalties.is_some() {
                    error!("--penalties should only be used without --type fastga");
                    std::process::exit(1);
                }
                "5,8,2,24,1".to_string()
            } else {
                penalties.unwrap_or_else(|| "5,8,2,24,1".to_string())
            };

            // Parse penalties
            let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = parse_penalties(&penalties_str)?;

            // Open the PAF file (or use stdin if "-" is provided).
            let paf_reader = get_paf_reader(&common.paf)?;

            // Process in chunks
            let chunk_size = std::cmp::max(common.threads * 100, 1000);
            let mut lines = Vec::with_capacity(chunk_size);
            for line_result in paf_reader.lines() {
                match line_result {
                    Ok(line) => {
                        if line.trim().is_empty() || line.starts_with('#') {
                            continue;
                        }

                        lines.push(line);

                        if lines.len() >= chunk_size {
                            // Process current chunk in parallel
                            process_decompress_chunk(
                                &lines,
                                &tp_type,
                                &query_fasta,
                                &target_fasta,
                                mismatch,
                                gap_open1,
                                gap_ext1,
                                gap_open2,
                                gap_ext2,
                                is_fastga,
                                trace_spacing,
                            );
                            lines.clear();
                        }
                    }
                    Err(e) => return Err(e.into()),
                }
            }

            // Process remaining lines
            if !lines.is_empty() {
                process_decompress_chunk(
                    &lines,
                    &tp_type,
                    &query_fasta,
                    &target_fasta,
                    mismatch,
                    gap_open1,
                    gap_ext1,
                    gap_open2,
                    gap_ext2,
                    is_fastga,
                    trace_spacing,
                );
            }
        }
        #[cfg(debug_assertions)]
        Args::Debug {
            paf,
            threads,
            query_fasta,
            target_fasta,
            penalties,
            max_diff,
            verbose,
        } => {
            setup_logger(verbose);
            info!("Debugging");

            let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = parse_penalties(&penalties)?;
            info!(
                "Penalties: {},{},{},{},{}",
                mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
            );

            if let (Some(paf), Some(query_fasta), Some(target_fasta)) =
                (paf, query_fasta, target_fasta)
            {
                info!("PAF file: {}", paf);
                info!("Query FASTA file: {}", query_fasta);
                info!("Target FASTA file: {}", target_fasta);

                // Count total lines for progress bar (if not stdin)
                let total_lines = if paf != "-" {
                    // Quick line count for progress estimation
                    let counter = get_paf_reader(&paf)?;
                    counter.lines().count()
                } else {
                    0 // Unknown size for stdin
                };

                // Create progress bar
                let progress_bar = if total_lines > 0 {
                    let pb = ProgressBar::new(total_lines as u64);
                    pb.set_style(
                        ProgressStyle::default_bar()
                            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
                            .unwrap()
                            .progress_chars("#>-")
                    );
                    pb.set_message("Processing PAF lines");
                    Some(pb)
                } else {
                    // For stdin, use a spinner
                    let pb = ProgressBar::new_spinner();
                    pb.set_style(
                        ProgressStyle::default_spinner()
                            .template(
                                "{spinner:.green} [{elapsed_precise}] {pos} lines processed {msg}",
                            )
                            .unwrap(),
                    );
                    pb.set_message("Processing PAF lines from stdin");
                    Some(pb)
                };

                // Open the PAF file (or use stdin if "-" is provided).
                let paf_reader = get_paf_reader(&paf)?;

                // Set the thread pool size
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()?;

                // Process in chunks
                let chunk_size = 1000; // Or make this configurable
                let mut lines = Vec::with_capacity(chunk_size);
                let mut processed_count = 0;

                for line_result in paf_reader.lines() {
                    match line_result {
                        Ok(line) => {
                            if line.trim().is_empty() || line.starts_with('#') {
                                continue;
                            }

                            lines.push(line);
                            processed_count += 1;

                            if lines.len() >= chunk_size {
                                // Process current chunk in parallel
                                process_debug_chunk(
                                    &lines,
                                    &query_fasta,
                                    &target_fasta,
                                    mismatch,
                                    gap_open1,
                                    gap_ext1,
                                    gap_open2,
                                    gap_ext2,
                                    max_diff,
                                );

                                // Update progress bar
                                if let Some(ref pb) = progress_bar {
                                    pb.set_position(processed_count as u64);
                                }

                                lines.clear();
                            }
                        }
                        Err(e) => return Err(e.into()),
                    }
                }

                // Process remaining lines
                if !lines.is_empty() {
                    process_debug_chunk(
                        &lines,
                        &query_fasta,
                        &target_fasta,
                        mismatch,
                        gap_open1,
                        gap_ext1,
                        gap_open2,
                        gap_ext2,
                        max_diff,
                    );

                    if let Some(ref pb) = progress_bar {
                        if total_lines > 0 {
                            pb.set_position(total_lines as u64);
                        } else {
                            pb.set_position(processed_count as u64);
                        }
                    }
                }

                // Finish progress bar
                if let Some(ref pb) = progress_bar {
                    pb.finish_with_message("Debug processing complete");
                }
            } else {
                // Fallback: run default example if no PAF/FASTA provided.
                info!("No PAF and FASTA files provided, running default example.");

                let query_seq = b"GAACAGAGAAATGGTGGAATTCAAATACAAAAAAACCGCAAAATTAAAAATCTTGCGGCTCTCTGAACTCATTTTCATGAGTGAATTTGGCGGAACGGACGGGACTCGAACCCGCGACCCCCTGCGTGACAGGCAGGTATTCTAACCGACTGAACTACCGCTCCGCCGTTGTGTTCCGTTGGGAACGGGCGAATATTACGGATTTGCCTCACCCTTCGTCAACGGTTTTTCTCATCTTTTGAATCGTTTGCTGCAAAAATCGCCCAAGCCGCTATTTTTAGCGCCTTTTACAGGTATTTATGCCCGCCAGAGGCAGCTTCCGCCCTTCTTCTCCACCAGATCAAGACGGGCTTCCTGAGCTGCAAGCTCTTCATCTGTCGCAAAAACAACGCGTAACTTACTTGCCTGACGTACAATGCGCTGAATTGTTGCTTCACCTTGTTGCTGCTGTGTCTCTCCTTCCATCGCAAAAGCCATCGACGTTTGACCACCGGTCATCG".to_owned();
                let target_seq = b"GAACAGAGAAATGGTGGAATTCAAATACAAAAAAACCGCAAAATTAACCCTTCGTCAACGGTTTTTCTCATCTTTTGAATCGTTTGCTGCAAAAATCGCCCAAGCCGCTATTTTTAGCGCCTTTTACAGGTATTTATGCCCGCCAGAGGCAGCTTCCGCCCTTCTTCTCCACCAGATCAAGACGGGCTTCCTGAGCTGCAAGCTCTTCATCTGTCGCAAAAACAACGCGTAACTTACTTGCCTGACGTACAATGCGCTGAATTGTTGCTTCACCTTGTTGCTGCTGTGTCTCTCCTTCCATCGCAAAAGCCATCGACGTTTGACCACCGGTCATCG".to_owned();

                let a_start = 0;
                let a_end = query_seq.len();
                let b_start = 0;
                let b_end = target_seq.len();

                // Create aligner and configure settings
                let mut aligner = AffineWavefronts::with_penalties_affine2p(
                    0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
                );
                let paf_cigar = align_sequences_wfa(
                    &query_seq[a_start..a_end],
                    &target_seq[b_start..b_end],
                    &mut aligner,
                );
                let paf_cigar = cigar_ops_to_cigar_string(&paf_cigar);

                let tracepoints = cigar_to_tracepoints(&paf_cigar, max_diff);

                let distance_mode = DistanceMode::Affine2p {
                    mismatch,
                    gap_open1,
                    gap_ext1,
                    gap_open2,
                    gap_ext2,
                };
                let cigar_from_tracepoints = tracepoints_to_cigar(
                    &tracepoints,
                    &query_seq,
                    &target_seq,
                    0,
                    0,
                    &distance_mode,
                );

                if false {
                    error!("CIGAR mismatch!");
                    error!("\t                         tracepoints: {:?}", tracepoints);
                    error!("\t                      CIGAR from PAF: {}", paf_cigar);
                    error!(
                        "\t              CIGAR from tracepoints: {}",
                        cigar_from_tracepoints
                    );
                    error!(
                        "\t                      bounds CIGAR from PAF: {:?}",
                        get_cigar_diagonal_bounds(&paf_cigar)
                    );
                    error!(
                        "\t              bounds CIGAR from tracepoints: {:?}",
                        get_cigar_diagonal_bounds(&cigar_from_tracepoints)
                    );

                    let (deviation, d_min, d_max, max_gap) =
                        compute_deviation(&cigar_from_tracepoints);
                    error!(
                        "\t                      deviation CIGAR from PAF: {:?}",
                        compute_deviation(&paf_cigar)
                    );
                    error!(
                        "\t              deviation CIGAR from tracepoints: {:?}",
                        (deviation, d_min, d_max, max_gap)
                    );
                    error!("=> Try using --wfa-heuristic=banded-static --wfa-heuristic-parameters=-{},{}\n", std::cmp::max(max_gap, -d_min), std::cmp::max(max_gap, d_max));
                }
            }
        }
    }

    Ok(())
}

/// Calculate alignment score based on edit distance from a CIGAR string
/// Alignment score = -(mismatches + inserted_bp + deleted_bp)
fn calculate_alignment_score_edit_distance(cigar: &str) -> i32 {
    let mut mismatches = 0;
    let mut inserted_bp = 0;
    let mut deleted_bp = 0;

    let mut num_buffer = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num_buffer.push(c);
        } else {
            let len = num_buffer.parse::<usize>().unwrap_or(0);
            num_buffer.clear();

            match c {
                'X' => mismatches += len,
                'I' => inserted_bp += len,
                'D' => deleted_bp += len,
                _ => {}
            }
        }
    }

    let edit_distance = mismatches + inserted_bp + deleted_bp;
    -(edit_distance as i32)
}

/// Calculate gap-compressed identity and block identity from a CIGAR string
fn calculate_identity_stats(cigar: &str) -> (f64, f64) {
    let mut matches = 0;
    let mut mismatches = 0;
    let mut insertions = 0;
    let mut inserted_bp = 0;
    let mut deletions = 0;
    let mut deleted_bp = 0;

    let mut num_buffer = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num_buffer.push(c);
        } else {
            let len = num_buffer.parse::<usize>().unwrap_or(0);
            num_buffer.clear();

            match c {
                'M' | '=' => matches += len,
                'X' => mismatches += len,
                'I' => {
                    insertions += 1;
                    inserted_bp += len;
                }
                'D' => {
                    deletions += 1;
                    deleted_bp += len;
                }
                'S' | 'H' | 'P' | 'N' => {}
                _ => {}
            }
        }
    }

    let gap_compressed_identity = if matches + mismatches + insertions + deletions > 0 {
        (matches as f64) / (matches + mismatches + insertions + deletions) as f64
    } else {
        0.0
    };

    let edit_distance = mismatches + inserted_bp + deleted_bp;
    let block_identity = if matches + edit_distance > 0 {
        (matches as f64) / (matches + edit_distance) as f64
    } else {
        0.0
    };

    (gap_compressed_identity, block_identity)
}

/// Process a chunk of lines in parallel for debugging
#[cfg(debug_assertions)]
fn process_debug_chunk(
    lines: &[String],
    query_fasta_path: &str,
    target_fasta_path: &str,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
    max_diff: usize,
) {
    lines.par_iter().for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            error!(
                "{}",
                message_with_truncate_paf_file("Skipping malformed PAF line", line)
            );
            std::process::exit(1);
        }

        let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
            error!(
                "{}",
                message_with_truncate_paf_file("Skipping CIGAR-less PAF line", line)
            );
            std::process::exit(1);
        };
        let paf_cigar = &cg_field[5..];

        // Parse mandatory PAF fields.
        let query_name = fields[0];
        let query_start: usize = fields[2].parse().unwrap_or_else(|_| {
            error!("Invalid query_start in PAF line");
            std::process::exit(1);
        });
        let query_end: usize = fields[3].parse().unwrap_or_else(|_| {
            error!("Invalid query_end in PAF line");
            std::process::exit(1);
        });
        let strand = fields[4];
        let target_name = fields[5];
        let target_start: usize = fields[7].parse().unwrap_or_else(|_| {
            error!("Invalid target_start in PAF line");
            std::process::exit(1);
        });
        let target_end: usize = fields[8].parse().unwrap_or_else(|_| {
            error!("Invalid target_end in PAF line");
            std::process::exit(1);
        });

        // Create thread-local FASTA readers
        let query_fasta_reader = FastaReader::from_path(query_fasta_path).unwrap_or_else(|e| {
            error!("Error reading query FASTA file: {}", e);
            std::process::exit(1);
        });
        let target_fasta_reader = FastaReader::from_path(target_fasta_path).unwrap_or_else(|e| {
            error!("Error reading target FASTA file: {}", e);
            std::process::exit(1);
        });

        // Fetch query sequence from query FASTA
        let query_seq = if strand == "+" {
            match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut seq_vec = seq.to_vec();
                    unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) };
                    seq_vec
                        .iter_mut()
                        .for_each(|byte| *byte = byte.to_ascii_uppercase());
                    seq_vec
                }
                Err(e) => {
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
                }
            }
        } else {
            match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut rc = reverse_complement(&seq.to_vec());
                    unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) };
                    rc.iter_mut()
                        .for_each(|byte| *byte = byte.to_ascii_uppercase());
                    rc
                }
                Err(e) => {
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
                }
            }
        };

        // Fetch target sequence from target FASTA
        let target_seq = match target_fasta_reader.fetch_seq(
            target_name,
            target_start,
            target_end - 1,
        ) {
            Ok(seq) => {
                let mut seq_vec = seq.to_vec();
                unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) };
                seq_vec
                    .iter_mut()
                    .for_each(|byte| *byte = byte.to_ascii_uppercase());
                seq_vec
            }
            Err(e) => {
                error!("Failed to fetch target sequence: {}", e);
                std::process::exit(1);
            }
        };

        // Create thread-local aligner
        let mut aligner = AffineWavefronts::with_penalties_affine2p(
            0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
        );
        let realn_cigar = align_sequences_wfa(&query_seq, &target_seq, &mut aligner);
        let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);
        // let paf_cigar = &realn_cigar;

        // Convert CIGAR to tracepoints using query (A) and target (B) coordinates.
        let tracepoints = cigar_to_tracepoints(&paf_cigar, max_diff);
        let variable_tracepoints = cigar_to_variable_tracepoints(paf_cigar, max_diff);
        
        // Also convert using raw functions
        let tracepoints_raw = cigar_to_tracepoints_raw(&paf_cigar, max_diff);
        let variable_tracepoints_raw = cigar_to_variable_tracepoints_raw(&paf_cigar, max_diff);
        
        // Convert using diagonal distance functions
        let tracepoints_diagonal = cigar_to_tracepoints_diagonal(&paf_cigar, max_diff);
        let mixed_tracepoints_diagonal = cigar_to_mixed_tracepoints_diagonal(&paf_cigar, max_diff);
        let variable_tracepoints_diagonal = cigar_to_variable_tracepoints_diagonal(&paf_cigar, max_diff);

        // Compare tracepoints (allowing variable tracepoints to have None for second coordinate)
        if tracepoints
            .iter()
            .zip(variable_tracepoints.iter())
            .any(|(a, b)| {
                a.0 != b.0 || (b.1.is_some() && Some(a.1) != b.1)
            })
        {
            println!("Tracepoints mismatch! {}", line);
            println!("\t         tracepoints: {:?}", tracepoints);
            println!("\tvariable_tracepoints: {:?}", variable_tracepoints);
            println!("\t     tracepoints_raw: {:?}", tracepoints_raw);
            println!("\tvariable_tracepoints_raw: {:?}", variable_tracepoints_raw);
            println!("\t tracepoints_diagonal: {:?}", tracepoints_diagonal);
            println!("\t mixed_tracepoints_diagonal: {:?}", mixed_tracepoints_diagonal);
            println!("\t variable_tracepoints_diagonal: {:?}", variable_tracepoints_diagonal);
            std::process::exit(1);
        }

        // Reconstruct the CIGAR from tracepoints.
        let distance_mode = DistanceMode::Affine2p {
            mismatch,
            gap_open1,
            gap_ext1,
            gap_open2,
            gap_ext2,
        };
        let cigar_from_tracepoints = tracepoints_to_cigar(
            &tracepoints,
            &query_seq,
            &target_seq,
            0,
            0,
            &distance_mode,
        );
        let cigar_from_variable_tracepoints = variable_tracepoints_to_cigar(
            &variable_tracepoints,
            &query_seq,
            &target_seq,
            0,
            0,
            &distance_mode,
        );
        
        // Reconstruct CIGAR from raw tracepoints
        let cigar_from_tracepoints_raw = tracepoints_to_cigar(
            &tracepoints_raw,
            &query_seq,
            &target_seq,
            0,
            0,
            &distance_mode,
        );
        let cigar_from_variable_tracepoints_raw = variable_tracepoints_to_cigar(
            &variable_tracepoints_raw,
            &query_seq,
            &target_seq,
            0,
            0,
            &distance_mode,
        );

        // Reconstruct CIGAR from diagonal tracepoints
        let cigar_from_tracepoints_diagonal = lib_tracepoints::tracepoints_to_cigar_diagonal(
            &tracepoints_diagonal,
            &query_seq,
            &target_seq,
            0,
            0,
            &distance_mode,
        );
        let cigar_from_mixed_tracepoints_diagonal = lib_tracepoints::mixed_tracepoints_to_cigar_diagonal(
            &mixed_tracepoints_diagonal,
            &query_seq,
            &target_seq,
            0,
            0,
            &distance_mode,
        );
        let cigar_from_variable_tracepoints_diagonal = lib_tracepoints::variable_tracepoints_to_cigar_diagonal(
            &variable_tracepoints_diagonal,
            &query_seq,
            &target_seq,
            0,
            0,
            &distance_mode,
        );

        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, paf_gap_compressed_id, paf_block_id) = calculate_cigar_stats(&paf_cigar);
        let (tracepoints_matches, tracepoints_mismatches, tracepoints_insertions, tracepoints_inserted_bp, tracepoints_deletions, tracepoints_deleted_bp, tracepoints_gap_compressed_id, tracepoints_block_id) = calculate_cigar_stats(&cigar_from_tracepoints);
        let (variable_tracepoints_matches, variable_tracepoints_mismatches, variable_tracepoints_insertions, variable_tracepoints_inserted_bp, variable_tracepoints_deletions, variable_tracepoints_deleted_bp, variable_tracepoints_gap_compressed_id, variable_tracepoints_block_id) = calculate_cigar_stats(&cigar_from_variable_tracepoints);
        let (tracepoints_raw_matches, tracepoints_raw_mismatches, tracepoints_raw_insertions, tracepoints_raw_inserted_bp, tracepoints_raw_deletions, tracepoints_raw_deleted_bp, tracepoints_raw_gap_compressed_id, tracepoints_raw_block_id) = calculate_cigar_stats(&cigar_from_tracepoints_raw);
        let (variable_tracepoints_raw_matches, variable_tracepoints_raw_mismatches, variable_tracepoints_raw_insertions, variable_tracepoints_raw_inserted_bp, variable_tracepoints_raw_deletions, variable_tracepoints_raw_deleted_bp, variable_tracepoints_raw_gap_compressed_id, variable_tracepoints_raw_block_id) = calculate_cigar_stats(&cigar_from_variable_tracepoints_raw);
        // Calculate stats for diagonal CIGAR reconstructions
        let (tracepoints_diagonal_matches, tracepoints_diagonal_mismatches, tracepoints_diagonal_insertions, tracepoints_diagonal_inserted_bp, tracepoints_diagonal_deletions, tracepoints_diagonal_deleted_bp, tracepoints_diagonal_gap_compressed_id, tracepoints_diagonal_block_id) = calculate_cigar_stats(&cigar_from_tracepoints_diagonal);
        let (mixed_tracepoints_diagonal_matches, mixed_tracepoints_diagonal_mismatches, mixed_tracepoints_diagonal_insertions, mixed_tracepoints_diagonal_inserted_bp, mixed_tracepoints_diagonal_deletions, mixed_tracepoints_diagonal_deleted_bp, mixed_tracepoints_diagonal_gap_compressed_id, mixed_tracepoints_diagonal_block_id) = calculate_cigar_stats(&cigar_from_mixed_tracepoints_diagonal);
        let (variable_tracepoints_diagonal_matches, variable_tracepoints_diagonal_mismatches, variable_tracepoints_diagonal_insertions, variable_tracepoints_diagonal_inserted_bp, variable_tracepoints_diagonal_deletions, variable_tracepoints_diagonal_deleted_bp, variable_tracepoints_diagonal_gap_compressed_id, variable_tracepoints_diagonal_block_id) = calculate_cigar_stats(&cigar_from_variable_tracepoints_diagonal);
        let (realign_matches, realign_mismatches, realign_insertions, realign_inserted_bp, realign_deletions, realign_deleted_bp, realign_gap_compressed_id, realign_block_id) = calculate_cigar_stats(&realn_cigar);

        let score_from_realign = compute_alignment_score_from_cigar(&realn_cigar, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let score_from_paf = compute_alignment_score_from_cigar(&paf_cigar, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let score_from_tracepoints = compute_alignment_score_from_cigar(&cigar_from_tracepoints, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let score_from_variable_tracepoints = compute_alignment_score_from_cigar(&cigar_from_variable_tracepoints, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let score_from_tracepoints_raw = compute_alignment_score_from_cigar(&cigar_from_tracepoints_raw, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let score_from_variable_tracepoints_raw = compute_alignment_score_from_cigar(&cigar_from_variable_tracepoints_raw, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        // Calculate scores for diagonal CIGAR reconstructions
        let score_from_tracepoints_diagonal = compute_alignment_score_from_cigar(&cigar_from_tracepoints_diagonal, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let score_from_mixed_tracepoints_diagonal = compute_alignment_score_from_cigar(&cigar_from_mixed_tracepoints_diagonal, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let score_from_variable_tracepoints_diagonal = compute_alignment_score_from_cigar(&cigar_from_variable_tracepoints_diagonal, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);

        if tracepoints.len() > tracepoints_diagonal.len()
        //if variable_tracepoints.len() != variable_tracepoints_diagonal.len()
        //if cigar_from_tracepoints != cigar_from_variable_tracepoints || (paf_cigar != cigar_from_tracepoints && score_from_paf != score_from_tracepoints) //&& paf_gap_compressed_id != tracepoints_gap_compressed_id
        {
            println!("CIGAR mismatch! {}", line);
            println!("\t seqa: {}", String::from_utf8(query_seq.clone()).unwrap());
            println!("\t seqb: {}", String::from_utf8(target_seq.clone()).unwrap());
            println!("\t                      CIGAR from realign: {}", realn_cigar);
            println!("\t                          CIGAR from PAF: {}", paf_cigar);
            println!("\t                  CIGAR from tracepoints: {}", cigar_from_tracepoints);
            println!("\t              CIGAR from tracepoints_raw: {}", cigar_from_tracepoints_raw);
            //println!("\t         CIGAR from variable_tracepoints: {}", cigar_from_variable_tracepoints);
            println!("\t     CIGAR from variable_tracepoints_raw: {}", cigar_from_variable_tracepoints_raw);
            println!("\t         CIGAR from tracepoints_diagonal: {}", cigar_from_tracepoints_diagonal);
            //println!("\t   CIGAR from mixed_tracepoints_diagonal: {}", cigar_from_mixed_tracepoints_diagonal);
            //println!("\tCIGAR from variable_tracepoints_diagonal: {}", cigar_from_variable_tracepoints_diagonal);
            println!("\t                      CIGAR score from realign: {}", score_from_realign);
            println!("\t                          CIGAR score from PAF: {}", score_from_paf);
            println!("\t                  CIGAR score from tracepoints: {}", score_from_tracepoints);
            println!("\t              CIGAR score from tracepoints_raw: {}", score_from_tracepoints_raw);
            //println!("\t         CIGAR score from variable tracepoints: {}", score_from_variable_tracepoints);
            //println!("\t     CIGAR score from variable_tracepoints_raw: {}", score_from_variable_tracepoints_raw);
            println!("\t         CIGAR score from tracepoints_diagonal: {}", score_from_tracepoints_diagonal);
            //println!("\t   CIGAR score from mixed_tracepoints_diagonal: {}", score_from_mixed_tracepoints_diagonal);
            //println!("\tCIGAR score from variable_tracepoints_diagonal: {}", score_from_variable_tracepoints_diagonal);
            println!("\t              cigar stats from realign: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
                realign_matches, realign_mismatches, realign_insertions, realign_inserted_bp, realign_deletions, realign_deleted_bp, realign_gap_compressed_id, realign_block_id);
            println!("\t                           cigar stats from PAF: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
                matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, paf_gap_compressed_id, paf_block_id);
            println!("\t                   cigar stats from tracepoints: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
                tracepoints_matches, tracepoints_mismatches, tracepoints_insertions, tracepoints_inserted_bp, tracepoints_deletions, tracepoints_deleted_bp, tracepoints_gap_compressed_id, tracepoints_block_id);
            println!("\t               cigar stats from tracepoints_raw: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
                tracepoints_raw_matches, tracepoints_raw_mismatches, tracepoints_raw_insertions, tracepoints_raw_inserted_bp, tracepoints_raw_deletions, tracepoints_raw_deleted_bp, tracepoints_raw_gap_compressed_id, tracepoints_raw_block_id);
            //println!("\t          cigar stats from variable_tracepoints: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
            //    variable_tracepoints_matches, variable_tracepoints_mismatches, variable_tracepoints_insertions, variable_tracepoints_inserted_bp, variable_tracepoints_deletions, variable_tracepoints_deleted_bp, variable_tracepoints_gap_compressed_id, variable_tracepoints_block_id);
            //println!("\t      cigar stats from variable_tracepoints_raw: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
            //    variable_tracepoints_raw_matches, variable_tracepoints_raw_mismatches, variable_tracepoints_raw_insertions, variable_tracepoints_raw_inserted_bp, variable_tracepoints_raw_deletions, variable_tracepoints_raw_deleted_bp, variable_tracepoints_raw_gap_compressed_id, variable_tracepoints_raw_block_id);
            println!("\t          cigar stats from tracepoints_diagonal: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
                tracepoints_diagonal_matches, tracepoints_diagonal_mismatches, tracepoints_diagonal_insertions, tracepoints_diagonal_inserted_bp, tracepoints_diagonal_deletions, tracepoints_diagonal_deleted_bp, tracepoints_diagonal_gap_compressed_id, tracepoints_diagonal_block_id);
            //println!("\t    cigar stats from mixed_tracepoints_diagonal: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
            //    mixed_tracepoints_diagonal_matches, mixed_tracepoints_diagonal_mismatches, mixed_tracepoints_diagonal_insertions, mixed_tracepoints_diagonal_inserted_bp, mixed_tracepoints_diagonal_deletions, mixed_tracepoints_diagonal_deleted_bp, mixed_tracepoints_diagonal_gap_compressed_id, mixed_tracepoints_diagonal_block_id);
            //println!("\t cigar stats from variable_tracepoints_diagonal: matches: {}, mismatches: {}, insertions: {}, inserted_bp: {}, deletions: {}, deleted_bp: {}, gap_compressed_id: {:.12}, block_id: {:.12}",
            //    variable_tracepoints_diagonal_matches, variable_tracepoints_diagonal_mismatches, variable_tracepoints_diagonal_insertions, variable_tracepoints_diagonal_inserted_bp, variable_tracepoints_diagonal_deletions, variable_tracepoints_diagonal_deleted_bp, variable_tracepoints_diagonal_gap_compressed_id, variable_tracepoints_diagonal_block_id);
            println!("\t                   tracepoints: {:?}", tracepoints);
            //println!("\t          variable_tracepoints: {:?}", variable_tracepoints);
            println!("\t               tracepoints_raw: {:?}", tracepoints_raw);
            //println!("\t      variable_tracepoints_raw: {:?}", variable_tracepoints_raw);
            println!("\t          tracepoints_diagonal: {:?}", tracepoints_diagonal);
            //println!("\t    mixed_tracepoints_diagonal: {:?}", mixed_tracepoints_diagonal);
            //println!("\t variable_tracepoints_diagonal: {:?}", variable_tracepoints_diagonal);
            println!("\t                           bounds CIGAR from PAF: {:?}", get_cigar_diagonal_bounds(&paf_cigar));
            println!("\t                   bounds CIGAR from tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_tracepoints));
            println!("\t               bounds CIGAR from tracepoints_raw: {:?}", get_cigar_diagonal_bounds(&cigar_from_tracepoints_raw));
            println!("\t          bounds CIGAR from tracepoints_diagonal: {:?}", get_cigar_diagonal_bounds(&cigar_from_tracepoints_diagonal));
            //println!("\t    bounds CIGAR from mixed_tracepoints_diagonal: {:?}", get_cigar_diagonal_bounds(&cigar_from_mixed_tracepoints_diagonal));
            //println!("\t bounds CIGAR from variable_tracepoints_diagonal: {:?}", get_cigar_diagonal_bounds(&cigar_from_variable_tracepoints_diagonal));
            //println!("\tbounds CIGAR from variable_tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_variable_tracepoints));
            //println!("\tbounds CIGAR from variable_tracepoints_raw: {:?}", get_cigar_diagonal_bounds(&cigar_from_variable_tracepoints_raw));

            //let (deviation, d_min, d_max, max_gap) = compute_deviation(&cigar_from_tracepoints);
            //println!("\t                 deviation CIGAR from PAF: {:?}", compute_deviation(&paf_cigar));
            //println!("\t         deviation CIGAR from tracepoints: {:?}", (deviation, d_min, d_max, max_gap));
            //println!("\t     deviation CIGAR from tracepoints_raw: {:?}", compute_deviation(&cigar_from_tracepoints_raw));
            //println!("\tdeviation CIGAR from variable_tracepoints: {:?}", compute_deviation(&cigar_from_variable_tracepoints));
            //println!("\tdeviation CIGAR from variable_tracepoints_raw: {:?}", compute_deviation(&cigar_from_variable_tracepoints_raw));
            // println!("=> Try using --wfa-heuristic=banded-static --wfa-heuristic-parameters=-{},{}\n", std::cmp::max(max_gap, -d_min), std::cmp::max(max_gap, d_max));
            println!("");
        }
    });
}

/// Calculate gap compressed identity and block identity from a CIGAR string
#[cfg(debug_assertions)]
fn calculate_cigar_stats(cigar: &str) -> (usize, usize, usize, usize, usize, usize, f64, f64) {
    let mut matches = 0;
    let mut mismatches = 0;
    let mut insertions = 0; // Number of insertion events
    let mut inserted_bp = 0; // Total inserted base pairs
    let mut deletions = 0; // Number of deletion events
    let mut deleted_bp = 0; // Total deleted base pairs

    // Parse CIGAR string
    let mut num_buffer = String::new();

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buffer.push(c);
        } else {
            // Get the count
            let len = num_buffer.parse::<usize>().unwrap_or(0);
            num_buffer.clear();

            match c {
                'M' => {
                    // Assuming 'M' represents matches for simplicity (as in your code)
                    matches += len;
                }
                '=' => {
                    matches += len;
                }
                'X' => {
                    mismatches += len;
                }
                'I' => {
                    insertions += 1; // One insertion event
                    inserted_bp += len; // Total inserted bases
                }
                'D' => {
                    deletions += 1; // One deletion event
                    deleted_bp += len; // Total deleted bases
                }
                'S' | 'H' | 'P' | 'N' => {
                    // Skip soft clips, hard clips, padding, and skipped regions
                }
                _ => {
                    // Unknown operation, skip
                }
            }
        }
    }

    // Calculate gap compressed identity
    let gap_compressed_identity = if matches + mismatches + insertions + deletions > 0 {
        (matches as f64) / (matches + mismatches + insertions + deletions) as f64
    } else {
        0.0
    };

    // Calculate block identity
    let edit_distance = mismatches + inserted_bp + deleted_bp;
    let block_identity = if matches + edit_distance > 0 {
        (matches as f64) / (matches + edit_distance) as f64
    } else {
        0.0
    };

    (
        matches,
        mismatches,
        insertions,
        inserted_bp,
        deletions,
        deleted_bp,
        gap_compressed_identity,
        block_identity,
    )
}

#[cfg(debug_assertions)]
fn compute_alignment_score_from_cigar(
    cigar: &str,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
) -> i32 {
    let mut score = 0i32;
    let mut num_buffer = String::new();

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buffer.push(c);
        } else {
            // Get the count
            let len = num_buffer.parse::<i32>().unwrap_or(0);
            num_buffer.clear();

            match c {
                '=' => {
                    // Matches - no penalty (score 0)
                    // score += 0;
                }
                'M' => {
                    // For 'M' operations, we'd need the actual sequences to determine matches vs mismatches
                    // For now, we'll treat 'M' as matches (you may want to adjust this)
                    // score += 0;
                    eprintln!(
                        "Warning: 'M' in CIGAR requires sequences to determine match/mismatch"
                    );
                }
                'X' => {
                    // Mismatches
                    score -= mismatch * len;
                }
                'I' | 'D' => {
                    // Gaps - using dual affine model
                    // Calculate both penalty options and take the minimum (best score)
                    let score1 = gap_open1 + gap_ext1 * len;
                    let score2 = gap_open2 + gap_ext2 * len;
                    let gap_penalty = std::cmp::min(score1, score2);
                    score -= gap_penalty;
                }
                'S' | 'H' | 'P' | 'N' => {
                    // Soft clips, hard clips, padding, and skipped regions
                    // These typically don't contribute to the alignment score
                }
                _ => {
                    eprintln!("Unknown CIGAR operation: {}", c);
                }
            }
        }
    }

    score
}
/// Initialize logger based on verbosity
fn setup_logger(verbosity: u8) {
    env_logger::Builder::new()
        .filter_level(match verbosity {
            0 => log::LevelFilter::Warn,  // Errors and warnings
            1 => log::LevelFilter::Info,  // Errors, warnings, and info
            _ => log::LevelFilter::Debug, // Errors, warnings, info, and debug
        })
        .init();
}

fn get_paf_reader(paf: &str) -> io::Result<Box<dyn BufRead>> {
    if paf == "-" {
        Ok(Box::new(BufReader::new(std::io::stdin())))
    } else if paf.ends_with(".gz") || paf.ends_with(".bgz") {
        let file = File::open(paf)?;
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        let file = File::open(paf)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Process a chunk of lines in parallel for compression
fn process_compress_chunk(lines: &[String], tp_type: &TracepointType, max_diff: usize) {
    lines.par_iter().for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            error!(
                "{}",
                message_with_truncate_paf_file("Skipping malformed PAF line", line)
            );
            std::process::exit(1);
        }

        let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
            error!(
                "{}",
                message_with_truncate_paf_file("Skipping CIGAR-less PAF line", line)
            );
            std::process::exit(1);
        };
        let cigar = &cg_field[5..];

        // Calculate identity stats from the CIGAR
        let (gap_compressed_identity, block_identity) = calculate_identity_stats(cigar);

        // Calculate alignment score based on edit distance
        let alignment_score = calculate_alignment_score_edit_distance(cigar);

        // Check for existing df:i: field
        let existing_df = fields.iter()
            .find(|&&s| s.starts_with("df:i:"))
            .map(|&s| s[5..].parse::<usize>().ok())
            .flatten();

        // Handle FastGA tracepoints with potential overflow splitting
        if matches!(tp_type, TracepointType::Fastga) {
            process_fastga_with_overflow(
                &fields,
                cigar,
                max_diff,
                gap_compressed_identity,
                block_identity,
                alignment_score,
                existing_df,
                cg_field,
            );
        } else {
            // Handle other tracepoint types (no overflow splitting needed)
            process_single_record(
                &fields,
                cigar,
                tp_type,
                max_diff,
                gap_compressed_identity,
                block_identity,
                alignment_score,
                existing_df,
                cg_field,
            );
        }
    });
}

fn process_fastga_with_overflow(
    fields: &[&str],
    cigar: &str,
    max_diff: usize,
    gap_compressed_identity: f64,
    block_identity: f64,
    alignment_score: i32,
    existing_df: Option<usize>,
    cg_field: &str,
) {
    // Parse coordinates
    let query_len = fields[1].parse().unwrap_or_else(|_| {
        error!("Invalid query_len in PAF line");
        std::process::exit(1);
    });
    let query_start: usize = fields[2].parse().unwrap_or_else(|_| {
        error!("Invalid query_start in PAF line");
        std::process::exit(1);
    });
    let query_end: usize = fields[3].parse().unwrap_or_else(|_| {
        error!("Invalid query_end in PAF line");
        std::process::exit(1);
    });
    let target_len: usize = fields[6].parse().unwrap_or_else(|_| {
        error!("Invalid target_length in PAF line");
        std::process::exit(1);
    });
    let target_start: usize = fields[7].parse().unwrap_or_else(|_| {
        error!("Invalid target_start in PAF line");
        std::process::exit(1);
    });
    let target_end: usize = fields[8].parse().unwrap_or_else(|_| {
        error!("Invalid target_end in PAF line");
        std::process::exit(1);
    });
    let complement = fields[4] == "-";

    // Process with overflow handling
    let segments = cigar_to_tracepoints_fastga(
        cigar, 
        max_diff, 
        query_start, 
        query_end, 
        query_len, 
        target_start, 
        target_end, 
        target_len, 
        complement
    );

    // Output each segment as a separate PAF record
    for (segment_idx, (tracepoints, (seg_query_start, seg_query_end, seg_target_start, seg_target_end))) in segments.iter().enumerate() {
        // Calculate segment-specific stats
        let sum_of_differences: usize = tracepoints.iter().map(|(diff, _)| diff).sum();
        let tracepoints_str = format_tracepoints(tracepoints);
        
        // Create modified fields for this segment
        let mut new_fields: Vec<String> = Vec::new();
        
        for (i, field) in fields.iter().enumerate() {
            match i {
                2 => new_fields.push(seg_query_start.to_string()), // query_start
                3 => new_fields.push(seg_query_end.to_string()),   // query_end
                7 => {
                    // Handle target coordinates based on strand
                    if complement {
                        // For reverse strand, transform back to PAF coordinates
                        new_fields.push((target_len - seg_target_end).to_string());
                    } else {
                        new_fields.push(seg_target_start.to_string());
                    }
                }
                8 => {
                    // Handle target coordinates based on strand
                    if complement {
                        // For reverse strand, transform back to PAF coordinates
                        new_fields.push((target_len - seg_target_start).to_string());
                    } else {
                        new_fields.push(seg_target_end.to_string());
                    }
                }
                9 => {
                    // Update residue_matches for this segment
                    let alignment_length = seg_query_end - seg_query_start;
                    let matches = alignment_length.saturating_sub(sum_of_differences);
                    new_fields.push(matches.to_string());
                }
                10 => {
                    // Update alignment_block_length for this segment
                    let alignment_length = seg_query_end - seg_query_start;
                    new_fields.push(alignment_length.to_string());
                }
                _ => {
                    if field.starts_with("cg:Z:") {
                        // Add identity stats before tracepoints
                        new_fields.push(format!("gi:f:{:.12}", gap_compressed_identity));
                        new_fields.push(format!("bi:f:{:.12}", block_identity));

                        // Add df fields
                        if let Some(old_df) = existing_df {
                            new_fields.push(format!("dfold:i:{}", old_df));
                        }
                        new_fields.push(format!("df:i:{}", sum_of_differences));

                        // Add alignment score
                        new_fields.push(format!("sc:i:{}", alignment_score));

                        // Add tracepoints
                        new_fields.push(format!("tp:Z:{}", tracepoints_str));
                    } else if field.starts_with("df:i:") || field.starts_with("gi:f:") || 
                              field.starts_with("bi:f:") || field.starts_with("sc:i:") {
                        // Skip existing fields as we handle them above
                        continue;
                    } else {
                        new_fields.push(field.to_string());
                    }
                }
            }
        }

        // Add segment marker if this record was split
        if segments.len() > 1 {
            new_fields.push(format!("sg:i:{}", segment_idx));
        }

        println!("{}", new_fields.join("\t"));
    }
}

fn process_single_record(
    fields: &[&str],
    cigar: &str,
    tp_type: &TracepointType,
    max_diff: usize,
    gap_compressed_identity: f64,
    block_identity: f64,
    alignment_score: i32,
    existing_df: Option<usize>,
    cg_field: &str,
) {
    // Convert CIGAR based on tracepoint type
    let (tracepoints_str, df_value) = match tp_type {
        TracepointType::Mixed => {
            let tp = cigar_to_mixed_tracepoints(cigar, max_diff);
            (format_mixed_tracepoints(&tp), None::<usize>)
        }
        TracepointType::Variable => {
            let tp = cigar_to_variable_tracepoints(cigar, max_diff);
            (format_variable_tracepoints(&tp), None)
        }
        TracepointType::Standard => {
            let tp = cigar_to_tracepoints(cigar, max_diff);
            (format_tracepoints(&tp), None)
        }
        TracepointType::Fastga => {
            // This should not happen as FastGA is handled separately
            unreachable!("FastGA should be handled by process_fastga_with_overflow")
        }
    };

    // Build the new line
    let mut new_fields: Vec<String> = Vec::new();

    for field in fields.iter() {
        if field.starts_with("cg:Z:") {
            // Add identity stats before tracepoints
            new_fields.push(format!("gi:f:{:.12}", gap_compressed_identity));
            new_fields.push(format!("bi:f:{:.12}", block_identity));

            // Add df field if needed
            if let Some(new_df) = df_value {
                if let Some(old_df) = existing_df {
                    new_fields.push(format!("dfold:i:{}", old_df));
                }
                new_fields.push(format!("df:i:{}", new_df));
            }

            // Add alignment score
            new_fields.push(format!("sc:i:{}", alignment_score));

            // Add tracepoints
            new_fields.push(format!("tp:Z:{}", tracepoints_str));
        } else if field.starts_with("df:i:") || field.starts_with("gi:f:") || 
                  field.starts_with("bi:f:") || field.starts_with("sc:i:") {
            // Skip existing fields as we handle them above
            continue;
        } else {
            new_fields.push(field.to_string());
        }
    }
    
    println!("{}", new_fields.join("\t"));
}

/// Process a chunk of lines in parallel for decompression
fn process_decompress_chunk(
    lines: &[String],
    tp_type: &TracepointType,
    query_fasta_path: &str,
    target_fasta_path: &str,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
    fastga: bool,
    trace_spacing: usize,
) {
    lines.par_iter().for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            error!(
                "{}",
                message_with_truncate_paf_file("Skipping malformed PAF line", line)
            );
            std::process::exit(1);
        }

        let Some(tp_field) = fields.iter().find(|&&s| s.starts_with("tp:Z:")) else {
            error!(
                "{}",
                message_with_truncate_paf_file("Skipping tracepoints-less PAF line", line)
            );
            std::process::exit(1);
        };
        let tracepoints_str = &tp_field[5..]; // Direct slice instead of strip_prefix("tp:Z:")

        // Parse mandatory PAF fields
        let query_name = fields[0];
        let query_start: usize = fields[2].parse().unwrap_or_else(|_| {
            error!(
                "{}",
                message_with_truncate_paf_file("Invalid query_start in PAF line", line)
            );
            std::process::exit(1);
        });
        let query_end: usize = fields[3].parse().unwrap_or_else(|_| {
            error!(
                "{}",
                message_with_truncate_paf_file("Invalid query_end in PAF line", line)
            );
            std::process::exit(1);
        });
        let strand = fields[4];
        let target_name = fields[5];
        let target_len: usize = fields[6].parse().unwrap_or_else(|_| {
            error!(
                "{}",
                message_with_truncate_paf_file("Invalid target_length in PAF line", line)
            );
            std::process::exit(1);
        });
        let target_start: usize = fields[7].parse().unwrap_or_else(|_| {
            error!(
                "{}",
                message_with_truncate_paf_file("Invalid target_start in PAF line", line)
            );
            std::process::exit(1);
        });
        let target_end: usize = fields[8].parse().unwrap_or_else(|_| {
            error!(
                "{}",
                message_with_truncate_paf_file("Invalid target_end in PAF line", line)
            );
            std::process::exit(1);
        });

        // Create thread-local FASTA readers for query and target
        let query_fasta_reader = FastaReader::from_path(query_fasta_path).unwrap_or_else(|e| {
            error!("Failed to create query FASTA reader: {}", e);
            std::process::exit(1);
        });

        let target_fasta_reader = FastaReader::from_path(target_fasta_path).unwrap_or_else(|e| {
            error!("Failed to create target FASTA reader: {}", e);
            std::process::exit(1);
        });

        // Fetch query sequence from query FASTA (with FASTGA we always use forward strand for the query)
        let query_seq = if strand == "+" || fastga {
            debug!("Fetching query sequence {}:{}-{} on + strand", query_name, query_start, query_end);

            match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut seq_vec = seq.to_vec();
                    unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) }; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                    seq_vec
                        .iter_mut()
                        .for_each(|byte| *byte = byte.to_ascii_uppercase());
                    seq_vec
                }
                Err(e) => {
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
                }
            }
        } else {
            debug!("Fetching query sequence {}:{}-{} on - strand", query_name, query_start, query_end);

            match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut rc = reverse_complement(&seq.to_vec());
                    unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) }; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                    rc.iter_mut()
                        .for_each(|byte| *byte = byte.to_ascii_uppercase());
                    rc
                }
                Err(e) => {
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
                }
            }
        };

        // Fetch target sequence from target FASTA (with FASTGA the strand affects the target sequence)
        let target_seq = if strand == "+" || !fastga {
            debug!("Fetching target sequence {}:{}-{} on + strand", target_name, target_start, target_end);

            match target_fasta_reader.fetch_seq(target_name, target_start, target_end - 1) {
                Ok(seq) => {
                    let mut seq_vec = seq.to_vec();
                    unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) }; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                    seq_vec
                        .iter_mut()
                        .for_each(|byte| *byte = byte.to_ascii_uppercase());
                    seq_vec
                }
                Err(e) => {
                    error!("Failed to fetch target sequence: {}", e);
                    std::process::exit(1);
                }
            }
        } else {
            //let (target_start, target_end) = (target_len - target_end, target_len - target_start);
            debug!("Fetching target sequence {}:{}-{} on - strand", target_name, target_start, target_end);

            match target_fasta_reader.fetch_seq(target_name, target_start, target_end - 1) {
                Ok(seq) => {
                    let mut rc = reverse_complement(&seq.to_vec());
                    unsafe { libc::free(seq.as_ptr() as *mut std::ffi::c_void) }; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                    rc.iter_mut()
                        .for_each(|byte| *byte = byte.to_ascii_uppercase());
                    rc
                }
                Err(e) => {
                    error!("Failed to fetch target sequence: {}", e);
                    std::process::exit(1);
                }
            }
        };

        // Use specified tracepoint type
        let distance_mode = DistanceMode::Affine2p {
            mismatch,
            gap_open1,
            gap_ext1,
            gap_open2,
            gap_ext2,
        };
        let cigar = if fastga {
            // FastGA decoding
            let tracepoints = parse_tracepoints(tracepoints_str);
            let query_start: usize = fields[2].parse().unwrap_or_else(|_| {
                error!(
                    "{}",
                    message_with_truncate_paf_file("Invalid query_start in PAF line", line)
                );
                std::process::exit(1);
            });
            let target_start: usize = fields[7].parse().unwrap_or_else(|_| {
                error!(
                    "{}",
                    message_with_truncate_paf_file("Invalid target_start in PAF line", line)
                );
                std::process::exit(1);
            });
            tracepoints_to_cigar_fastga(
                &tracepoints,
                trace_spacing,
                &query_seq,
                &target_seq,
                query_start,
                if  strand == "+" { target_len - target_end } else { target_start },
                strand == "-",
            )
        } else {
            match tp_type {
                TracepointType::Mixed => {
                    // Mixed representation
                    let mixed_tracepoints = parse_mixed_tracepoints(tracepoints_str);
                    mixed_tracepoints_to_cigar(
                        &mixed_tracepoints,
                        &query_seq,
                        &target_seq,
                        0,
                        0,
                        &distance_mode,
                    )
                }
                TracepointType::Variable => {
                    // Variable tracepoints representation
                    let variable_tracepoints = parse_variable_tracepoints(tracepoints_str);
                    variable_tracepoints_to_cigar(
                        &variable_tracepoints,
                        &query_seq,
                        &target_seq,
                        0,
                        0,
                        &distance_mode,
                    )
                }
                TracepointType::Standard => {
                    // Standard tracepoints
                    let tracepoints = parse_tracepoints(tracepoints_str);
                    tracepoints_to_cigar(
                        &tracepoints,
                        &query_seq,
                        &target_seq,
                        0,
                        0,
                        &distance_mode,
                    )
                }
                TracepointType::Fastga => {
                    // Should not reach here - fastga is handled above
                    error!("Fastga should not be handled in non-fastga path");
                    std::process::exit(1);
                }
            }
        };

        // Calculate identity stats from the reconstructed CIGAR
        let (gap_compressed_identity, block_identity) = calculate_identity_stats(&cigar);

        // Calculate alignment score based on edit distance
        let alignment_score = calculate_alignment_score_edit_distance(&cigar);

        // Check for existing gi:f:, bi:f:, and sc:i: fields
        let existing_gi = fields.iter().find(|&&s| s.starts_with("gi:f:"));
        let existing_bi = fields.iter().find(|&&s| s.starts_with("bi:f:"));
        let existing_sc = fields.iter().find(|&&s| s.starts_with("sc:i:"));

        // Build the new line
        let mut new_fields: Vec<String> = Vec::new();

        for field in fields.iter() {
            if field.starts_with("tp:Z:") {
                // If there were existing gi/bi/sc fields, rename them with 'old' prefix
                if let Some(old_gi) = existing_gi {
                    new_fields.push(format!("giold:f:{}", &old_gi[5..]));
                }
                if let Some(old_bi) = existing_bi {
                    new_fields.push(format!("biold:f:{}", &old_bi[5..]));
                }
                if let Some(old_sc) = existing_sc {
                    new_fields.push(format!("scold:i:{}", &old_sc[5..]));
                }

                // Add new identity stats and alignment score before the CIGAR
                new_fields.push(format!("gi:f:{:.12}", gap_compressed_identity));
                new_fields.push(format!("bi:f:{:.12}", block_identity));
                new_fields.push(format!("sc:i:{}", alignment_score));

                // Replace tracepoints with CIGAR
                new_fields.push(format!("cg:Z:{}", cigar));
            } else if field.starts_with("gi:f:") || field.starts_with("bi:f:") || field.starts_with("sc:i:") {
                // Skip existing gi, bi, and sc fields - we've already handled them
                continue;
            } else {
                new_fields.push(field.to_string());
            }
        }
        
        println!("{}", new_fields.join("\t"));
    });
}

/// Combines a message with the first 9 columns of a PAF line.
fn message_with_truncate_paf_file(message: &str, line: &str) -> String {
    let truncated_line = line.split('\t').take(9).collect::<Vec<&str>>().join("\t");
    format!("{}: {} ...", message, truncated_line)
}

fn format_tracepoints(tracepoints: &[(usize, usize)]) -> String {
    tracepoints
        .iter()
        .map(|(a, b)| format!("{},{}", a, b))
        .collect::<Vec<String>>()
        .join(";")
}
fn format_mixed_tracepoints(mixed_tracepoints: &[MixedRepresentation]) -> String {
    mixed_tracepoints
        .iter()
        .map(|tp| match tp {
            MixedRepresentation::Tracepoint(a, b) => format!("{},{}", a, b),
            MixedRepresentation::CigarOp(len, op) => format!("{}{}", len, op),
        })
        .collect::<Vec<String>>()
        .join(";")
}

fn format_variable_tracepoints(variable_tracepoints: &[(usize, Option<usize>)]) -> String {
    variable_tracepoints
        .iter()
        .map(|(a, b_opt)| match b_opt {
            Some(b) => format!("{},{}", a, b),
            None => format!("{}", a),
        })
        .collect::<Vec<String>>()
        .join(";")
}
fn parse_penalties(
    penalties: &str,
) -> Result<(i32, i32, i32, i32, i32), Box<dyn std::error::Error>> {
    let tokens: Vec<&str> = penalties.split(',').collect();
    if tokens.len() != 5 {
        error!(
            "Error: penalties must be provided as mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2"
        );
        std::process::exit(1);
    }

    let mismatch: i32 = tokens[0].parse()?;
    let gap_open1: i32 = tokens[1].parse()?;
    let gap_ext1: i32 = tokens[2].parse()?;
    let gap_open2: i32 = tokens[3].parse()?;
    let gap_ext2: i32 = tokens[4].parse()?;

    Ok((mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2))
}

fn parse_tracepoints(tp_str: &str) -> Vec<(usize, usize)> {
    tp_str
        .split(';')
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            Some((parts[0].parse().unwrap(), parts[1].parse().unwrap()))
        })
        .collect()
}

fn parse_mixed_tracepoints(tp_str: &str) -> Vec<MixedRepresentation> {
    tp_str
        .split(';')
        .filter_map(|s| {
            if s.contains(',') {
                // This is a tracepoint
                let parts: Vec<&str> = s.split(',').collect();
                Some(MixedRepresentation::Tracepoint(
                    parts[0].parse().unwrap(),
                    parts[1].parse().unwrap(),
                ))
            } else {
                // This is a cigar operation
                let mut chars = s.chars();
                let mut len_str = String::new();

                // Read digits
                while let Some(c) = chars.next() {
                    if c.is_digit(10) {
                        len_str.push(c);
                    } else {
                        // Found operator character
                        let len = len_str.parse().unwrap();
                        return Some(MixedRepresentation::CigarOp(len, c));
                    }
                }

                // If we get here, parsing failed
                None
            }
        })
        .collect()
}

fn parse_variable_tracepoints(tp_str: &str) -> Vec<(usize, Option<usize>)> {
    tp_str
        .split(';')
        .filter_map(|s| {
            if s.contains(',') {
                // This has both coordinates
                let parts: Vec<&str> = s.split(',').collect();
                Some((parts[0].parse().unwrap(), Some(parts[1].parse().unwrap())))
            } else {
                // This has only first coordinate
                match s.parse() {
                    Ok(a) => Some((a, None)),
                    Err(_) => None,
                }
            }
        })
        .collect()
}

#[cfg(debug_assertions)]
fn get_cigar_diagonal_bounds(cigar: &str) -> (i64, i64) {
    let mut current_diagonal = 0; // Current diagonal position
    let mut min_diagonal = 0; // Lowest diagonal reached
    let mut max_diagonal = 0; // Highest diagonal reached

    // Parse CIGAR string with numerical counts
    let mut num_buffer = String::new();

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buffer.push(c);
        } else {
            // Get the count
            let count = num_buffer.parse::<i64>().unwrap();
            num_buffer.clear();

            match c {
                'M' | '=' | 'X' => {
                    // Matches stay on same diagonal
                }
                'D' => {
                    // Deletions move down diagonal by count amount
                    current_diagonal -= count;
                    min_diagonal = min_diagonal.min(current_diagonal);
                }
                'I' => {
                    // Insertions move up diagonal by count amount
                    current_diagonal += count;
                    max_diagonal = max_diagonal.max(current_diagonal);
                }
                _ => panic!("Invalid CIGAR operation: {}", c),
            }
        }
    }

    (min_diagonal, max_diagonal)
}

#[cfg(debug_assertions)]
fn compute_deviation(cigar: &str) -> (i64, i64, i64, i64) {
    let mut deviation = 0;
    let mut d_max = -10000;
    let mut d_min = 10000;
    let mut max_gap = 0;

    // Parse CIGAR string with numerical counts
    let mut num_buffer = String::new();

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buffer.push(c);
        } else {
            // Get the count
            let count = num_buffer.parse::<i64>().unwrap();
            num_buffer.clear();

            match c {
                'M' | '=' | 'X' => {
                    // Matches stay on same diagonal
                }
                'D' => {
                    // Deletions move down diagonal by count amount
                    deviation -= count;
                    max_gap = std::cmp::max(max_gap, count);
                }
                'I' => {
                    deviation += count;
                    max_gap = std::cmp::max(max_gap, count);
                }
                _ => panic!("Invalid CIGAR operation: {}", c),
            }

            d_max = std::cmp::max(d_max, deviation);
            d_min = std::cmp::min(d_min, deviation);
        }
    }

    (deviation, d_min, d_max, max_gap)
}

/// Returns the reverse complement of a DNA sequence
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            _ => b'N', // Convert any unexpected bases to N
        })
        .collect()
}

// /// Calculate alignment coordinates from a CIGAR string and starting positions
// /// Returns (query_end, query_len, target_end, target_len)
// fn calculate_alignment_coordinates(
//     cigar: &str,
//     query_start: usize,
//     target_start: usize,
// ) -> (usize, usize, usize, usize) {
//     let ops = cigar_str_to_cigar_ops(cigar);

//     let mut query_len = 0;
//     let mut target_len = 0;

//     // Calculate total lengths by checking which operations consume query/target bases
//     for &(len, op) in &ops {
//         if consumes_a(op) {
//             query_len += len;
//         }
//         if consumes_b(op) {
//             target_len += len;
//         }
//     }

//     // Calculate end positions by adding consumed lengths to start positions
//     let query_end = query_start + query_len;
//     let target_end = target_start + target_len;

//     (query_end, query_len, target_end, target_len)
// }