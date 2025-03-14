use clap::Parser;
use std::io::{self, BufRead, BufReader};
use std::fs::File;
use flate2::read::MultiGzDecoder;
use rust_htslib::faidx::Reader as FastaReader;
use lib_tracepoints::{cigar_to_tracepoints, cigar_to_single_band_tracepoints, cigar_to_double_band_tracepoints, cigar_to_mixed_tracepoints, MixedRepresentation, tracepoints_to_cigar, single_band_tracepoints_to_cigar, double_band_tracepoints_to_cigar, mixed_tracepoints_to_cigar, align_sequences_wfa, cigar_ops_to_cigar_string};
use lib_wfa2::affine_wavefront::{AffineWavefronts};
use log::{info, error};
use rayon::prelude::*;

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

        /// Banding strategy: none, single, or double
        #[arg(short = 'b', long = "band", default_value = "none")]
        band: String,
        
        /// Use mixed representation (preserves S/H/P/N CIGAR operations)
        #[arg(short = 'm', long = "mixed", default_value_t = false)]
        mixed: bool,

        /// Max-diff value for tracepoints
        #[arg(long, default_value = "32")]
        max_diff: usize,
    },
    /// Decompression of alignments
    Decompress {
        #[clap(flatten)]
        common: CommonOpts,
        
        /// FASTA file for query sequences
        #[arg(short = 'q', long = "query-fasta")]
        query_fasta: String,

        /// FASTA file for target sequences
        #[arg(short = 't', long = "target-fasta")]
        target_fasta: String,
        
        /// Gap penalties in the format mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2
        #[arg(long, default_value = "3,4,2,24,1")]
        penalties: String,
    },
    /// Run debugging mode (only available in debug builds)
    #[cfg(debug_assertions)]
    Debug {
        /// PAF file for alignments (use "-" to read from standard input)
        #[arg(short = 'p', long = "paf")]
        paf: Option<String>,

        /// FASTA file for query sequences
        #[arg(short = 'q', long = "query-fasta")]
        query_fasta: Option<String>,

        /// FASTA file for target sequences
        #[arg(short = 't', long = "target-fasta")]
        target_fasta: Option<String>,

        /// Gap penalties in the format mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2
        #[arg(long, default_value = "3,4,2,24,1")]
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
        Args::Compress { common, band, mixed, max_diff } => {
            setup_logger(common.verbose);

            // Validate and convert band type
            let band_type : u8 = match band.to_lowercase().as_str() {
                "none" => 0,
                "single" => 1,
                "double" => 2,
                _ => {
                    error!("Invalid banding strategy '{}'. Supported options: none, single, double, mixed", band);
                    std::process::exit(1);
                }
            };
            
            // If mixed is true, validate that band is "none"
            if mixed && band_type != 0 {
                error!("The --mixed option cannot be used with a banding strategy. Use --band none with --mixed.");
                std::process::exit(1);
            }

            info!("Converting CIGAR to {} tracepoints", 
                match band_type {
                    0 => "no-band",
                    1 => "single-band",
                    2 => "double-band",
                    _ => unreachable!()
                }
            );

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
                            process_compress_chunk(&lines, band_type, mixed, max_diff);
                            lines.clear();
                        }
                    },
                    Err(e) => return Err(e.into()),
                }
            }
            
            // Process remaining lines
            if !lines.is_empty() {
                process_compress_chunk(&lines, band_type, mixed, max_diff);
            }
        },
        Args::Decompress { common, query_fasta, target_fasta, penalties } => {
            setup_logger(common.verbose);
            info!("Converting tracepoints to CIGAR");

            // Set the thread pool size
            rayon::ThreadPoolBuilder::new()
                .num_threads(common.threads)
                .build_global()?;

            // Parse penalties
            let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = parse_penalties(&penalties)?;

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
                                &query_fasta,
                                &target_fasta,
                                mismatch, 
                                gap_open1, 
                                gap_ext1, 
                                gap_open2, 
                                gap_ext2
                            );
                            lines.clear();
                        }
                    },
                    Err(e) => return Err(e.into()),
                }
            }
            
            // Process remaining lines
            if !lines.is_empty() {
                process_decompress_chunk(
                    &lines, 
                    &query_fasta,
                    &target_fasta,
                    mismatch, 
                    gap_open1, 
                    gap_ext1, 
                    gap_open2, 
                    gap_ext2
                );
            }
        },
        #[cfg(debug_assertions)]
        Args::Debug { paf, query_fasta, target_fasta, penalties, max_diff, verbose } => {
            setup_logger(verbose);
            info!("Debugging");

            let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = parse_penalties(&penalties)?;
            info!("Penalties: {},{},{},{},{}", mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);

            if let (Some(paf), Some(query_fasta), Some(target_fasta)) = (paf, query_fasta, target_fasta) {
                info!("PAF file: {}", paf);
                info!("Query FASTA file: {}", query_fasta);
                info!("Target FASTA file: {}", target_fasta);

                // Open the FASTA files
                let query_fasta_reader = FastaReader::from_path(&query_fasta).expect("Error reading query FASTA file");
                let target_fasta_reader = FastaReader::from_path(&target_fasta).expect("Error reading target FASTA file");

                // Open the PAF file (or use stdin if "-" is provided).
                let paf_reader = get_paf_reader(&paf)?;

                for line in paf_reader.lines() {
                    let line = line?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 12 {
                        error!("{}", message_with_truncate_paf_file("Skipping malformed PAF line", &line));
                        std::process::exit(1);
                    }
                    
                    let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
                        error!("{}", message_with_truncate_paf_file("Skipping CIGAR-less PAF line", &line));
                        std::process::exit(1);
                    };
                    let _paf_cigar = &cg_field[5..]; // Direct slice instead of strip_prefix("cg:Z:")

                    // Parse mandatory PAF fields.
                    let query_name = fields[0];
                    //let query_len: usize = fields[1].parse()?;
                    let query_start: usize = fields[2].parse()?;
                    let query_end: usize = fields[3].parse()?;
                    let strand = fields[4];
                    let target_name = fields[5];
                    //let target_len: usize = fields[6].parse()?;
                    let target_start: usize = fields[7].parse()?;
                    let target_end: usize = fields[8].parse()?;
                    
                    // Fetch query sequence from query FASTA
                    let query_seq = if strand == "+" {
                        match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                            Ok(seq) => {
                                let mut seq_vec = seq.to_vec();
                                unsafe {libc::free(seq.as_ptr() as *mut std::ffi::c_void)}; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                                seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                                seq_vec
                            },
                            Err(e) => {
                                error!("Failed to fetch query sequence: {}", e);
                                std::process::exit(1);
                            }
                        }
                    } else {
                        match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                            Ok(seq) => {
                                let mut rc = reverse_complement(&seq.to_vec());
                                unsafe {libc::free(seq.as_ptr() as *mut std::ffi::c_void)}; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                                rc.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                                rc
                            },
                            Err(e) => {
                                error!("Failed to fetch query sequence: {}", e);
                                std::process::exit(1);
                            }
                        }
                    };

                    // Fetch target sequence from target FASTA
                    let target_seq = match target_fasta_reader.fetch_seq(target_name, target_start, target_end - 1) {
                        Ok(seq) => {
                            let mut seq_vec = seq.to_vec();
                            unsafe {libc::free(seq.as_ptr() as *mut std::ffi::c_void)}; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                            seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                            seq_vec
                        },
                        Err(e) => {
                            error!("Failed to fetch target sequence: {}", e);
                            std::process::exit(1);
                        }
                    };

                    let mut aligner = AffineWavefronts::with_penalties_affine2p(0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
                    let realn_cigar = align_sequences_wfa(&query_seq, &target_seq, &mut aligner);
                    let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);
                    let paf_cigar = &realn_cigar;

                    // Convert CIGAR to tracepoints using query (A) and target (B) coordinates.
                    let tracepoints = cigar_to_tracepoints(paf_cigar, max_diff);
                    let single_band_tracepoints = cigar_to_single_band_tracepoints(paf_cigar, max_diff);
                    let double_band_tracepoints = cigar_to_double_band_tracepoints(paf_cigar, max_diff);

                    // Compare tracepoints and banded_tracepoints, but only the first 2 elements of the tuple. if different, print error
                    if tracepoints.iter().zip(single_band_tracepoints.iter()).any(|(a, b)| a.0 != b.0 || a.1 != b.1) {
                        error!("Tracepoints mismatch! {}", line);
                        error!("\t            tracepoints: {:?}", tracepoints);
                        error!("\tsingle_band_tracepoints: {:?}", single_band_tracepoints);
                        std::process::exit(1);
                    }
                    if tracepoints.iter().zip(double_band_tracepoints.iter()).any(|(a, b)| a.0 != b.0 || a.1 != b.1) {
                        error!("Tracepoints mismatch! {}", line);
                        error!("\t            tracepoints: {:?}", tracepoints);
                        error!("\tdouble_band_tracepoints: {:?}", double_band_tracepoints);
                        std::process::exit(1);
                    }

                    // Reconstruct the CIGAR from tracepoints.
                    let cigar_from_tracepoints = tracepoints_to_cigar(
                        &tracepoints,
                        &query_seq,
                        &target_seq,
                        0,
                        0,
                        mismatch,
                        gap_open1,
                        gap_ext1,
                        gap_open2,
                        gap_ext2
                    );
                    let cigar_from_single_band_tracepoints = single_band_tracepoints_to_cigar(
                        &single_band_tracepoints,
                        &query_seq,
                        &target_seq,
                        0,
                        0,
                        mismatch,
                        gap_open1,
                        gap_ext1,
                        gap_open2,
                        gap_ext2
                    );
                    let cigar_from_double_band_tracepoints = double_band_tracepoints_to_cigar(
                        &double_band_tracepoints,
                        &query_seq,
                        &target_seq,
                        0,
                        0,
                        mismatch,
                        gap_open1,
                        gap_ext1,
                        gap_open2,
                        gap_ext2
                    );

                    // let (query_end_variable, query_len_variable, target_end_variable, target_len_variable) = calculate_alignment_coordinates(&cigar_from_tracepoints, query_start, target_start);
                    // let (query_end_variable2, query_len_variable2, target_end_variable2, target_len_variable2) = calculate_alignment_coordinates(&cigar_from_double_band_tracepoints, query_start, target_start);
                    // let (query_end_paf, query_len_paf, target_end_paf, target_len_paf) = calculate_alignment_coordinates(paf_cigar, query_start, target_start);
                    // if (query_len_paf != query_len_variable) || (target_len_paf != target_len_variable) || (query_end_paf != query_end_variable) || (target_end_paf != target_end_variable) {
                    //     info!("cigar_from_tracepoints {:?}", (query_end_variable, query_len_variable, target_end_variable, target_len_variable));
                    //     info!("           paf_cigar {:?}", (query_end_paf, query_len_paf, target_end_paf, target_len_paf) );
                    //     error!("Line {}: seq. len. mismatch!", line);
                    // }
                    // if (query_len_paf != query_len_variable2) || (target_len_paf != target_len_variable2) || (query_end_paf != query_end_variable2) || (target_end_paf != target_end_variable2) {
                    //     info!("cigar_from_banded_tracepoints {:?}", (query_end_variable2, query_len_variable2, target_end_variable2, target_len_variable2));
                    //     info!("           p af_cigar {:?}", (query_end_paf, query_len_paf, target_end_paf, target_len_paf) );
                    //     error!("Line {}: seq. len. mismatch!", line);
                    // }

                    if cigar_from_tracepoints != cigar_from_single_band_tracepoints || cigar_from_tracepoints != cigar_from_double_band_tracepoints {
                        error!("CIGAR mismatch! {}", line);
                        error!("\t                        tracepoints: {:?}", tracepoints);
                        error!("\t            single_band_tracepoints: {:?}", single_band_tracepoints);
                        error!("\t            double_band_tracepoints: {:?}", double_band_tracepoints);
                        error!("\t                     CIGAR from PAF: {}", paf_cigar);
                        error!("\t             CIGAR from tracepoints: {}", cigar_from_tracepoints);
                        error!("\t CIGAR from single_band_tracepoints: {}", cigar_from_single_band_tracepoints);
                        error!("\t CIGAR from double_band_tracepoints: {}", cigar_from_double_band_tracepoints);
                        error!("\t seqa: {}", String::from_utf8(query_seq).unwrap());
                        error!("\t seqb: {}", String::from_utf8(target_seq).unwrap());
                        error!("                    bounds CIGAR from PAF: {:?}", get_cigar_diagonal_bounds(&paf_cigar));
                        error!("            bounds CIGAR from tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_tracepoints));
                        error!("bounds CIGAR from single_band_tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_single_band_tracepoints));
                        error!("bounds CIGAR from double_band_tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_double_band_tracepoints));

                        let (deviation, d_min, d_max, max_gap) = compute_deviation(&cigar_from_tracepoints);
                        error!("                    deviation CIGAR from PAF: {:?}", compute_deviation(&paf_cigar));
                        error!("            deviation CIGAR from tracepoints: {:?}", (deviation, d_min, d_max, max_gap));
                        error!("deviation CIGAR from single_band_tracepoints: {:?}", compute_deviation(&cigar_from_single_band_tracepoints));
                        error!("deviation CIGAR from double_band_tracepoints: {:?}", compute_deviation(&cigar_from_double_band_tracepoints));
                        error!("=> Try using --wfa-heuristic=banded-static --wfa-heuristic-parameters=-{},{}\n", std::cmp::max(max_gap, -d_min), std::cmp::max(max_gap, d_max));
                        std::process::exit(1);
                    }
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
                let mut aligner = AffineWavefronts::default();
                let paf_cigar = align_sequences_wfa(&query_seq[a_start..a_end], &target_seq[b_start..b_end], &mut aligner);
                let paf_cigar = cigar_ops_to_cigar_string(&paf_cigar);

                let tracepoints = cigar_to_tracepoints(&paf_cigar, max_diff);
                let recon_cigar_variable = tracepoints_to_cigar(
                    &tracepoints,
                    &query_seq,
                    &target_seq,
                    0,
                    0,
                    mismatch,
                    gap_open1,
                    gap_ext1,
                    gap_open2,
                    gap_ext2
                );

                //let realn_cigar = align_sequences_wfa(&query_seq, &target_seq);
                //let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);

                info!("\t            tracepoints: {:?}", tracepoints);
                info!("\t CIGAR from tracepoints_variable: {}", recon_cigar_variable);
                info!("\t              CIGAR from the PAF: {}", paf_cigar);
                println!("get_cigar_diagonal_bounds from the PAF {:?}", get_cigar_diagonal_bounds(&paf_cigar));
                //info!("\t CIGAR from realignment: {}", realn_cigar);

                assert!(paf_cigar == recon_cigar_variable);
            }
        },
    }

    Ok(())
}

/// Initialize logger based on verbosity
fn setup_logger(verbosity: u8) {
    env_logger::Builder::new()
        .filter_level(match verbosity {
            0 => log::LevelFilter::Warn,    // Errors and warnings
            1 => log::LevelFilter::Info,    // Errors, warnings, and info
            _ => log::LevelFilter::Debug,   // Errors, warnings, info, and debug
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
fn process_compress_chunk(lines: &[String], band_type: u8, mixed: bool, max_diff: usize) {
    lines.par_iter().for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            error!("{}", message_with_truncate_paf_file("Skipping malformed PAF line", line));
            std::process::exit(1);
        }
        
        let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
            error!("{}", message_with_truncate_paf_file("Skipping CIGAR-less PAF line", line));
            std::process::exit(1);
        };
        let cigar = &cg_field[5..]; // Direct slice instead of strip_prefix("cg:Z:")
        
        // Convert CIGAR based on options and add type prefix
        let tracepoints_str = if mixed {
            // Use mixed representation (type 3)
            let tp = cigar_to_mixed_tracepoints(cigar, max_diff);
            format!("D{}", format_mixed_tracepoints(&tp))
        } else {
            // Use standard tracepoints with optional banding
            match band_type {
                0 => {
                    // No-band tracepoints (type 0)
                    let tp = cigar_to_tracepoints(cigar, max_diff);
                    format!("A{}", format_tracepoints(&tp))
                },
                1 => {
                    // Single-band tracepoints (type 1)
                    let tp = cigar_to_single_band_tracepoints(cigar, max_diff);
                    format!("B{}", format_single_band_tracepoints(&tp))
                },
                2 => {
                    // Double-band tracepoints (type 2)
                    let tp = cigar_to_double_band_tracepoints(cigar, max_diff);
                    format!("C{}", format_double_band_tracepoints(&tp))
                },
                _ => unreachable!()
            }
        };
       
        // Print the result
        let new_line = line.replace(cg_field, &format!("tp:Z:{}", tracepoints_str));
        println!("{}", new_line); // It is thread-safe by default
    });
}

/// Process a chunk of lines in parallel for decompression
fn process_decompress_chunk(
    lines: &[String], 
    query_fasta_path: &str,
    target_fasta_path: &str, 
    mismatch: i32, 
    gap_open1: i32, 
    gap_ext1: i32, 
    gap_open2: i32, 
    gap_ext2: i32
) {
    lines.par_iter().for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            error!("{}", message_with_truncate_paf_file("Skipping malformed PAF line", line));
            std::process::exit(1);
        }
        
        let Some(tp_field) = fields.iter().find(|&&s| s.starts_with("tp:Z:")) else {
            error!("{}", message_with_truncate_paf_file("Skipping tracepoints-less PAF line", line));
            std::process::exit(1);
        };
        let tracepoints_with_type = &tp_field[5..]; // Direct slice instead of strip_prefix("tp:Z:")
        // Extract type and tracepoints string
        let tp_type = tracepoints_with_type.chars().next().unwrap();
        let tracepoints_str = &tracepoints_with_type[1..]; // Skip the first character
        
        // Parse mandatory PAF fields
        let query_name = fields[0];
        let query_start: usize = fields[2].parse().unwrap_or_else(|_| {
            error!("{}", message_with_truncate_paf_file("Invalid query_start in PAF line", line));
            std::process::exit(1);
        });
        let query_end: usize = fields[3].parse().unwrap_or_else(|_| {
            error!("{}", message_with_truncate_paf_file("Invalid query_end in PAF line", line));
            std::process::exit(1);
        });
        let strand = fields[4];
        let target_name = fields[5];
        //let target_len: usize = fields[6].parse()?;
        let target_start: usize = fields[7].parse().unwrap_or_else(|_| {
            error!("{}", message_with_truncate_paf_file("Invalid target_start in PAF line", line));
            std::process::exit(1);
        });
        let target_end: usize = fields[8].parse().unwrap_or_else(|_| {
            error!("{}", message_with_truncate_paf_file("Invalid target_end in PAF line", line));
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
        
        // Fetch query sequence from query FASTA
        let query_seq = if strand == "+" {
            match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut seq_vec = seq.to_vec();
                    unsafe {libc::free(seq.as_ptr() as *mut std::ffi::c_void)}; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                    seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                    seq_vec
                },
                Err(e) => {
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
                }
            }
        } else {
            match query_fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut rc = reverse_complement(&seq.to_vec());
                    unsafe {libc::free(seq.as_ptr() as *mut std::ffi::c_void)}; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                    rc.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                    rc
                },
                Err(e) => {
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
                }
            }
        };

        // Fetch target sequence from target FASTA
        let target_seq = match target_fasta_reader.fetch_seq(target_name, target_start, target_end - 1) {
            Ok(seq) => {
                let mut seq_vec = seq.to_vec();
                unsafe {libc::free(seq.as_ptr() as *mut std::ffi::c_void)}; // Free up memory (bug https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
                seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                seq_vec
            },
            Err(e) => {
                error!("Failed to fetch target sequence: {}", e);
                std::process::exit(1);
            }
        };

        // Convert to CIGAR based on the type prefix
        let cigar = match tp_type {
            'A' => {
                // No-band tracepoints
                let tracepoints = parse_tracepoints(tracepoints_str);
                tracepoints_to_cigar(
                    &tracepoints,
                    &query_seq,
                    &target_seq,
                    0, 0,
                    mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
                )
            },
            'B' => {
                // Single-band tracepoints
                let tracepoints = parse_single_band_tracepoints(tracepoints_str);
                single_band_tracepoints_to_cigar(
                    &tracepoints,
                    &query_seq,
                    &target_seq,
                    0, 0,
                    mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
                )
            },
            'C' => {
                // Double-band tracepoints
                let tracepoints = parse_double_band_tracepoints(tracepoints_str);
                double_band_tracepoints_to_cigar(
                    &tracepoints,
                    &query_seq,
                    &target_seq,
                    0, 0,
                    mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
                )
            },
            'D' => {
                // Mixed representation
                let mixed_tracepoints = parse_mixed_tracepoints(tracepoints_str);
                mixed_tracepoints_to_cigar(
                    &mixed_tracepoints,
                    &query_seq,
                    &target_seq,
                    0, 0,
                    mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
                )
            },
            _ => {
                error!("Invalid tracepoint type '{}' in PAF line", tp_type);
                std::process::exit(1);
            }
        };

        // Print the original line, replacing the tracepoints tag with the CIGAR string
        let new_line = line.replace(tp_field, &format!("cg:Z:{}", cigar));
        println!("{}", new_line);
    });
}

/// Combines a message with the first 9 columns of a PAF line.
fn message_with_truncate_paf_file(message: &str, line: &str) -> String {
    let truncated_line = line
        .split('\t')
        .take(9)
        .collect::<Vec<&str>>()
        .join("\t");
    format!("{}: {} ...", message, truncated_line)
}

fn format_tracepoints(tracepoints: &[(usize, usize)]) -> String {
    tracepoints.iter()
        .map(|(a, b)| format!("{},{}", a, b))
        .collect::<Vec<String>>()
        .join(";")
}
fn format_single_band_tracepoints(tracepoints: &[(usize, usize, usize)]) -> String {
    tracepoints.iter()
        .map(|(a, b, c)| format!("{},{},{}", a, b, c))
        .collect::<Vec<String>>()
        .join(";")
}
fn format_double_band_tracepoints(tracepoints: &[(usize, usize, (isize, isize))]) -> String {
    tracepoints.iter()
        .map(|(a, b, (c, d))| format!("{},{},{},{}", a, b, c, d))
        .collect::<Vec<String>>()
        .join(";")
}
fn format_mixed_tracepoints(mixed_tracepoints: &[MixedRepresentation]) -> String {
    mixed_tracepoints.iter()
        .map(|tp| match tp {
            MixedRepresentation::Tracepoint(a, b) => format!("{},{}", a, b),
            MixedRepresentation::CigarOp(len, op) => format!("{}{}", len, op),
        })
        .collect::<Vec<String>>()
        .join(";")
}

fn parse_penalties(penalties: &str) -> Result<(i32, i32, i32, i32, i32), Box<dyn std::error::Error>> {
    let tokens: Vec<&str> = penalties.split(',').collect();
    if tokens.len() != 5 {
        error!("Error: penalties must be provided as mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2");
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
    tp_str.split(';')
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            Some((
                parts[0].parse().unwrap(),
                parts[1].parse().unwrap(),
            ))
        }).collect()
}
fn parse_single_band_tracepoints(tp_str: &str) -> Vec<(usize, usize, usize)> {
    tp_str.split(';')
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            Some((
                parts[0].parse().unwrap(),
                parts[1].parse().unwrap(),
                parts[2].parse().unwrap(),
            ))
        }).collect()
}
fn parse_double_band_tracepoints(tp_str: &str) -> Vec<(usize, usize, (isize, isize))> {
    tp_str.split(';')
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            Some((
                parts[0].parse().unwrap(),
                parts[1].parse().unwrap(),
                (parts[2].parse().unwrap(), parts[3].parse().unwrap()),
            ))
        }).collect()
}
fn parse_mixed_tracepoints(tp_str: &str) -> Vec<MixedRepresentation> {
    tp_str.split(';')
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

fn get_cigar_diagonal_bounds(cigar: &str) -> (i64, i64) {
    let mut current_diagonal = 0;  // Current diagonal position
    let mut min_diagonal = 0;      // Lowest diagonal reached
    let mut max_diagonal = 0;      // Highest diagonal reached
    
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
                },
                'D' => {
                    // Deletions move down diagonal by count amount
                    current_diagonal -= count;
                    min_diagonal = min_diagonal.min(current_diagonal);
                },
                'I' => {
                    // Insertions move up diagonal by count amount
                    current_diagonal += count;
                    max_diagonal = max_diagonal.max(current_diagonal);
                },
                _ => panic!("Invalid CIGAR operation: {}", c)
            }
        }
    }
    
    (min_diagonal, max_diagonal)
}

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
                },
                'D' => {
                    // Deletions move down diagonal by count amount
                    deviation -= count;
                    max_gap = std::cmp::max(max_gap, count);
                },
                'I' => {
                    deviation += count;
                    max_gap = std::cmp::max(max_gap, count);
                },
                _ => panic!("Invalid CIGAR operation: {}", c)
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
            _ => b'N',  // Convert any unexpected bases to N
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
