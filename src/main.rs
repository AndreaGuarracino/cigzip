use clap::Parser;
use std::io::{self, BufRead, BufReader};
use std::fs::File;
use flate2::read::GzDecoder;
use rust_htslib::faidx::Reader as FastaReader;
use lib_tracepoints::{cigar_to_tracepoints, cigar_to_banded_tracepoints, tracepoints_to_cigar, banded_tracepoints_to_cigar};
use lib_wfa2::affine_wavefront::{AffineWavefronts};
use log::{info, warn, error};
use rayon::prelude::*;

/// Common options shared between all commands
#[derive(Parser, Debug)]
struct CommonOpts {
    /// PAF file for alignments (use "-" to read from standard input)
    #[arg(short = 'p', long = "paf")]
    paf: String,

    /// Number of threads to use (default: 2)
    #[arg(short = 't', long = "threads", default_value_t = 2)]
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

        /// Use banded tracepoints instead of regular tracepoints
        #[arg(short = 'b', long = "banded")]
        banded: bool,
        
        /// Max-diff value for tracepoints
        #[arg(long, default_value = "128")]
        max_diff: usize,
    },
    /// Decompression of alignments
    Decompress {
        #[clap(flatten)]
        common: CommonOpts,
        
        /// FASTA file for sequences
        #[arg(short = 'f', long = "fasta")]
        fasta: String,
        
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

        /// FASTA file for sequences
        #[arg(short = 'f', long = "fasta")]
        fasta: Option<String>,

        /// Gap penalties in the format mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2
        #[arg(long, default_value = "3,4,2,24,1")]
        penalties: String,

        /// Max-diff value for tracepoints
        #[arg(long, default_value = "128")]
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
        Args::Compress { common, banded, max_diff } => {
            setup_logger(common.verbose);
            info!("Converting CIGAR to {} tracepoints", if banded { "banded" } else { "basic" });

            // Set the thread pool size
            rayon::ThreadPoolBuilder::new()
                .num_threads(common.threads)
                .build_global()?;

            // Open the PAF file (or use stdin if "-" is provided).
            let paf_reader = get_paf_reader(&common.paf)?;

            // Process in chunks
            const CHUNK_SIZE: usize = 1000;
            let mut lines = Vec::with_capacity(CHUNK_SIZE);
            
            for line_result in paf_reader.lines() {
                match line_result {
                    Ok(line) => {
                        if line.trim().is_empty() || line.starts_with('#') {
                            continue;
                        }
                        
                        lines.push(line);
                        
                        if lines.len() >= CHUNK_SIZE {
                            // Process current chunk in parallel
                            process_compress_chunk(&lines, banded, max_diff);
                            lines.clear();
                        }
                    },
                    Err(e) => return Err(e.into()),
                }
            }
            
            // Process remaining lines
            if !lines.is_empty() {
                process_compress_chunk(&lines, banded, max_diff);
            }
        },
        Args::Decompress { common, fasta, penalties } => {
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
            const CHUNK_SIZE: usize = 1000;
            let mut lines = Vec::with_capacity(CHUNK_SIZE);

            for line_result in paf_reader.lines() {
                match line_result {
                    Ok(line) => {
                        if line.trim().is_empty() || line.starts_with('#') {
                            continue;
                        }
                        
                        lines.push(line);
                        
                        if lines.len() >= CHUNK_SIZE {
                            // Process current chunk in parallel
                            process_decompress_chunk(
                                &lines, 
                                &fasta, 
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
                    &fasta, 
                    mismatch, 
                    gap_open1, 
                    gap_ext1, 
                    gap_open2, 
                    gap_ext2
                );
            }
        },
        #[cfg(debug_assertions)]
        Args::Debug { paf, fasta, penalties, max_diff, verbose } => {
            setup_logger(verbose);
            info!("Debugging");

            let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = parse_penalties(&penalties)?;
            info!("Penalties: {},{},{},{},{}", mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);

            if let (Some(paf), Some(fasta)) = (paf, fasta) {
                info!("PAF file: {}", paf);
                info!("FASTA file: {}", fasta);
                info!("Penalties: {}", penalties);

                // Open the FASTA file
                let fasta_reader = FastaReader::from_path(&fasta).expect("Error reading FASTA file");

                // Open the PAF file (or use stdin if "-" is provided).
                let reader: Box<dyn BufRead> = if paf == "-" {
                    Box::new(BufReader::new(std::io::stdin()))
                } else {
                    Box::new(BufReader::new(File::open(&paf)?))
                };

                for line in reader.lines() {
                    let line = line?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 12 {
                        warn!("Skipping malformed PAF line: {}", line);
                        continue;
                    }
                    let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
                        warn!("Skipping CIGAR-less PAF line: {}", line);
                        continue;
                    };
                    let paf_cigar = &cg_field[5..]; // Direct slice instead of strip_prefix("cg:Z:")

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
                    
                    // Fetch query sequence from FASTA.
                    let query_seq = if strand == "+" {
                        fasta_reader.fetch_seq(query_name, query_start, query_end - 1)?.to_vec()
                    } else {
                        // For reverse strand, fetch the sequence and reverse complement it
                        reverse_complement(&fasta_reader.fetch_seq(query_name, query_start, query_end - 1)?.to_vec())
                    };
                    // Fetch target sequence from FASTA.
                    let target_seq = fasta_reader.fetch_seq(target_name, target_start, target_end - 1)?.to_vec();
                    
                    // let mut aligner = AffineWavefronts::with_penalties_affine2p(0, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
                    // let realn_cigar = align_sequences_wfa(&query_seq, &target_seq, &mut aligner);
                    // let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);
                    // let paf_cigar = &realn_cigar;

                    // Convert CIGAR to tracepoints using query (A) and target (B) coordinates.
                    let tracepoints = cigar_to_tracepoints(paf_cigar, max_diff);
                    let banded_tracepoints = cigar_to_banded_tracepoints(paf_cigar, max_diff);

                    // Compare tracepoints and banded_tracepoints, but only the first 2 elements of the tuple. if different, print error
                    if tracepoints.iter().zip(tracepoints.iter()).any(|(a, b)| a.0 != b.0 || a.1 != b.1) {
                        error!("Tracepoints mismatch! {}", line);
                        error!("\t       tracepoints: {:?}", tracepoints);
                        error!("\tbanded_tracepoints: {:?}", banded_tracepoints);
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
                    let cigar_from_banded_tracepoints = banded_tracepoints_to_cigar(
                        &banded_tracepoints,
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
                    // let (query_end_variable2, query_len_variable2, target_end_variable2, target_len_variable2) = calculate_alignment_coordinates(&cigar_from_banded_tracepoints, query_start, target_start);
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

                    if cigar_from_tracepoints != cigar_from_banded_tracepoints {
                        error!("CIGAR mismatch! {}", line);
                        info!("{}", line);
                        info!("\t                   tracepoints: {:?}", tracepoints);
                        info!("\t            banded_tracepoints: {:?}", banded_tracepoints);
                        info!("\t                CIGAR from PAF: {}", paf_cigar);
                        info!("\t        CIGAR from tracepoints: {}", cigar_from_tracepoints);
                        info!("\t CIGAR from banded_tracepoints: {}", cigar_from_banded_tracepoints);
                        info!("\t seqa: {}", String::from_utf8(query_seq).unwrap());
                        info!("\t seqb: {}", String::from_utf8(target_seq).unwrap());
                        info!("               bounds CIGAR from PAF: {:?}", get_cigar_diagonal_bounds(&paf_cigar));
                        info!("       bounds CIGAR from tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_tracepoints));
                        info!("bounds CIGAR from banded_tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_banded_tracepoints));
                    }
                }
            } else {
                // Fallback: run default example if no PAF/FASTA provided.
                info!("No PAF and FASTA provided, running default example.");

                let a = get_cigar_diagonal_bounds("514M1D366M1X104M1X92M1X197M1X261M3I664M1X126M1X17M1X572M1X525M1X234M5D320M1X584M1X221M1X64M1X231M1X971M1X586M1X324M1X305M1X344M1X38M1X14M1X109M1X283M1X211M1X1011M1X77M1X33M1X136M1X1190M1X779M1X302M1X420M1X539M1X317M1X2M1X782M1X425M1X496M1X35M1X723M1X146M1X77M1X454M1X16M2D8M1X146M2X14M1X208M158I15M1X6M1X4M2X5M1X2M1X3M2X3M2X1M1X22M1I2M1X6M2X2M1X1M1X7M3X6M2X1M1X4M1X10M2X25M1D6M5D10M1X1M2I4M2X1M235D1M");
                println!("from_paf {:?}", a);

                let query_seq: String = "GAACAGAGAAATGGTGGAATTCAAATACAAAAAAACCGCAAAATTAAAAATCTTGCGGCTCTCTGAACTCATTTTCATGAGTGAATTTGGCGGAACGGACGGGACTCGAACCCGCGACCCCCTGCGTGACAGGCAGGTATTCTAACCGACTGAACTACCGCTCCGCCGTTGTGTTCCGTTGGGAACGGGCGAATATTACGGATTTGCCTCACCCTTCGTCAACGGTTTTTCTCATCTTTTGAATCGTTTGCTGCAAAAATCGCCCAAGCCGCTATTTTTAGCGCCTTTTACAGGTATTTATGCCCGCCAGAGGCAGCTTCCGCCCTTCTTCTCCACCAGATCAAGACGGGCTTCCTGAGCTGCAAGCTCTTCATCTGTCGCAAAAACAACGCGTAACTTACTTGCCTGACGTACAATGCGCTGAATTGTTGCTTCACCTTGTTGCTGCTGTGTCTCTCCTTCCATCGCAAAAGCCATCGACGTTTGACCACCGGTCATCG".to_owned();
                let target_seq: String = "GAACAGAGAAATGGTGGAATTCAAATACAAAAAAACCGCAAAATTAACCCTTCGTCAACGGTTTTTCTCATCTTTTGAATCGTTTGCTGCAAAAATCGCCCAAGCCGCTATTTTTAGCGCCTTTTACAGGTATTTATGCCCGCCAGAGGCAGCTTCCGCCCTTCTTCTCCACCAGATCAAGACGGGCTTCCTGAGCTGCAAGCTCTTCATCTGTCGCAAAAACAACGCGTAACTTACTTGCCTGACGTACAATGCGCTGAATTGTTGCTTCACCTTGTTGCTGCTGTGTCTCTCCTTCCATCGCAAAAGCCATCGACGTTTGACCACCGGTCATCG".to_owned();

                let a_start = 0;
                let a_end = query_seq.len();
                let b_start = 0;
                let b_end = target_seq.len();

                // Create aligner and configure settings
                // let mut aligner = AffineWavefronts::default();
                // let paf_cigar = align_sequences_wfa(&query_seq[a_start..a_end], &target_seq[b_start..b_end], &mut aligner);
                // let paf_cigar = cigar_ops_to_cigar_string(&paf_cigar);

                // let tracepoints = cigar_to_tracepoints_variable(&paf_cigar, max_diff);

                // let recon_cigar_variable = tracepoints_to_cigar_variable(
                //     &tracepoints_vartracepointsable,
                //     &query_seq,
                //     &target_seq,
                //     0,
                //     0,
                //     mismatch,
                //     gap_open1,
                //     gap_ext1,
                //     gap_open2,
                //     gap_ext2
                // );

                // //let realn_cigar = align_sequences_wfa(&query_seq, &target_seq);
                // //let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);

                // info!("\t            tracepoints: {:?}", tracepoints);
                // info!("\t CIGAR from tracepoints_variable: {}", recon_cigar_variable);
                // info!("\t              CIGAR from the PAF: {}", paf_cigar);
                // //info!("\t CIGAR from realignment: {}", realn_cigar);

                // assert!(paf_cigar == recon_cigar_variable);
            }
        },
    }

    Ok(())
}

/// Initialize logger based on verbosity
fn setup_logger(verbosity: u8) {
    env_logger::Builder::new()
        .filter_level(match verbosity {
            0 => log::LevelFilter::Error,
            1 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        })
        .init();
}

fn get_paf_reader(paf: &str) -> io::Result<Box<dyn BufRead>> {
    if paf == "-" {
        Ok(Box::new(BufReader::new(std::io::stdin())))
    } else if paf.ends_with(".gz") || paf.ends_with(".bgz") {
        let file = File::open(paf)?;
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        let file = File::open(paf)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Process a chunk of lines in parallel for compression
fn process_compress_chunk(lines: &[String], banded: bool, max_diff: usize) {
    lines.par_iter().for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            warn!("Skipping malformed PAF line: {}", line);
            return;
        }
        
        let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
            warn!("Skipping CIGAR-less PAF line: {}", line);
            return;
        };
        let cigar = &cg_field[5..]; // Direct slice instead of strip_prefix("cg:Z:")
        
        // Convert CIGAR to tracepoints
        let tracepoints_str = if banded {
            let tp = cigar_to_banded_tracepoints(cigar, max_diff);
            format_banded_tracepoints(&tp)
        } else {
            let tp = cigar_to_tracepoints(cigar, max_diff);
            format_tracepoints(&tp)
        };
        
        // Print the result
        let new_line = line.replace(cg_field, &format!("tp:Z:{}", tracepoints_str));
        println!("{}", new_line); // It is thread-safe by default
    });
}
/// Process a chunk of lines in parallel for decompression
fn process_decompress_chunk(
    lines: &[String], 
    fasta_path: &str, 
    mismatch: i32, 
    gap_open1: i32, 
    gap_ext1: i32, 
    gap_open2: i32, 
    gap_ext2: i32
) {
    lines.par_iter().for_each(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            warn!("Skipping malformed PAF line: {}", line);
            return;
        }
        
        let Some(tp_field) = fields.iter().find(|&&s| s.starts_with("tp:Z:")) else {
            warn!("Skipping tracepoints-less PAF line: {}", line);
            return;
        };
        let tracepoints_str = &tp_field[5..]; // Direct slice instead of strip_prefix("tp:Z:")
        
        // Parse mandatory PAF fields
        let query_name = fields[0];
        let query_start: usize = match fields[2].parse() {
            Ok(val) => val,
            Err(_) => {
                warn!("Invalid query_start in PAF line: {}", line);
                return;
            }
        };
        let query_end: usize = match fields[3].parse() {
            Ok(val) => val,
            Err(_) => {
                warn!("Invalid query_end in PAF line: {}", line);
                return;
            }
        };
        let strand = fields[4];
        let target_name = fields[5];
        let target_start: usize = match fields[7].parse() {
            Ok(val) => val,
            Err(_) => {
                warn!("Invalid target_start in PAF line: {}", line);
                return;
            }
        };
        let target_end: usize = match fields[8].parse() {
            Ok(val) => val,
            Err(_) => {
                warn!("Invalid target_end in PAF line: {}", line);
                return;
            }
        };

        // Create a thread-local FASTA reader
        let fasta_reader = match FastaReader::from_path(fasta_path) {
            Ok(reader) => reader,
            Err(e) => {
                error!("Failed to create FASTA reader: {}", e);
                return;
            }
        };

        // Fetch query sequence
        let query_seq = if strand == "+" {
            match fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut seq_vec = seq.to_vec();
                    seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                    seq_vec
                },
                Err(e) => {
                    warn!("Failed to fetch query sequence: {}", e);
                    return;
                }
            }
        } else {
            match fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut rc = reverse_complement(&seq.to_vec());
                    rc.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                    rc
                },
                Err(e) => {
                    warn!("Failed to fetch query sequence: {}", e);
                    return;
                }
            }
        };

        // Fetch target sequence
        let target_seq = match fasta_reader.fetch_seq(target_name, target_start, target_end - 1) {
            Ok(seq) => {
                let mut seq_vec = seq.to_vec();
                seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                seq_vec
            },
            Err(e) => {
                warn!("Failed to fetch target sequence: {}", e);
                return;
            }
        };

        // Check if the tracepoints are banded or not
        let is_banded = tracepoints_str.split(';').next().unwrap().split(',').count() == 4;
        
        // Convert tracepoints to CIGAR
        let cigar = if is_banded {
            let tracepoints = parse_banded_tracepoints(tracepoints_str);
            banded_tracepoints_to_cigar(
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
            )
        } else {
            let tracepoints = parse_tracepoints(tracepoints_str);
            tracepoints_to_cigar(
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
            )
        };

        // Print the original line, replacing the tracepoints tag with the CIGAR string
        let new_line = line.replace(tp_field, &format!("cg:Z:{}", cigar));
        println!("{}", new_line);
    });
}

fn format_tracepoints(tracepoints: &[(usize, usize)]) -> String {
    tracepoints.iter()
        .map(|(a, b)| format!("{},{}", a, b))
        .collect::<Vec<String>>()
        .join(";")
}
fn format_banded_tracepoints(tracepoints: &[(usize, usize, (isize, isize))]) -> String {
    tracepoints.iter()
        .map(|(a, b, (c, d))| format!("{},{},{},{}", a, b, c, d))
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
fn parse_banded_tracepoints(tp_str: &str) -> Vec<(usize, usize, (isize, isize))> {
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
            // Get the count (or 1 if no number specified)
            let count = if num_buffer.is_empty() {
                1
            } else {
                num_buffer.parse::<i64>().unwrap()
            };
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
