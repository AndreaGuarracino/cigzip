use core::panic;
use clap::Parser;
use std::io::{BufRead, BufReader};
use std::fs::File;
use rust_htslib::faidx::Reader as FastaReader;

use lib_tracepoints::{cigar_to_tracepoints, cigar_to_banded_tracepoints, cigar_from_tracepoints, cigar_from_banded_tracepoints};
use lib_wfa2::affine_wavefront::{AffineWavefronts};
use log::{info, warn, error};

/// Command-line arguments parsed with Clap.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "CIGAR Alignment and Tracepoint Reconstruction",
    long_about = None
)]
struct Args {
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
    #[arg(long, default_value = "128")]//, value_parser = clap::value_parser!(usize).range(1..usize::MAX))]
    max_diff: usize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "0")]
    verbose: u8,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments.
    let args = Args::parse();

    // Initialize logger based on verbosity
    env_logger::Builder::new()
        .filter_level(match args.verbose {
            0 => log::LevelFilter::Error,
            1 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        })
        .init();

    let tokens: Vec<&str> = args.penalties.split(',').collect();
    if tokens.len() != 5 {
        error!("Error: penalties must be provided as mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2");
        std::process::exit(1);
    }
    let mismatch: i32 = tokens[0].parse()?;
    let gap_open1: i32 = tokens[1].parse()?;
    let gap_ext1: i32 = tokens[2].parse()?;
    let gap_open2: i32 = tokens[3].parse()?;
    let gap_ext2: i32 = tokens[4].parse()?;
    info!("Penalties: {},{},{},{},{}", mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);

    let max_diff = args.max_diff;

    // If both a PAF file and FASTA file are provided, process each PAF record.
    if let (Some(paf_path), Some(fasta_path)) = (args.paf, args.fasta) {
        info!("Using FASTA file: {}", fasta_path);
        info!("Processing PAF file: {}", paf_path);

        // Open the FASTA file
        let fasta_reader = FastaReader::from_path(&fasta_path)
            .expect("Error reading FASTA file");

        // Open the PAF file (or use stdin if "-" is provided).
        let reader: Box<dyn BufRead> = if paf_path == "-" {
            Box::new(BufReader::new(std::io::stdin()))
        } else {
            Box::new(BufReader::new(File::open(&paf_path)?))
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
                String::from_utf8(fasta_reader.fetch_seq(query_name, query_start, query_end - 1)?.to_vec())?
            } else {
                // For reverse strand, fetch the sequence and reverse complement it
                reverse_complement(&String::from_utf8(
                    fasta_reader.fetch_seq(query_name, query_start, query_end - 1)?.to_vec()
                )?)
            };
            // Fetch target sequence from FASTA.
            let target_seq = String::from_utf8(fasta_reader.fetch_seq(target_name, target_start, target_end - 1)?.to_vec())?;
            
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
            let recon_cigar_from_tracepoints = cigar_from_tracepoints(
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
            let recon_cigar_from_banded_tracepoints = cigar_from_banded_tracepoints(
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

            // let (query_end_variable, query_len_variable, target_end_variable, target_len_variable) = calculate_alignment_coordinates(&recon_cigar_from_tracepoints, query_start, target_start);
            // let (query_end_variable2, query_len_variable2, target_end_variable2, target_len_variable2) = calculate_alignment_coordinates(&recon_cigar_from_banded_tracepoints, query_start, target_start);
            // let (query_end_paf, query_len_paf, target_end_paf, target_len_paf) = calculate_alignment_coordinates(paf_cigar, query_start, target_start);
            // if (query_len_paf != query_len_variable) || (target_len_paf != target_len_variable) || (query_end_paf != query_end_variable) || (target_end_paf != target_end_variable) {
            //     info!("recon_cigar_from_tracepoints {:?}", (query_end_variable, query_len_variable, target_end_variable, target_len_variable));
            //     info!("           paf_cigar {:?}", (query_end_paf, query_len_paf, target_end_paf, target_len_paf) );
            //     error!("Line {}: seq. len. mismatch!", line);
            // }
            // if (query_len_paf != query_len_variable2) || (target_len_paf != target_len_variable2) || (query_end_paf != query_end_variable2) || (target_end_paf != target_end_variable2) {
            //     info!("recon_cigar_from_banded_tracepoints {:?}", (query_end_variable2, query_len_variable2, target_end_variable2, target_len_variable2));
            //     info!("           p af_cigar {:?}", (query_end_paf, query_len_paf, target_end_paf, target_len_paf) );
            //     error!("Line {}: seq. len. mismatch!", line);
            // }

            if recon_cigar_from_tracepoints != recon_cigar_from_banded_tracepoints {
                error!("CIGAR mismatch! {}", line);
                info!("{}", line);
                info!("\t                   tracepoints: {:?}", tracepoints);
                info!("\t            banded_tracepoints: {:?}", banded_tracepoints);
                info!("\t                CIGAR from PAF: {}", paf_cigar);
                info!("\t        CIGAR from tracepoints: {}", recon_cigar_from_tracepoints);
                info!("\t CIGAR from banded_tracepoints: {}", recon_cigar_from_banded_tracepoints);
                info!("\t seqa: {}", query_seq);
                info!("\t seqb: {}", target_seq);
                info!("               bounds CIGAR from PAF: {:?}", get_cigar_diagonal_bounds(&paf_cigar));
                info!("       bounds CIGAR from tracepoints: {:?}", get_cigar_diagonal_bounds(&recon_cigar_from_tracepoints));
                info!("bounds CIGAR from banded_tracepoints: {:?}", get_cigar_diagonal_bounds(&recon_cigar_from_banded_tracepoints));
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

    Ok(())
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
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' => 'N',
            _ => c  // Keep other characters unchanged
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
