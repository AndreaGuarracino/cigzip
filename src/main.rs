use clap::Parser;
use std::io::{self, BufRead, BufReader};
use std::fs::File;
use flate2::read::GzDecoder;
use rust_htslib::faidx::Reader as FastaReader;
use lib_tracepoints::{cigar_to_tracepoints, cigar_to_banded_tracepoints, tracepoints_to_cigar, banded_tracepoints_to_cigar, align_sequences_wfa, cigar_ops_to_cigar_string};
use lib_wfa2::affine_wavefront::{AffineWavefronts};
use log::{info, error};
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

                // Open the FASTA file
                let fasta_reader = FastaReader::from_path(&fasta).expect("Error reading FASTA file");

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
                    
                    // Fetch query sequence
                    let query_seq = if strand == "+" {
                        match fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                            Ok(seq) => {
                                let mut seq_vec = seq.to_vec();
                                seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                                seq_vec
                            },
                            Err(e) => {
                                error!("Failed to fetch query sequence: {}", e);
                                std::process::exit(1);
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
                                error!("Failed to fetch query sequence: {}", e);
                                std::process::exit(1);
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
                    let banded_tracepoints = cigar_to_banded_tracepoints(paf_cigar, max_diff);

                    // Compare tracepoints and banded_tracepoints, but only the first 2 elements of the tuple. if different, print error
                    if tracepoints.iter().zip(tracepoints.iter()).any(|(a, b)| a.0 != b.0 || a.1 != b.1) {
                        error!("Tracepoints mismatch! {}", line);
                        error!("\t       tracepoints: {:?}", tracepoints);
                        error!("\tbanded_tracepoints: {:?}", banded_tracepoints);
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
                        error!("\t                   tracepoints: {:?}", tracepoints);
                        error!("\t            banded_tracepoints: {:?}", banded_tracepoints);
                        error!("\t                CIGAR from PAF: {}", paf_cigar);
                        error!("\t        CIGAR from tracepoints: {}", cigar_from_tracepoints);
                        error!("\t CIGAR from banded_tracepoints: {}", cigar_from_banded_tracepoints);
                        error!("\t seqa: {}", String::from_utf8(query_seq).unwrap());
                        error!("\t seqb: {}", String::from_utf8(target_seq).unwrap());
                        error!("               bounds CIGAR from PAF: {:?}", get_cigar_diagonal_bounds(&paf_cigar));
                        error!("       bounds CIGAR from tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_tracepoints));
                        error!("bounds CIGAR from banded_tracepoints: {:?}", get_cigar_diagonal_bounds(&cigar_from_banded_tracepoints));

                        let (deviation, d_min, d_max, max_gap) = compute_deviation(&cigar_from_tracepoints);
                        error!("               deviation CIGAR from PAF: {:?}", compute_deviation(&paf_cigar));
                        error!("       deviation CIGAR from tracepoints: {:?}", (deviation, d_min, d_max, max_gap));
                        error!("deviation CIGAR from banded_tracepoints: {:?}", compute_deviation(&cigar_from_banded_tracepoints));
                        error!("=> Try using --wfa-heuristic=banded-static --wfa-heuristic-parameters=-{},{}\n", std::cmp::max(max_gap, -d_min), std::cmp::max(max_gap, d_max));
                        std::process::exit(1);
                    }
                }
            } else {
                // Fallback: run default example if no PAF/FASTA provided.
                info!("No PAF and FASTA provided, running default example.");

                let a = get_cigar_diagonal_bounds("1I5=1I10=4I11=1X3=3I25=1I3=1I16=1D6=2X3=1D7=1I8=3I19=7I32=1X1=15I2=1X15=2X11=1D3X10=1X3=1X7=1D14=1X15=1X1=1X2=1X8=1I14=1X5=1X15=1X3=5I50=1X1=1X11=1X5=1X18=1X6=2X11=1D4=1X16=1X5=1X2=1D4=1D1X6=1X2=1X19=2X2=11D4=1X4=1X20=1X2=2I8=1X29=1I1=1X22=1X6=3X8=1D2=1X16=1D20=1X41=1X5=1X22=1X2=1X31=1X28=1X7=1X2=1X32=1X10=1X12=1X37=1X17=2X8=1X1=1X7=1X4=3X13=1I22=1X4=1X7=1X19=1X12=1X21=1X3=1X2=1X2=1X11=1X26=1D4=1X3=1X6=1X5=1X3=1X3=2X1=1X29=2D12=1X15=1X1=1X44=1X1=1X3=1X10=1X14=1X1=1X7=1X4=1X1=1X4=1X12=1X3=1X2=1X1=1X12=1X20=1X1=1X4=2I1X1=1X53=1X5=1X8=1X12=1X20=1X4=4D2X4=1X3=1X21=1X2=1D1X2=1X1=1D1X1=1X2=1X1=3X1=1X1=1X1=1X1=1X1=1I2=1X1=1X2=1X1=2X2=2X1=9I1=1X2=5X1=2X1=1X3=1X1=1X5=1I1=2X1=1X5=2I3X2=1X2=2X2=1X2=1D1=1X7=2X1=1X1=3X1=3X1=1X1=1X1=1I4=9I1X3=1X1=1X1=1X2=2X2=1D1X8=1X6=8D4=2X2=7X3=3X2=1X1=3D4=1X1=2X4=1X1=1D1=1X1=3X1=2X2=3X1=1X2=3X1=4X2=3X1=2X1=2X1=1I3=1I5=4I1X6=5D4=1X2=3D3=1X2=9D4=1D1X2=2X2=4D2X1=2X6=2I4=3X4=1X1=2X1=4X2=1D1X3=1X3=2X2=1X1=2D2X1=1X1=1X1=1X1=1X3=1X3=1X1=1X1=11X2=1X5=2D3=1X2=2X2=4X3=1I2=3I1X4=3X1=1X1=1X4=1X2=9I6=3X1=1X3=1X1=4I1X3=3X3=1X1=1X2=9I2=1X2=2I3=1X2=1D1X1=2X6=2X1=3I2=1X1=2X2=2X1=2X2=3X1=3X1=1X2=2X3=2X1=3D3X4=2X2=2X1=2D1X1=1X1=1X2=3X1=2X1=1I3X4=2X1=1X1=3X2=4X1=2X3=2X1=1X1=5X1=1X1=5X3=1I2X4=2X2=3I4=1X2=1X2=1I2X4=1X3=5I1=1X1=2X4=2X1=2X1=1X3=5X2=1D5=4X1=3X1=4X1I2=1X1=1X1=1X1=1X2=1X3=1D1=1X2=1X1=1X1=3X4=1X1=1D2=2D4X2=2X4=1X1=2D1=1X2=5D1=1X3=1X2=2X2=1X2=5D1=2X1=1X3=2X1=2X2=2X2=3X2=1X2=1I1X1=1X3=1X2=2X1=3I2=3I4=1X3=1X1=1X2=5I3=1X1=2X1I5=1X2=1X3I1=1X4=2I2=2X1=1X1=1X2=2X1=1X1=3X1=1X3=1X2=4D1X2=3X1=2X1=1X4=2X1=2I1X3=1X2=1X2=3D1=1X1=1X1=1X2=3X1=1X3=4D2=2X2=3X5=3D1X3=4D2=2X4=1D7=2I1=2X3=7D3=5D2=1X2=1X3=1X2=3X1=4X4=2X2=2I1X1=1X3=3X2=1X1=1X2=2I1=2X1=2X3=1X2=2D1=2X3=1I1=1X1=1X3=2X1=2X1=1X3=2I4=4I1X4=1X4=2D1=1X1=3X2=1X1=2D1=1X1=2X4=4D1=2X2=2X1=6X5=2X1=3X2=2D2=2X2=5D2X2=1X4=1I1X1=1X4=2D1=1X1=2X1=1X1=1X2=2X1=1X2=4D2=1X3=3X1=2X3=2X1=1X3=9D3=3X2=2X1=4X2=2X1=3I1=4X1=3X2=1X2=2X1=5X1=2X1=3X3=2X2=2X2=2X1=1X1=2X1=5I4=3X2=3X1=1X1=3X1=1X1=1X2=1X1=1X2=2X1=2X2=2X2=2X1=1X5=3I3X2=2X1=1X1=2X1=1X1=5X2=1X2=1X2=1X1=3X1=4X4=6D2=5D1=1X2=4X3=2X1=2X2=1D2=1X1=1X1=1I3=2X4=2X3=2I1=1X2=1X2=1X1=1X2=1X1=2X2=1X1=1X1=1X1=3X1=1X2=2I6=1I3=2X1=1X1=2X3=1X1=6X5=11D2=3D3=2X3=1D5=1X1=1X1=2I4X2=1X1=2X1=1I1=1X5=1X1=1X2=1X2=4X1=1I1X1=1X5=1X1=1D2=1X1=1X2=2X2=1X2=1I1X4=2I1=1X3=1D3=5I1X1=1X4=1X1=1X1=1X3=4I1X4=3X2=1I1=1X1=1X5=1X1=1X1=5D1X2=1X3=3X1=1X1=1X2=1I4=3X1=3X4=1D1X1=1X1=1X2=1X1=6D3=2X1=1X1=1X1=2X1=2X1=2I6=1X1=1I2=1I2=1X1=1X2=1X1=1X1=2X1=1I1X2=1X1=1X1=2X1=3X5=5X1=4X3=3I1=3X5=4I2=4X2=1X4=1X2=1D1X2=1X1=4X1=1X1=1X7=1I1=1X5=3X4=2X1=1X1=3X1=1X3=2I4=1X1=5X4=1X2=11D1=1X2=1X1=1X1=1X1=1X4=2X2=1X2=2X2=1X2=1X1=6D2=2X1=2X3=1X3=3D1=2X1=1X3=1I2=1X1=1X2=4X4=4X2=2X3=2I1=1X1=1X3=1X2=1I1=1X1=1X1=1X1=1X2=2X1=2D1=1X1=1X1=1X2=1X1=2X3=1X2=1I2X3=1X1=2D1X1=5X5=1X1=3X1=1X2=2X4=1X1=2X4=3I3=1D2=1X3=1X1=2X1=8D2=1D3=1X1=1X2=2X1=2X1=2X1=1D2X4=1X1=1X1=1D1X2=1X2=2X2=1X1=1X1=1X1=1X2=1D2=2D3=1X1=1X2=2X1=1X2=2D1X1=1X1=2X1=4X1=1X4=2X3=1D2X1=1X1=2X1=4X1=3X1=1X1=1D1=1X3=1X2=10D6=1D1=4X1=1X2=2X3=1X1=4X1=1X2=2X1=3X1=1X1=1X1=3D2=1X3=1X2=1X1=1X2=5X2=2X1=1X1=3X1=5X4=1X1=1X2=3D2=3X3=4X1=1X2=2X1=2X1=3X3=3X2=2X1=1X1=1X4=1X1=1X1=2X6=4X1=2X1=1X1=2X2=1I3X3=3X3=2X3=2I2=1X4=1I1=2X5=3X1=1I1X1=2X1=1X2=2X1=1X2=1X1=2X1=1X1=1X1=3X3=1X3=1X1=3D3=4D4X5=2D1=1X1=5X2=1X2=3X2=1X1=1X3=2I5=1X1=1I2X2=1X1=3X1=1X1=1X1=3X2=1X2=5D1=1X2=3D1X3=1X3=6X4=3X1=1X1=3X2=1X1=1X3=1X1=2X3=1X1=1X2=2D1=1X1=1X1=2X1=3X2=2X2=3X1=1X1=1X2=5I2=1D1=1X5=1X3=1X3=3I3=1I1=1X2=1X1=3X1=3X1=1X3=1I2=1X2=1I2=3X2=1X1=1X2=2I1X1=1X1=1X1=2X2=7X2=2X2=1X2=4X3=7X7=4D2=1X1=1X2=2X2=1X1=1I1X1=1X2=2X3=4I1=2X2=1X3=2I5=3X2=2I3=1X1=1X1=2D3X2=3X1=1X1=1X2=2X3=2X2=1X1=4X3=2D1=1X1=1X2=4X2=1X2=1X2=1D1X1=1X2=1X4=2D1X3=4X1=1X1=2X2=2X1=2X2=1X1=2X2=1X1=2X2=1X2=2D3=1X2=2X1=1X2=1X1=1D3X4=1X2=2I1=4X2=4X3=1X4=8D1X1=2X2=2X2=2X7=1X1=1X1=1D1=1X1=1X2=5X6=8D1=1X3=2D1X4=2D3=1X1=4D2=4X3=4X1=2X1=1X1=1X1=1X1=1X3=2X1=3X1=4X1=1X1=1X3=4X2=7X1=1X3=2X1=1X2=1X2=1X2=3X1=3X3=1X1=2I1=1X1=2X1=5X1=1X1=1X1=4X4=1I4=3X1=1X2=1X1=3D1=2X2=1X3=2X1=1X2=1X1=2I1=2X1=1X4=2X1=2I1X2=2X2=1X2=2X2=1X2=2I1=1X1=1X1=1X3=1X1=1X1=4I5=1X1=1I2X1=2X4=1D1=2X1=1X1=1X1=1X1=3X4=1X1=1X1=1X1=2I4=1I2=3X3=4X2=1X1=4D2=1X1=2X1=2X1=1X2=2X2=3X4=1X1=8D3=1X1=3X2=3X2=1X4=1X2=2X2=3X1=2X1=1X2=1X2=3X1=1X2=1X1=1X1=1X1=2I3=1X1=2X3=4D2X4=4X4=1X2=3X2=3X1=2X1=2X2=1D2=1X3=1X1=2D2=1X1=2X1=6D1=2X1=2X2=3X5=2X3=1X1=6D2X2=2X1=2X5=1X1=1X2=4I1X3=2X1=3I3=4X1=1X2=2X2=1X2=1X1=1X1=4X1=2X2=1X3=1X1=1X2=5D1=1X1=1X3=2X1=3X2=2X4=3I4=1X2=1X1=4I1=2X3=1X1=1X1=6I6=13I2=1X1=1X1=1X1=1X1=2X2=2X3=1X1=1X1=5I1X1=1X2=1X3=1X1=2X1=1X1=3X3=1X3=1X1=2X2=4I1X2=1X6=5I1=1X1=1X4=1X1=1X1=3X1=1X1=3X3=3I2X2=2X3=2D2X1=1X2=7X2=3X5=9I1=1X4=1X1=1X4=1X5=1X1=4D2=2X1=2X1=2X2=1X1=1X1=1X1=3X1=2X1=1X1=1X1=4X1=1X1=1X2=2X1=2I6=3X2=1X2=2X2=3I1=1X2=4X1=1X1=1X1=2X1=1X1=1X2=1I1=1X2=2X1=4I4=4X1=2X1=1X1=2X2=2X2=4X1=1X1=2X1=1X1=1X1=2X1=3X2=2X5=1X1=1X1=1X2=1I1X2=1X1=1X1=3I1=2X2=2I1X2=2X1=2X2=2X2=1X1=1X1=1X1=7X1=1X1=2X2=1X1=1X3=5D4=1X1=3X2=1X2=2X3=3X2=2X2=1X3=2I3=1X1=1X1=1X3=2I1X1=1X1=2X1=3X1=1X4=2D2=1X1=2X3=6D1=1X1=1X1=1X3=1X1=3X1=5X2=3X1=1X3=2D6=1X1=2X3=1X2=1X1=4X1=1X1=1X2=2X1=2X3=2X1=2X2=5D3=1X3=4X1=2X1=2X1=3X2=2X1=5X2=1X2=3X1=4I1=1X2=2I3=1I4=1X1=1X1=1I3X2=2X2=1X2=3X1=3X1=3X2=2X3=1X2=3X3=4I1=1X4=1X1=1D2X1=1X1=2X1=3X2=1X1=2X2=2X1=1X3=2X4=1X2=3X1=1I1=1X2=1X2=1X2=2I2X1=1X3=1X1=4D1X2=4X2=4X5=1X2=2I2=5X3=1X2=1X2=3X3=1I2=1X2=1X1=3I2X3=2X2=3X1=1X1=2D5=3D2=2X1=2X2=1X1=1X1=1X1=9X2=1X1=3X5=2X1=1X3=2I1=2X3=1X1=1X1=2X1=1X2=3X2=2X1=1X2=2X4=10I1=1X2=1X2=1I1=1X1=1X1=2I3=3X1=1X2=2X1=1X1=2D1=3X3=1X2=1X2=3I4=3I4=2I5=2X4=1I1X5=6I1=2X3=1X1=4I1=1X1=2X2=1X2=1X1=2I2=1I2X3=1X1=2D4=1X1=1X1=1X2=2D2X4=4X1=1X1=1X1=2X1=1X2=6X2=1X1=2D2=1X1=1X3=4X2=2X1=2X1=1X1=2X1=2X1=1X3=2D2X2=1X1=1X1=1X2=1X2=2X2=3I6=1I1X1=4X1=1X3=2X2=4X3=5X2=6X2=3X3=1X3=2D2=4X1=1X1=1X1=2X1=1X1=2X2=1X3=1D2=1X1=1X1=2I3=1X2=1D3=1X2=1D3=3X2=1X1=1X1=1D1X1=3X1=1X1=1X2=2I1X1=1X1=1X1=1X2=2I1=1X4=1D2=6I1X5=2X2=1X2=2X2=1X3=3X1=3D1X3=1X1=3X1=2X1=4X1=2X4=1X5=2X3=1X1=2X1=1X2=1X1=1X1=4X1=1I1=1X1=1X1=3X3=4X2=2I2=1X2=1X1=4X5=1D2=2X1=2X1=1X3=4X1=2X4=3D1=2X1=1X3=1X2=1X1=2X2=1X1=2X1=1X1=2X1=1X3=1D2=1X1=3X1=1X2=2X1=5X1=6X2=1X2=1X2=6D2X1=1X3=1X1=1X1=1X2=1X1=3D4=1I2=2X2=3X1=1X1=1X1=1X1=1D1X3=1X3=3X2=1D1=1X3=3X1=1I2=2X1=1X1=4I2=1X2=2X2=2X1=2X1=1I3X1=2X4=2X1=2X1=2X3=1I1X1=2X1=1X3=4X1=1X4=5I1=1X1=1X1=1X3=2X3=1X2=3I1X2=3X2=2X1=1X1=1X2=1X2=2I4=6I5=2X2=1X2=1X1=1X1=1X1=3I1X1=1X3=1X2=9X1=1X2=2X1=1X2=2X1=2I5=1X1=3D1X2=2X3=2X3=1X2=2I1=1X3=1X2=3X1=1X1=1D4=2X1=1X2=4X1=1X1=3X1=3X1=3X2=3X1=1X1=2X2=4I3=1X1=1X1=1I2=1X2=2D1X2=2X2=1D1X1=2X1=1X1=1X2=1X1=4D2=1X5=1D1=1X1=4X5=2X5=3X3=2I1=2X1=1X1=1X1=3X1=6X4=4D1=1X2=2X1=1X2=1X1=1I1X1=1X1=1X1=1I2X2=1X1=1X4=1I3=3X1=3I2=1I1=1X2=1X1=1X2=2X4=2D1X3=2X1=2X2=6D1=1X2=2X2=2X5=1I1=1X1=1X1=4X5=3I2=1X2=2X1=4X1=5X1=2X1=1X1=2X1=2X1=2I2X1=1X1=1X1=1X2=1X2=1D1X6=2X2=1X2=2X1=3X1=1X7=1D3X1=3X1=2X2=1X1=1X2=2X1=1X1=2D1X2=1X3=1X1=4D1=2X4=1X1=2X1=2X1=3X1=5X2=3X1=1I1=3X1=3X1=1X4=1X3=5I1X2=1X1=1X4=4X2=2X2=1X1=1D1X2=1X1=1X1=2X1=1D1X1=1X2=3I2=2X1=6X1=1X2=1X2=1I3=2X1=3X1=1X1=1D4=2D1X1=3X2=3X3=2X2=5D3=3X2=1X1=2X1=2X1=1I3=2X2=1X1=2D1=1X1=1X1=3X2=2X2=2X1=3I1X1=1X1=3X1=1X1=1X3=3I2=3X2=1X1=2X1=1X3=3X1=1X1=3X1=2D1=4X2=1I3=2X2=1X1=1X1=3X1=1X2=1I1=2X1=2X2=2X1=1X2=5I1X1=1X2=3X4=3X1=1X1=3X2=1X1=1X2=1I2X3=1X3=1X1=1D1=1X1=2X4=1X2=1X1=5D2X1=1X2=1X1=2X4=3D1=2X2=2X2=1X1=1X3=1X2=3X1=1X1=1X1=2X1=1X2=2X4=1X1=7I1=1X3=1X3=1X2=2I1X3=2X1=1X1=1X2=1X3=3I2=2X1=2X1=2X1=1X1=1X3=4X1=2X2=3D2=1X1=3X1=1X1=3X3=1D2=1X2=1X1=1I2X1=2X1=2X1=4X5=5X1=1X2=3X1=2X1=1X1=2X3=1X3=5X2=2I1X1=1X2=1X4=2X1=3D2=2X2=1X1=2D1X2=1X2=1D2=2X4=1X3=2X1=1X2=6X4=1X1=3D1X1=3X1=1X1=1X3=3X3=1X1=1X1=1X2=2D2X1=2X2=2X2=2X2=2X4=3I1X1=1X4=1D4=3I3X6=4D2X3=2X3=2X1=1X1=1X3=4D2=1X1=1X1=1I2=2X1=4X1=1X2=1X2=4X1=4X1=2X1=1X4=3D4=1D1X1=2X2=2X1=2X2=2X2=1X2=5D4=1X1=1X2=6X1=2X3=1X1=1X2=2I2=1D1=1X4=4X2=2X1=1X1=3X2=4X1=1X1=1X1=1X1=2X1=1X1=1X2=2D1=1X3=3X1=3X1=2X2=1X1=1I2=1X2=1X1=2I1=1X3=1X1=7I2=1X2=2X1=2X1=3I6=2X2=4X3=1X1=1X1=6I3=2X2=1I1=2X3=1X1=4I3=2X4=1X2=1D1X3=1X1=2X3=5I1=1X3=1X1=1X5=1D3=1D1X2=1X1=1X3=1X1=1X1=1I1=1X2=2X1=1I4=1X1=1X1=2D2=2D2=2X1=1D1=2X1=1X1=1X7=6I1=1X3=1I3=1I1X2=1X2=3X3=1X2=2X1=2X2=8D2=1X2=1X1=1D2X1=1X1=2X4=2X2=2X1=3X3=3D1X2=1X1=2X2=5D2=1X4=3X1=1X3=5D1X1=1X2=4X3=2D1X2=1X2=3X2=1X1=1X1=1X1=2X1=1X1=1X3=2D4=2X2=2X1=3X3=1X1=1X1=1I2=6X3=1X2=1D1=1X1=5X1=1X3=1X1=1X1=4X1=4X3=2X2=1X1=1X3=1X2=3X1=3X4=2X3=5D2X1=2X1=1X2=1X1=2X2=1X2=1I1=1X1=2X3=1X2=2X1=1X3=2X2=1D1=1X1=2X1=2X2=3X1=1X2=1X1=1X3=1X1=10I1=2X4=1X1=3X2=5X1=2X2=1D1=1X1=1X1=2X1=1X1=1X2=3X2=1X3=1X2=2X1=4X2=1D7=4I1=2X1=1I1=1X2=1X2=1X1=1I1X1=3X1=1X3=2X2=3X2=2X2=2X3=2X2=2X1=2X1=1X2=1X1=2I2X4=1X1=1I3=3X3=3X1=5X1=2X3=1X2=1X4=1D1=2X1=2X2=2X1=5X4=1X1=1I1X1=1X2=1I4=9I1=1X5=3X1=1X3=1X1=3X1=1X4=3I4X1=1X1=1X5=1X1=1X3=1D1=1X1=2X3=2X2=1X1=1I1X3=1X1=3X5=1X1=1X1=1X1=1X3=5D2=1X3=2X1=1X2=1D3X2=4X4=1X1=1X1=5X3=1X1=3D2X2=1X3=1X3=1I1=3X4=3D3=2X1=1X1=2D1=2X1=1X1=1X1=1X2=1X1=3X1=2X1=4X1=2X1=1X1=2D3=1X1=2X3=1X4=1X1=1I2=3X2=2X1=1X1=2X1=3X4=1X2=2X1=3X1=1X2=3X2=1X1=1X1=6X3=3D5=2X3=4X3=7D2=1X1=2X2=3X4=1X1=1X3=2X1=2I2X5=1I1=1X2=1X2=1X1=3D3X3=6X2=2X3=2X1=2X1=2X1=1X3=1X1=3I1X2=4X3=2X4=1X2=5D5=1X1=3D2X3=2X2=1X1=2D3=1X2=1X2=3D4=2D2X1=1X8=5D1X3=2X2=2I3=1X2=4I1X3=1X1=1X2=2X1=1X4=4X2=1D2=2X1=4X1=1X1=1X2=1X1=2D2X1=3X4=1D2=2X1=2D5=6X2=1X1=3X1=2X2=2X1=1X2=5X1=1X2=3X1=2X1=2X2=1D3X5=2X1=3X1=2X4=2X1=2X1=2X3=5D2=3X1=1X3=2X2=2X1=2X1=1X1=2X2=2X3=1X1=2X1=2X1=1X1=1X1=3I3=5X2=3X3=2X1=1X4=3I1X2=7X3=1X3=1X2=2D1=1X2=1X2=2I1X1=2X3=1X5=1X2=3I4=1X2=1X3=1X2=8I1=2X1=2X1=1X5=4I8=3I5=9I5=2X1=1X1=1D1=1X1=1X2=1X1=6I1=1X3=1D3=1X3=1X1=3X1=1I2X1=4X2=1X2=2X1=2X3=2X2=1X1=2X2=1X3=1X1=1X2=1I1=1X1=1X4=5D2X2=2X1=3X1=1X1=1X2=4D3=1X2=2X1=3X1=1X1=1X2=4X4=5D3=1X3=2X1=1D2X2=1X1=7X1=1X1=5X1=2X7=3D1=1X2=1X2=3D3=6D2=1X4=6X2=2X1=4X1=2X4=1X3=3X1=4X1=2X3=3X1=1X1=1X2=2X4=1D1=1X3=1D2=1X1=2X3=4X2=3I1=1X2=4X2=3X2=3X1=3X1=1X4=1D6=2X3=1X1=2X3=2D1X2=1X2=3X1=1D1X2=1X1=2X3=5D2=1X2=1X3=1I1X2=5X2=2X1=1X2=1I2=1X1=2X1=1X1=2X1=1X2=2X2=1X1=1X1=2X1=2X1=1X2=1X1=1X2=3I1X1=3X1=1X2=1X2=1X4=1X1=3I2=1X2=1X5=2X1=7X2=1X3=3D1=1X2=1X1=2X2=1D2X3=1D1=2X1=1X2=3X2=2X2=2X1=1X1=10X3=2X1=2X1=5X4=2I2=5X1=1X4=3X1=1X5=4D1X1=1X2=1X2=1X2=2I2=1X1=1X1=1X1=1I1=1X2=4X3=1D3=1X3=7I2=1X1=1I1X3=5X2=2X2=1D2=3X2=1X2=2X1=7X2=1I1X8=1I3X1=2X5=3X1=1X3=2X1=1X3=2I5=2X1=3X3=1X3=1X1=3X4=2X1=9X1=1D1=1X1=2X1=5D6X1=1X4=1X1=2D1=1X5=4X2=2X1=1I2X2=2X1=4X1=6X3=1X2=2X2=3X5=1D4=3I1=3X4=1X1=1I3X4=1X2=2D3=2X3=7D1=3X3=12D3=3X2=1X2=2D5=7X1=1X6=2X2=5D1=3X2=1X2=1X1=2I2=1X2=2X2=3D2=2X1=1X2=2D1X1=1X1=3X2=2X5=1X1=3X1=12D1=1X1=1X4=1X1=1I5=4I1=1X4=2X2=1X1=1X2=2X1=1D2=1X3=1X1=1X2=1I1X1=3X1=8X2=1X2=7D2=2X2=1X3=2D1X1=1X1=1X1=1X1=2X1=4D2X2=1X3=1X1=4X1=3X1=1X4=2X4=1I1=2X2=2X1=2X2=2D5=5D1X2=5X4=3X1=3X3=2X2=1I1=2X2=1X1=4X2=1X4=2X1=1I2=1X2=1X1=1X1=4X1=1X4=7X1=2X1=2X1=3X4=5D1X3=1X2=1I1=1X1=1X1=2X1=2X1=1X1=1X1=1X1=3X1=2I1X1=2X1=1X2=1X3=1X4=2I1X2=2X1=3X4=1I1=1X3=1X1=2X1=1X1=3X2=2I1=5X1=3X1=1X2=2X1=1D1=1X4=3X1=2X1=1X1=5X2=1X1=2X1=2X1=3X1=2X1=1X1=2D1X2=1X2=1D3=1X1=3D3=2X2=1X1=4X1=2X2=2X1=2D1=1X3=2X1=2X2=2D6=2X1=3D3=1X1=2X2=2X3=2X2=1X2=2D2X3=2D1X1=4X1=2X1=1X2=4X1=1X1=1X1=2X2=1X2=4D2=1X7=2X1=4X2=1X2=3D1X4=1X2=1X1=2X2=2X4=1I3X1=1X3=1X1=1X1=1X3=2X3=3X1=7I2=2X3=1I1=1X1=1X1=1X3=2X2=1X1=2X1=2X1=1X1=1X2=4X2=1X1=2X4=3I2=1X2=2X1=2X1=1X2=2X3=2D1X1=5X1=1X1=1X1=3X1=1X2=1X2=3I1X1=1X4=1X1=1X1=1X5=1X1=2D1=1X2=1X2=1D2=2X2=4X2=2X1=2D2=1X2=3X1=2I3=3I2=1X6=4I1=2X1=1X2=1I1=1X1=4X3=2X2=4X1=2X2=1X1=1X1=1D4=2X3=1X1=1X1=3X2=1X1=4D1X3=2D2=1X1=3X1=2X1=1X2=1X2=5X1=1I4=1X1=2X1=1D4=1X1=3X2=5D1=1X6=1X2=5X1=1X1=1I1X2=1X4=1X1=1X4=2D1=1X2=2D1X1=1X2=3X3=1X1=4D1X1=2X5=1X1=1X1=2D2=1X1=1X2=4D4=3X2=2D1=1X3=1D2X3=1X3=3X2=2D3=1D1=1X1=2D1=1X2=3X1=2X3=1D1X2=1X4=1X1=2D1X1=1X2=1X1=1I2=1X3=2X2=1D2X2=4X2=1X2=1X2=1D1=2X2=2X1=1X1=3X4=3X1=3X1=1X1=1X1=1X1=1I4=2I1=1X3=2X1=3X1=3X4=1X1=3X1=1X2=1X2=1X2=2X1=6I3=1X1=1X1=1X1=3X1=1I6=2X1=1D4X1=2X1=1X1=1X3=2X1=3D1X8=2X2=1X2=2D1=3X1=1X3=1X2=1X2=1D1=1X3=3X1=2D4=1X3=1X2=4I1X2=1X2=1X1=2X1=6X1=1D1X1=1X2=1X4=1D1=1X1=1X3=1X1=1D1X4=5X1=1X3=1X3=5I3X1=3X2=1X2=5X1=2X3=1X1=1X1=3X3=2X1=3X1=1X1=3X1=1X1=2X3=2X1=2X2=1X1=2I1=1X4=1I1X1=1X4=1I3=1X1=1I5=2X3=4X1=1X1=1X1=3I1=2X1=1X2=3X2=1X1=2X3=1I1=2X2=1X1=1X1=1X2=1X2=1X2=4X4=1X1=1X2=2X1=1X2=1X1=1X2=4X3=2D1=4X1=1X5=2X2=9I3=2X4=3I2=1X1=2I1=1X3=1X3=2X3=1X1=1I3=1X3=2X5=1X1=9D1X2=3X1=1X4=1D1X1=1X1=1X3=5X2=1X1=2X1=1X3=1X1=6D4=3X1=1X3=1X1=1X1=1X1=1X1=2X1=4D1=1X4=2X1=1X1=1X1=2X1=1X1=2D1=1X4=1X1=1X2=4X1=2X1=1X3=2X1=1X3=8I2X1=2X3=3I1X2=2X4=1I2X3=1X1=5X2=1D1X2=1X2=4X2=2X1=4X1=1X1=1X4=3I2X1=1X1=1X1=2X1=4I1X5=1X1=2X1=1X2=2X1=1X1=2X1=2D1=1X1=1X3=2X1=2X1=1X1=1X2=2I1=2X1=1X1=1X1=1X4=2X1=2X1=3X4=2X2=6D1X3=1X1=3X2=3D1X1=2X3=1X1=1X1=6X2=1X2=2X1=1X2=2X1=2X2=2X2=2X1=1X1=4X4=2I1=2X1=1X2=2X2=1X1=3X2=2I2=1X3=1X1=2I1X1=1X1=1X1=1X1=1X3=1X5=5I1=1X4=2X1=2X2=3X3=1D2=1X1=3X1=3X1=4X1=3X4=5I2=1X2=1X2=1X1=1X2=1D4=1X2=1X4=2X1=1X1=3I2=1X1=1X3=1I1X1=4X1=1X1=3X1=3X2=1X2=3X1=2D1X1=3X7=2D2=1X1=1D2X1=1X2=1X2=2X6=1X2=5I1=1X2=2X1=2X3=2X6=3D2X2=4X2=1X1=1X2=1X1=1X3=2D1X2=2X1=1X2=2X2=3I3X1=2X3=1X1=1X1=1D3=1X2=1D1=1X1=1X1=2X1=2X2=1X4=4X1=1X2=3X1=1X3=5D1X6=1D1X1=1X1=3X2=2X2=1X4=2X1=3I1X1=3X3=1X2=2X2=1X1=2X1=1X1=2X2=2X2=2X2=2D1=1X1=1X1=1X1=5D1X2=2X2=3X1=2X2=2X1=1X4=7D4=2X1=2X2=1X2=1X1=5D1=2X1=2X1=3X2=3X1=1X1=4D2=1X3=1X1=1X1=1X1=3X1=1X1=1I3=3X1=2X1=2X2=6X1=1X5=1X1=3D1X2=1X2=1X2=6I1=1X1=1X3=2X3=2X2=1X1=2I1X1=2X5=1I3=2I2=1X3=3X2=1I1X2=2X4=2I2=1X1=1I4=1I2X3=2X1=1X1=3X1=1X2=4X1=2D2=1X1=1I1X2=1X1=1X1=1X1=1X3=1I2=1X2=3D5=1D1=2X2=2X2=4X1=5X2=1X2=7D1=2X5=2X1=1X1=1X5=1X1=3D1=1X1=1X3=3D1=2X3=1X1=2X4=1X1=4D1X2=1X1=3X2=1X3=4X2=1X2=1X1=6I2=2X1=1X3=4I1X1=2X1=3X1=2X1=1X1=1X1=1X1=1X7=2X3=3X2=5X1=2X1=2X2=2X2=2X2=2X2=11D3=4D3=1X2=1X4=8X1=1X1=3D1=2X2=1X1=2X2=1X1=1X1=1X1=1X1=1X1=1X1=3X2=1I1=1X2=2D2=1X1=5X6=1D1X1=3X3=2X4=2D4=1X2=1X2=2X2=9D1=1X1=3X1=1X4=1D2X1=1X3=1X2=1X1=2X2=2X2=1D2=1X1=2X2=3D1X1=2X1=1X2=2X3=1X2=2X4=1X4=2X1=2X1=1X2=4X2=2X1=2X2=2X1=2X1=1X1=1X1=1X1=2X1=3I2X1=5X3=1X5=7X1=1X1=2X1=1X2=1X1=1X1=1X3=2I2X2=1X1=3X3=3I2=1I4=2X1=1X3=1X1=1X4=2X1=3X1=1X1=2X1=2X2=5X1=1X1=1X1=4X1=1X3=1X1=1X1=1X1=8I1=1X1=1X1=3X6=4X1=1X2=1X3=2X1=3X2=1X1=8X3=3X1=1X1=1X3=4I1=3X3=2X2=1X1=2I2X1=1X2=2X3=1I1=1X2=1X1=1X1=1D3X3=2X1=1X2=5X1=1X1=1X1=1X1=1I2=2X2=2X2=3X2=1X1=1X1=1X1=5X1=1X2=1X1=1X2=3X2=1X2=4X3=1I2X5=4X3=2X3=1X1=1X2=1X1=5D1=1X4=1X2=1X2=1I2=1X1=2X1=2D1=1X1=2X2=4X1=2X1=1X1=6X4=5D3=1D5=2X1=1X2=2X1=1X1=2X2=1D1X3=1X1=3X2=1X1=2X3=1D1X1=1X1=2X1=2X3=3X3=4I1X1=1X1=1X3=1I1=4X3=1X1=3X4=2X1=2X1=3X1=1X1=1X1=1X3=1I1X3=7X2=2X1=1X2=1D2=1X5=1X1=9I1=1X1=1X1=3X6=1X2=1X4=1X1=1X2=4D1X1=1X4=7X1=1X3=1X1=3X3=1D2=1X1=2I3X1=1X2=1X1=9X1=1X1=1X1=4X1=2X2=1X1=3X1=1X1=4X1=1X3=2X3=1X2=5I3=2X4=2D4=3X3=1X3=2X4=1X1=4X1=2X1=2X3=8X1=3X1=2X2=1X2=1D1X2=1X1=1X1=1X4=2X1=3X2=1X1=1X3=3D1=1X4=2X2=4D4=2D1X1=1X3=3X3=1I3=6X1=1X3=1I3=5X2=1X1=2X3=1I2=1X2=1I1=1X2=1X1=1X1=1X2=1X1=1X1=2D1X1=1X1=5X2=2X2=1I2=2X2=1I7=6I1X2=1X4=1I1=1X1=2X1=1X1=4X3=1I1X1=1X1=1X1=5X3=1X1=3D1=3X3=3X1=2X2=9X9=2X1=1X1=2X1=1X3=3D1X1=1X2=3X1=1X4=3D7=7D1=1X1=1X2=1X1=1X3=1X1=1X1=1X1=2X2=1X1=1X3=6D2=1X1=1X2=3D3=1X1=2X3=2X2=1X1=1X3=1D1=1X1=5X1=1X1=1X1=1X1=3X1=1I1=1X1=3X4=2X1=1X1=1X1=1X2=1X1=1I2X2=2X1=2X1=1X1=2I1X1=4X3=1X2=2X2=2X1=3X2=9I1=1X1=1X3=2I2=3X7=6D2=1X1=4X1=1X1=1X3=2X2=1I1=1X6=1I2=2X3=2I6=2X3=6X2=4X2=1X1=1X2=1X1=1D1=3X1=1X3=1X1=3X1=1X2=4X2=2I2X3=1X1=1X1=1X2=2I1=3X6=4I1=1X1=4X1=3X2=2X3=1X1=1X3=1X3=2X1=1I2=3I2X5=1X1=2X1=2X1=1X1=3X3=2X1=17I4=1X1=1X1=1X2=1I1=4X5=1I1X1=1X1=1X3=1I1X1=2X3=3X1=1X1=1X2=4X1=1X1=1X2=1X1=1I1=1X2=1X2=4I1X2=2X1=1X3=1X2=2X2=9I1X3=5I2=2X1=1X6=1X1=2X2=1X1=1I1=2X3=1X2=2X1=1X1=4X1=1X1=1X1=1X3=1X3=1I1=1X1=1X1=2X4=2X1=4X1=3X1=2X3=5D6=1X1=2X3=1I6=2X1=1X4=2X1=1X2=2D1=1X6=1X1=2X1=1X1=3X3=3X2=4I1X3=2X2=1I4=1X1=3D1=2X2=1X4=2I2X1=1X1=1X1=1X1=1X1=6X1=1D1=2X2=1X1=4X2=3X3=1X5=1X1=2D3=1X1=1X1=1X3=2X1=3I1X1=1X3=3X2=1X1=3I2X4=1D1X2=2X3=1X1=3X1=1X1=3X3=1I1=2X1=1X1=1X1=1X1=6I1=1X1=2X2=2X1=2I3=1X4=1I2X1=1X3=1X1=4X2=1X1=2X2=9X3=1I1X5=1X1=2X6=1X1=2X1=5D1X4=2X4=3X1=3X1=2X3=2X1=2X4=4X2=1X1=3X3=1I3=1X3=2D2=2X2=2D2=1X2=2X2=3X1=1X1=1D2X3=2X1=5X2=2X2=1X1=1X1=2X1=2X1=2I2=1X1=1X2=4X1=5I1=1X1=1X2=1X5=4D1X6=1I1=1X1=1X2=3X1=1X2=1X1=3X1=4X2=1X1=1X1=1X1=2X2=1X2=2I2X6=5I3=6X1=1X2=2X3=1I3=1X4=5I3=1X3=1X4=5I1=4X6=2X1=3I3X3=1X1=4X1=2X1=3X1=1X2=3X2=1I5=3D3=3X1=1X1=3X1=1X5=1X3=2X3=2X1=1X2=3D1=2X4=3X2=2X1=1D5=1X2=1X1=2D2X3=1X3=1X2=3D1X1=2X1=2X2=1X1=1X1=1X1=2X3=1X3=1I3X1=1X3=1D2=2X1=1X1=1X2=1X2=5I1=1X2=1I1=4X3=1X1=1X3=1X1=1D1=1X3=3X4=2X1=2X2=4X3=1I2=2X1=1X2=3D1=3X1=3X2=1X4=3I4=2X1=2X2=1I3X4=3X3=3D2=2X4=1X1=2X1=1X1=1X4=1X2=1X1=3D4=1X1=2D3=1X3=1X1=1X1=1X1=2X1=1X1=3D1=1X4=2X2=3X2=1X1=3X1=2X2=2I3=1X1=1X1=1D1=1X2=2X2=3X1=2I1=1X1=2X2=1X2=3I1=1X4=1X1=1X1=2X1=1X1=1X1=7D1=1X5=1I1=1X4=2X1=1X2=1X1=1X1=5D4=1X3=1X2=1X1=1I1=1X1=1X4=8D3=2X2=8D1=1X2=3X5=8D1X1=1X3=3X1=1X1=2X1=1X1=4D1X1=3X1=1X3=3X1=1X2=2X2=1X2=1D2=2X3=1X2=1X1=1X3=2D4=1X1=2X1=3D1X2=1X7=3X3=1X1=1X2=2I4=2I1X2=1X4=3X1=1I3=1X1=1X1=4X1=2X2=1X2=4X1=3X2=1I2=4X1=2X2=1D2=1X1=1X1=3X5=2X1=3I5=1D1X6=1X4=9I2=1X1=1X1=2X1=1X3=1X2=1X1=1X3=1X1=2X1=6D5=2X1=2X1=1D3=2X1=1X2=1X1=4D3=2X1=1X2=1X1=3X1=1X1=1X2=3X1=1X2=2X1=4I3=1X3=2X1=3X1=1X1=1I2=1X1=3X1=2X1=2X5=1X1=1I1X2=1I5=3X3=1X2=1I2X3=1X1=6I1=2X77=1X3=1X5=1D14=1X28=1X201=1X41=1X26=1X125=1X59=1X26=1X98=1X11=1X5=1X5=1X23=1X11=1X8=1X8=1X35=1X26=1X176=1X5=1X104=1X24=1X43=1X147=1D52=1X21=1X46=1X39=1X39=1X384=1X1=2X3=2X3=1X7=2X1=4X4=2X2=1X33=1X17=1X79=1X66=1X20=1X53=1X4=1X18=1X11=1X32=1X38=1X21=1X7=1X8=1X2=1X5=1X5=1X2=4X5=1X5=3X3=1X2=1X3=2X12=1X1=2X2=1X2=1X14=1X29=2X1=1X7=1X9=1X5=1X2=1X8=1X1=3X4=1X8=1X2=1X14=1X8=1X30=1X16=1X1=1X12=1X1=2X5=1X2=1X8=1X26=1X11=2X6=1X3=1X6=2X18=1X102=1X76=1X158=2I43=1X34=1X23=1X9=1X23=1X2=1X67=2X1=2X59=1X14=3I36=1I13=1X36=1X41=1X290=1X133=1X53=1X4=1X58=1X52=1X29=1X30=1X47=1X62=1X47=1X172=1X75=1X17=1X12=12D15=1X47=1X32=1X119=1X5=1X19=1X161=1X1=1X138=1X4=1X42=1X9=2X118=1X50=1X42=1X50=1X18=1X25=1X4=1X28=1X1=1X164=1X81=1X21=1X25=1X21=1X8=1X21=1X11=1X20=1X6=1X5=1X22=1X76=1X42=1X140=1X92=1X20=1X7=1X3=1X62=1X70=1X30=1X31=1X80=1X118=1X14=1X145=1X53=1X137=1X76=1X97=1X30=1X2=1X29=1X341=1X140=1X8=1X755=1X197=1X185=1X876=1X5=1X244=2X110=1X383=1X347=1X40=1X121=1X17=1X119=1X24=1X22=1X41=1X76=1X209=1X40=2X18=1X20=1X8=1X55=1X89=1X54=1X24=1X145=1X26=1X409=1X42=1X6=1X28=1X3=1X7=1X60=1X158=1X43=1I11=1X46=1X67=1X134=1X20=1X4=6D1X16=1X264=1X152=1X83=1X38=1X23=1X236=1X11=1X437=1X8=1X2=1X29=1X17=1X80=1X409=1X39=1X116=1X111=1X28=1X71=1X166=1X18=1X173=1X145=1X149=1X108=1X50=1X7=1I113=1X191=1X440=1X357=1X88=4D48=1X805=1X104=1X471=1X134=1X17=1X16=1X12=1X51=1X12=1X4=1X369=1X12=1X393=1X139=1X157=1X12=1X29=1X23=1X475=1X420=1X53=1X14=1X189=1X178=1X458=1X178=1X470=1X114=1X65=1X83=1X156=1X43=1X198=1X218=1X79=1X405=1X452=1X72=1X153=1X171=1X38=1X87=1X217=1X406=1X687=1X316=2D17=1D84=1X214=1X51=1X182=1X6=1X181=1X50=1X92=1X116=1X272=2D1489=1X704=1X21=1X24=1X48=1X14=1X74=1I38=1X29=1X74=1X420=1X424=1X121=1X297=1X346=1X35=1X17=2X175=1X25=1X81=1X180=1X272=1X75=1X129=1X98=1X9=1X736=1X215=1X16=");
                println!("from_paf {:?}", a);

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
            error!("{}", message_with_truncate_paf_file("Skipping malformed PAF line", line));
            std::process::exit(1);
        }
        
        let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
            error!("{}", message_with_truncate_paf_file("Skipping CIGAR-less PAF line", line));
            std::process::exit(1);
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
            error!("{}", message_with_truncate_paf_file("Skipping malformed PAF line", line));
            std::process::exit(1);
        }
        
        let Some(tp_field) = fields.iter().find(|&&s| s.starts_with("tp:Z:")) else {
            error!("{}", message_with_truncate_paf_file("Skipping tracepoints-less PAF line", line));
            std::process::exit(1);
        };
        let tracepoints_str = &tp_field[5..]; // Direct slice instead of strip_prefix("tp:Z:")
        
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

        // Create a thread-local FASTA reader
        let fasta_reader = FastaReader::from_path(fasta_path).unwrap_or_else(|e| {
            error!("Failed to create FASTA reader: {}", e);
            std::process::exit(1);
        });
        
        // Fetch query sequence
        let query_seq = if strand == "+" {
            match fasta_reader.fetch_seq(query_name, query_start, query_end - 1) {
                Ok(seq) => {
                    let mut seq_vec = seq.to_vec();
                    seq_vec.iter_mut().for_each(|byte| *byte = byte.to_ascii_uppercase());
                    seq_vec
                },
                Err(e) => {
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
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
                    error!("Failed to fetch query sequence: {}", e);
                    std::process::exit(1);
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
                error!("Failed to fetch target sequence: {}", e);
                std::process::exit(1);
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
