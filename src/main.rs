mod sequence_index;

use crate::sequence_index::{collect_sequence_paths, SequenceIndex};
use clap::{Parser, ValueEnum};
use flate2::read::MultiGzDecoder;
#[cfg(debug_assertions)]
use indicatif::ProgressBar;
#[cfg(debug_assertions)]
use indicatif::ProgressStyle;
#[cfg(debug_assertions)]
use lib_tracepoints::{
    align_sequences_wfa, cigar_ops_to_cigar_string, cigar_to_tracepoints_raw,
    cigar_to_variable_tracepoints_raw,
};
use lib_tracepoints::{
    cigar_to_mixed_tracepoints, cigar_to_tracepoints, cigar_to_tracepoints_fastga,
    cigar_to_variable_tracepoints, mixed_tracepoints_to_cigar,
    mixed_tracepoints_to_cigar_with_aligner, tracepoints_to_cigar, tracepoints_to_cigar_fastga,
    tracepoints_to_cigar_with_aligner, variable_tracepoints_to_cigar,
    variable_tracepoints_to_cigar_with_aligner, ComplexityMetric, MixedRepresentation,
};
use lib_wfa2::affine_wavefront::Distance;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::sync::Arc;

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

fn parse_complexity_metric(value: &str) -> Result<ComplexityMetric, String> {
    match value {
        "edit-distance" => Ok(ComplexityMetric::EditDistance),
        "diagonal-distance" => Ok(ComplexityMetric::DiagonalDistance),
        _ => Err(format!(
            "invalid complexity metric '{value}', expected 'edit-distance' or 'diagonal-distance'"
        )),
    }
}

fn complexity_metric_to_str(metric: ComplexityMetric) -> &'static str {
    match metric {
        ComplexityMetric::EditDistance => "edit-distance",
        ComplexityMetric::DiagonalDistance => "diagonal-distance",
    }
}

/// Distance model used for WFA re-alignment
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum DistanceChoice {
    /// Edit distance (unit costs)
    Edit,
    /// Gap-affine with a single gap penalty pair
    GapAffine,
    /// Gap-affine with dual gap penalty pairs
    #[value(alias = "gap-affine2p", alias = "gap-affine-2p")]
    GapAffine2p,
}

impl DistanceChoice {
    fn default_penalties(&self) -> Option<&'static str> {
        match self {
            DistanceChoice::Edit => None,
            DistanceChoice::GapAffine => Some("5,8,2"),
            DistanceChoice::GapAffine2p => Some("5,8,2,24,1"),
        }
    }
}

impl fmt::Display for DistanceChoice {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DistanceChoice::Edit => write!(f, "edit"),
            DistanceChoice::GapAffine => write!(f, "gap-affine"),
            DistanceChoice::GapAffine2p => write!(f, "gap-affine-2p"),
        }
    }
}

/// Common options shared between all commands
#[derive(Parser)]
struct CommonOpts {
    /// PAF file for alignments (use "-" to read from standard input)
    #[arg(short = 'p', long = "paf")]
    paf: String,

    /// Number of threads to use (default: 4)
    #[arg(short, long = "threads", default_value_t = 4)]
    threads: usize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[arg(short, long, default_value = "0")]
    verbose: u8,
}

#[derive(Parser)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Encode alignments into tracepoints
    Encode {
        #[clap(flatten)]
        common: CommonOpts,

        /// Tracepoint type
        #[arg(long = "type", default_value_t = TracepointType::Standard)]
        tp_type: TracepointType,

        /// Complexity metric for tracepoint segmentation
        #[arg(
            long = "complexity-metric",
            default_value = "edit-distance",
            value_parser = parse_complexity_metric,
            value_name = "METRIC"
        )]
        complexity_metric: ComplexityMetric,

        /// Maximum complexity value for tracepoint segmentation (default: 32; 100 if type is fastga)
        #[arg(long = "max-complexity")]
        max_complexity: Option<usize>,

        /// Skip adding optional fields (gi/bi/sc fields)
        #[arg(long = "minimal")]
        minimal: bool,
    },
    /// Decode tracepoints back to CIGAR
    Decode {
        #[clap(flatten)]
        common: CommonOpts,

        /// Tracepoint type
        #[arg(long = "type", default_value_t = TracepointType::Standard)]
        tp_type: TracepointType,

        /// Complexity metric for segmentation (must match what was used during compression)
        #[arg(
            long = "complexity-metric",
            default_value = "edit-distance",
            value_parser = parse_complexity_metric,
            value_name = "METRIC"
        )]
        complexity_metric: ComplexityMetric,

        /// FASTA files containing sequences referenced in the PAF (repeatable)
        #[arg(long = "sequence-files", value_name = "FASTA", num_args = 1..)]
        sequence_files: Vec<String>,

        /// File listing FASTA paths (one per line)
        #[arg(long = "sequence-list", value_name = "FILE")]
        sequence_list: Option<String>,

        /// Keep original gi/bi/sc/sc fields as giold/biold/scold when replacing
        #[arg(long = "keep-old-stats")]
        keep_old_stats: bool,

        /// Trace spacing for fastga (only used with fastga type, default: 100)
        #[arg(long)]
        trace_spacing: Option<usize>,

        /// Distance metric for realignment
        #[arg(long = "distance", default_value_t = DistanceChoice::Edit)]
        distance: DistanceChoice,

        /// Gap penalties (only for gap-affine distances; ignored with edit distance)
        #[arg(long)]
        penalties: Option<String>,

        /// Use static band heuristic during decompression (edit-distance metric only)
        #[arg(long)]
        heuristic: bool,

        /// Maximum complexity value (required when enabling heuristics)
        #[arg(long = "max-complexity")]
        max_complexity: Option<usize>,
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

        /// FASTA files containing sequences referenced in the PAF (repeatable)
        #[arg(long = "sequence-files", value_name = "FASTA", num_args = 1..)]
        sequence_files: Vec<String>,

        /// File listing FASTA paths (one per line)
        #[arg(long = "sequence-list", value_name = "FILE")]
        sequence_list: Option<String>,

        /// Gap penalties in the format mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2
        #[arg(long, default_value = "5,8,2,24,1")]
        penalties: String,

        /// Maximum complexity value for tracepoint segmentation
        #[arg(long = "max-complexity", default_value = "32")]
        max_complexity: usize,

        /// Verbosity level (0 = error, 1 = info, 2 = debug)
        #[arg(short, long, default_value = "0")]
        verbose: u8,
    },
}

impl fmt::Debug for CommonOpts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CommonOpts")
            .field("paf", &self.paf)
            .field("threads", &self.threads)
            .field("verbose", &self.verbose)
            .finish()
    }
}

impl fmt::Debug for Args {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Args::Encode {
                common,
                tp_type,
                complexity_metric,
                max_complexity,
                minimal,
            } => f
                .debug_struct("Args::Encode")
                .field("common", common)
                .field("tp_type", tp_type)
                .field(
                    "complexity_metric",
                    &complexity_metric_to_str(*complexity_metric),
                )
                .field("max_complexity", max_complexity)
                .finish(),
            Args::Decode {
                common,
                tp_type,
                complexity_metric,
                sequence_files,
                sequence_list,
                keep_old_stats,
                trace_spacing,
                distance,
                penalties,
                heuristic,
                max_complexity,
            } => f
                .debug_struct("Args::Decode")
                .field("common", common)
                .field("tp_type", tp_type)
                .field(
                    "complexity_metric",
                    &complexity_metric_to_str(*complexity_metric),
                )
                .field("sequence_files", sequence_files)
                .field("sequence_list", sequence_list)
                .field("keep_old_stats", keep_old_stats)
                .field("trace_spacing", trace_spacing)
                .field("distance", distance)
                .field("penalties", penalties)
                .field("heuristic", heuristic)
                .field("max_complexity", max_complexity)
                .finish(),
            #[cfg(debug_assertions)]
            Args::Debug {
                paf,
                threads,
                sequence_files,
                sequence_list,
                penalties,
                max_complexity,
                verbose,
            } => f
                .debug_struct("Args::Debug")
                .field("paf", paf)
                .field("threads", threads)
                .field("sequence_files", sequence_files)
                .field("sequence_list", sequence_list)
                .field("penalties", penalties)
                .field("max_complexity", max_complexity)
                .field("verbose", verbose)
                .finish(),
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments.
    let args = Args::parse();

    match args {
        Args::Encode {
            common,
            tp_type,
            max_complexity,
            complexity_metric,
            minimal,
        } => {
            setup_logger(common.verbose);

            // Determine max_complexity: use provided value, or default to 100 if fastga, else 32
            let is_fastga = matches!(tp_type, TracepointType::Fastga);
            let max_complexity = max_complexity.unwrap_or(if is_fastga { 100 } else { 32 });

            // Validate that complexity-metric is not used with fastga
            if is_fastga && matches!(complexity_metric, ComplexityMetric::DiagonalDistance) {
                error!("--complexity-metric cannot be used with --type fastga");
                error!("FastGA uses its own segmentation algorithm based on trace spacing");
                std::process::exit(1);
            }

            info!(
                "Converting CIGAR to {} tracepoints ({}={}, complexity-metric={})",
                tp_type,
                if is_fastga {
                    "trace_spacing"
                } else {
                    "max_complexity"
                },
                max_complexity,
                complexity_metric_to_str(complexity_metric)
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
                            process_compress_chunk(
                                &lines,
                                &tp_type,
                                max_complexity,
                                &complexity_metric,
                                minimal
                            );
                            lines.clear();
                        }
                    }
                    Err(e) => return Err(e.into()),
                }
            }

            // Process remaining lines
            if !lines.is_empty() {
                process_compress_chunk(&lines, &tp_type, max_complexity, &complexity_metric, minimal);
            }
        }
        Args::Decode {
            common,
            tp_type,
            sequence_files,
            sequence_list,
            penalties,
            trace_spacing,
            distance,
            complexity_metric,
            max_complexity,
            heuristic,
            keep_old_stats,
        } => {
            setup_logger(common.verbose);

            // Determine if we're using fastga based on tp_type
            let is_fastga = matches!(tp_type, TracepointType::Fastga);

            // Validate that complexity-metric is not used with fastga
            if is_fastga && matches!(complexity_metric, ComplexityMetric::DiagonalDistance) {
                error!("--complexity-metric cannot be used with --type fastga");
                error!("FastGA uses its own segmentation algorithm based on trace spacing");
                std::process::exit(1);
            }

            if heuristic && is_fastga {
                error!("--heuristic cannot be used with --type fastga");
                std::process::exit(1);
            }

            let heuristic_max_complexity = if heuristic {
                Some(max_complexity.unwrap_or_else(|| {
                    error!("--heuristic requires specifying --max-complexity");
                    std::process::exit(1);
                }))
            } else {
                if let Some(value) = max_complexity {
                    warn!(
                        "Ignoring --max-complexity={} because --heuristic was not requested",
                        value
                    );
                }
                None
            };

            if matches!(distance, DistanceChoice::Edit)
                && penalties.is_some()
                && !penalties.as_deref().unwrap_or("").is_empty()
            {
                error!("--penalties cannot be used with --distance edit (as it uses unit costs)");
                std::process::exit(1);
            }

            if is_fastga && penalties.is_some() && !penalties.as_deref().unwrap_or("").is_empty() {
                error!("--penalties should only be used without --type fastga");
                std::process::exit(1);
            }

            let sequence_paths = match collect_sequence_paths(sequence_files, sequence_list) {
                Ok(paths) if !paths.is_empty() => paths,
                Ok(_) => {
                    error!(
                        "At least one FASTA must be provided via --sequence-files or --sequence-list"
                    );
                    std::process::exit(1);
                }
                Err(msg) => {
                    error!("{}", msg);
                    std::process::exit(1);
                }
            };

            let sequence_index = match SequenceIndex::build(&sequence_paths) {
                Ok(index) => index,
                Err(msg) => {
                    error!("{}", msg);
                    std::process::exit(1);
                }
            };

            let penalties_value = match distance {
                DistanceChoice::Edit => None,
                _ => {
                    let default = distance.default_penalties().unwrap();
                    Some(penalties.clone().unwrap_or_else(|| default.to_string()))
                }
            };

            let sequence_index = Arc::new(sequence_index);

            let distance_mode = match parse_distance(distance, penalties_value.as_deref()) {
                Ok(dist) => dist,
                Err(msg) => {
                    error!("Invalid --distance configuration: {}", msg);
                    std::process::exit(1);
                }
            };

            if is_fastga && !matches!(distance, DistanceChoice::GapAffine2p) {
                warn!(
                    "--distance {} is ignored with --type fastga (uses edit distance)",
                    distance
                );
            }

            let penalties_summary = penalties_value.as_deref().unwrap_or("n/a");

            info!(
                "Converting {} tracepoints to CIGAR (complexity-metric={}, distance={}, penalties={}, heuristic={}{})",
                tp_type,
                complexity_metric_to_str(complexity_metric),
                distance,
                penalties_summary,
                if heuristic { "enabled" } else { "disabled" },
                heuristic_max_complexity.map(|mc| format!(", max-complexity={}", mc)).unwrap_or_default()
            );

            // Set the thread pool size
            rayon::ThreadPoolBuilder::new()
                .num_threads(common.threads)
                .build_global()?;

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
                                sequence_index.as_ref(),
                                &distance_mode,
                                is_fastga,
                                trace_spacing,
                                &complexity_metric,
                                heuristic,
                                heuristic_max_complexity,
                                keep_old_stats,
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
                    sequence_index.as_ref(),
                    &distance_mode,
                    is_fastga,
                    trace_spacing,
                    &complexity_metric,
                    heuristic,
                    heuristic_max_complexity,
                    keep_old_stats,
                );
            }
        }
        #[cfg(debug_assertions)]
        Args::Debug {
            paf,
            threads,
            sequence_files,
            sequence_list,
            penalties,
            max_complexity,
            verbose,
        } => {
            setup_logger(verbose);
            info!("Debugging");

            let (mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = parse_penalties(&penalties)?;
            info!(
                "Penalties: {},{},{},{},{}",
                mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2
            );

            let sequence_paths = match collect_sequence_paths(sequence_files, sequence_list) {
                Ok(paths) if !paths.is_empty() => paths,
                Ok(_) => {
                    error!(
                        "Debug mode requires at least one FASTA via --sequence-files or --sequence-list"
                    );
                    std::process::exit(1);
                }
                Err(msg) => {
                    error!("{}", msg);
                    std::process::exit(1);
                }
            };

            let sequence_index = match SequenceIndex::build(&sequence_paths) {
                Ok(index) => Arc::new(index),
                Err(msg) => {
                    error!("{}", msg);
                    std::process::exit(1);
                }
            };

            if let Some(paf) = paf {
                info!("PAF file: {}", paf);
                info!("Sequence FASTA files: {}", sequence_paths.join(", "));

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
                                    sequence_index.as_ref(),
                                    mismatch,
                                    gap_open1,
                                    gap_ext1,
                                    gap_open2,
                                    gap_ext2,
                                    max_complexity,
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
                        sequence_index.as_ref(),
                        mismatch,
                        gap_open1,
                        gap_ext1,
                        gap_open2,
                        gap_ext2,
                        max_complexity,
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
                let distance = Distance::GapAffine2p {
                    mismatch,
                    gap_opening1: gap_open1,
                    gap_extension1: gap_ext1,
                    gap_opening2: gap_open2,
                    gap_extension2: gap_ext2,
                };
                let aligner = distance.create_aligner(None);
                let paf_cigar = align_sequences_wfa(
                    &query_seq[a_start..a_end],
                    &target_seq[b_start..b_end],
                    &aligner,
                );
                let paf_cigar = cigar_ops_to_cigar_string(&paf_cigar);
                let tracepoints = cigar_to_tracepoints(
                    &paf_cigar,
                    max_complexity,
                    ComplexityMetric::EditDistance,
                );
                let cigar_from_tracepoints = tracepoints_to_cigar(
                    &tracepoints,
                    &query_seq,
                    &target_seq,
                    0,
                    0,
                    ComplexityMetric::EditDistance,
                    &distance,
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
    sequence_index: &SequenceIndex,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
    max_complexity: usize,
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

        let mut query_seq = sequence_index
            .fetch_sequence(query_name, query_start, query_end)
            .unwrap_or_else(|msg| {
                error!("{}", message_with_truncate_paf_file(&msg, line));
                std::process::exit(1);
            });
        if strand == "-" {
            query_seq = reverse_complement(&query_seq);
        }

        let target_seq = sequence_index
            .fetch_sequence(target_name, target_start, target_end)
            .unwrap_or_else(|msg| {
                error!("{}", message_with_truncate_paf_file(&msg, line));
                std::process::exit(1);
            });

        // Create thread-local aligner
        let distance = Distance::GapAffine2p {
            mismatch,
            gap_opening1: gap_open1,
            gap_extension1: gap_ext1,
            gap_opening2: gap_open2,
            gap_extension2: gap_ext2,
        };
        let aligner = distance.create_aligner(None);
        let realn_cigar = align_sequences_wfa(&query_seq, &target_seq, &aligner);
        let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);
        // let paf_cigar = &realn_cigar;

        // Convert CIGAR to tracepoints using query (A) and target (B) coordinates.
        let tracepoints =
            cigar_to_tracepoints(paf_cigar, max_complexity, ComplexityMetric::EditDistance);
        let variable_tracepoints = cigar_to_variable_tracepoints(
            paf_cigar,
            max_complexity,
            ComplexityMetric::EditDistance,
        );

        // Also convert using raw functions
        let tracepoints_raw =
            cigar_to_tracepoints_raw(paf_cigar, max_complexity, ComplexityMetric::EditDistance);
        let variable_tracepoints_raw = cigar_to_variable_tracepoints_raw(
            paf_cigar,
            max_complexity,
            ComplexityMetric::EditDistance,
        );

        // Convert using diagonal distance functions
        let tracepoints_diagonal = cigar_to_tracepoints(
            paf_cigar,
            max_complexity,
            ComplexityMetric::DiagonalDistance,
        );
        let mixed_tracepoints_diagonal = cigar_to_mixed_tracepoints(
            paf_cigar,
            max_complexity,
            ComplexityMetric::DiagonalDistance,
        );
        let variable_tracepoints_diagonal = cigar_to_variable_tracepoints(
            paf_cigar,
            max_complexity,
            ComplexityMetric::DiagonalDistance,
        );

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
        let distance = Distance::GapAffine2p {
            mismatch,
            gap_opening1: gap_open1,
            gap_extension1: gap_ext1,
            gap_opening2: gap_open2,
            gap_extension2: gap_ext2,
        };
        let cigar_from_tracepoints = tracepoints_to_cigar(
            &tracepoints,
            &query_seq,
            &target_seq,
            0,
            0,
            ComplexityMetric::EditDistance,
            &distance,
        );
        let cigar_from_variable_tracepoints = variable_tracepoints_to_cigar(
            &variable_tracepoints,
            &query_seq,
            &target_seq,
            0,
            0,
            ComplexityMetric::EditDistance,
            &distance,
        );

        // Reconstruct CIGAR from raw tracepoints
        let cigar_from_tracepoints_raw = tracepoints_to_cigar(
            &tracepoints_raw,
            &query_seq,
            &target_seq,
            0,
            0,
            ComplexityMetric::EditDistance,
            &distance,
        );
        let cigar_from_variable_tracepoints_raw = variable_tracepoints_to_cigar(
            &variable_tracepoints_raw,
            &query_seq,
            &target_seq,
            0,
            0,
            ComplexityMetric::EditDistance,
            &distance,
        );

        // Reconstruct CIGAR from diagonal tracepoints
        let cigar_from_tracepoints_diagonal = tracepoints_to_cigar(
            &tracepoints_diagonal,
            &query_seq,
            &target_seq,
            0,
            0,
            ComplexityMetric::DiagonalDistance,
            &distance,
        );
        let cigar_from_mixed_tracepoints_diagonal = mixed_tracepoints_to_cigar(
            &mixed_tracepoints_diagonal,
            &query_seq,
            &target_seq,
            0,
            0,
            ComplexityMetric::DiagonalDistance,
            &distance,
        );
        let cigar_from_variable_tracepoints_diagonal = variable_tracepoints_to_cigar(
            &variable_tracepoints_diagonal,
            &query_seq,
            &target_seq,
            0,
            0,
            ComplexityMetric::DiagonalDistance,
            &distance,
        );

        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, paf_gap_compressed_id, paf_block_id) = calculate_cigar_stats(paf_cigar);
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
        let score_from_paf = compute_alignment_score_from_cigar(paf_cigar, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
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
            println!("\t                           bounds CIGAR from PAF: {:?}", get_cigar_diagonal_bounds(paf_cigar));
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
            println!();
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
        if c.is_ascii_digit() {
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
        if c.is_ascii_digit() {
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
fn process_compress_chunk(
    lines: &[String],
    tp_type: &TracepointType,
    max_complexity: usize,
    complexity_metric: &ComplexityMetric,
    minimal: bool
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
        let cigar = &cg_field[5..];

        // Calculate identity stats from the CIGAR
        let (gap_compressed_identity, block_identity) = calculate_identity_stats(cigar);

        // Calculate alignment score based on edit distance
        let alignment_score = calculate_alignment_score_edit_distance(cigar);

        // Check for existing df:i: field
        let existing_df = fields
            .iter()
            .find(|&&s| s.starts_with("df:i:"))
            .and_then(|&s| s[5..].parse::<usize>().ok());

        // Handle FastGA tracepoints with potential overflow splitting
        if matches!(tp_type, TracepointType::Fastga) {
            process_fastga_with_overflow(
                &fields,
                cigar,
                max_complexity,
                gap_compressed_identity,
                block_identity,
                alignment_score,
                existing_df,
                complexity_metric,
                minimal
            );
        } else {
            // Handle other tracepoint types (no overflow splitting needed)
            process_single_record(
                &fields,
                cigar,
                tp_type,
                max_complexity,
                gap_compressed_identity,
                block_identity,
                alignment_score,
                existing_df,
                complexity_metric,
                minimal
            );
        }
    });
}

fn process_fastga_with_overflow(
    fields: &[&str],
    cigar: &str,
    max_complexity: usize,
    gap_compressed_identity: f64,
    block_identity: f64,
    alignment_score: i32,
    existing_df: Option<usize>,
    _complexity_metric: &ComplexityMetric,
    minimal: bool,
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
        max_complexity,
        query_start,
        query_end,
        query_len,
        target_start,
        target_end,
        target_len,
        complement,
    );

    // Debug: show all segments
    debug!("Total segments from lib_tracepoints: {}", segments.len());
    for (i, (tps, (qs, qe, ts, te))) in segments.iter().enumerate() {
        debug!("  Segment {}: q={}..{}, t={}..{}, {} tracepoints", i, qs, qe, ts, te, tps.len());
    }

    // For FASTGA mode, output MULTIPLE records (one per segment) like PAFtoALN.c does
    // Each segment has its own coordinates and tracepoints
    // Gaps between segments are implicit in coordinate discontinuities

    let mut valid_segment_count = 0;
    for (segment_idx, (tracepoints, (seg_query_start, seg_query_end, seg_target_start, seg_target_end)))
        in segments.iter().enumerate()
    {
        let seg_query_len = seg_query_end.saturating_sub(*seg_query_start);
        let seg_target_len = seg_target_end.saturating_sub(*seg_target_start);

        // Skip zero-length segments (gap markers with no alignment content)
        if seg_query_len == 0 && seg_target_len == 0 {
            info!("Skipping zero-length segment {}", segment_idx);
            continue;
        }

        // Warn if segment has no tracepoints (shouldn't happen for valid segments)
        if tracepoints.is_empty() {
            warn!(
                "Segment {} has coordinates (q:{}-{}, t:{}-{}) but no tracepoints",
                segment_idx, seg_query_start, seg_query_end, seg_target_start, seg_target_end
            );
            continue;
        }

        valid_segment_count += 1;

        // Convert target coordinates back to original space for output
        let (out_target_start, out_target_end) = if complement {
            // Reversed space -> Original space
            (target_len - seg_target_end, target_len - seg_target_start)
        } else {
            (*seg_target_start, *seg_target_end)
        };

        // Calculate segment-specific stats
        let sum_of_differences: usize = tracepoints.iter().map(|(diff, _)| diff).sum();
        let tracepoints_str = format_tracepoints(tracepoints);
        let alignment_length = seg_query_end - seg_query_start;
        let matches = alignment_length.saturating_sub(sum_of_differences);

        // Create output fields for this segment
        let mut new_fields: Vec<String> = Vec::new();

        for (i, field) in fields.iter().enumerate() {
            match i {
                2 => new_fields.push(seg_query_start.to_string()),
                3 => new_fields.push(seg_query_end.to_string()),
                7 => new_fields.push(out_target_start.to_string()),
                8 => new_fields.push(out_target_end.to_string()),
                9 => new_fields.push(matches.to_string()),
                10 => new_fields.push(alignment_length.to_string()),
                _ => {
                    if field.starts_with("cg:Z:") {
                        // Add identity stats before tracepoints
                        if !minimal {
                            new_fields.push(format!("gi:f:{:.12}", gap_compressed_identity));
                            new_fields.push(format!("bi:f:{:.12}", block_identity));

                            // Add df fields
                            if let Some(old_df) = existing_df {
                                new_fields.push(format!("dfold:i:{}", old_df));
                            }
                            new_fields.push(format!("df:i:{}", sum_of_differences));

                            // Add alignment score
                            new_fields.push(format!("sc:i:{}", alignment_score));
                        }

                        // Add tracepoints
                        new_fields.push(format!("tp:Z:{}", tracepoints_str));
                    } else if field.starts_with("df:i:")
                        || field.starts_with("gi:f:")
                        || field.starts_with("bi:f:")
                        || field.starts_with("sc:i:")
                    {
                        // Skip existing fields as we handle them above
                        continue;
                    } else {
                        new_fields.push(field.to_string());
                    }
                }
            }
        }

        // Add segment marker if there are multiple valid segments
        if segments.len() > 1 {
            new_fields.push(format!("sg:i:{}", valid_segment_count - 1));
        }

        println!("{}", new_fields.join("\t"));
    }
}

fn process_single_record(
    fields: &[&str],
    cigar: &str,
    tp_type: &TracepointType,
    max_complexity: usize,
    gap_compressed_identity: f64,
    block_identity: f64,
    alignment_score: i32,
    existing_df: Option<usize>,
    complexity_metric: &ComplexityMetric,
    minimal: bool,
) {
    // Convert CIGAR based on tracepoint type and complexity metric
    let (tracepoints_str, df_value) = match tp_type {
        TracepointType::Standard => {
            let tp = cigar_to_tracepoints(cigar, max_complexity, *complexity_metric);
            (format_tracepoints(&tp), None)
        }
        TracepointType::Mixed => {
            let tp = cigar_to_mixed_tracepoints(cigar, max_complexity, *complexity_metric);
            (format_mixed_tracepoints(&tp), None::<usize>)
        }
        TracepointType::Variable => {
            let tp = cigar_to_variable_tracepoints(cigar, max_complexity, *complexity_metric);
            (format_variable_tracepoints(&tp), None)
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
            if !minimal {
                new_fields.push(format!("gi:f:{:.12}", gap_compressed_identity));
                new_fields.push(format!("bi:f:{:.12}", block_identity));
                new_fields.push(format!("sc:i:{}", alignment_score));
                if let Some(new_df) = df_value {
                    if let Some(old_df) = existing_df {
                        new_fields.push(format!("dfold:i:{}", old_df));
                    }
                    new_fields.push(format!("df:i:{}", new_df));
                }
            }

            // Add tracepoints
            new_fields.push(format!("tp:Z:{}", tracepoints_str));
        } else if field.starts_with("df:i:")
            || field.starts_with("gi:f:")
            || field.starts_with("bi:f:")
            || field.starts_with("sc:i:")
        {
            // Skip existing fields as we handle them above
            continue;
        } else {
            new_fields.push(field.to_string());
        }
    }

    println!("{}", new_fields.join("\t"));
}

/// Calculate residue_matches and alignment_block_length from CIGAR string
fn calculate_paf_fields_from_cigar(cigar: &str) -> (usize, usize) {
    let mut matches = 0;
    let mut mismatches = 0;
    let mut insertions = 0;
    let mut deletions = 0;
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
                'I' => insertions += len,
                'D' => deletions += len,
                'S' | 'H' | 'P' | 'N' => {} // Skip clipping operations
                _ => {}
            }
        }
    }

    let residue_matches = matches;
    let alignment_block_length = matches + mismatches + insertions + deletions;
    
    (residue_matches, alignment_block_length)
}

/// Process a chunk of lines in parallel for decompression
fn process_decompress_chunk(
    lines: &[String],
    tp_type: &TracepointType,
    sequence_index: &SequenceIndex,
    distance_mode: &Distance,
    fastga: bool,
    trace_spacing: usize,
    complexity_metric: &ComplexityMetric,
    heuristic: bool,
    heuristic_max_complexity: Option<usize>,
    keep_old_stats: bool,
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

        debug!(
            "Fetching query sequence {}:{}-{} on {} strand",
            query_name,
            query_start,
            query_end,
            if strand == "-" && !fastga { "-" } else { "+" }
        );
        let mut query_seq = sequence_index
            .fetch_sequence(query_name, query_start, query_end)
            .unwrap_or_else(|msg| {
                error!("{}", message_with_truncate_paf_file(&msg, line));
                std::process::exit(1);
            });
        if strand == "-" && !fastga {
            query_seq = reverse_complement(&query_seq);
        }

        debug!(
            "Fetching target sequence {}:{}-{} on {} strand",
            target_name,
            target_start,
            target_end,
            if strand == "-" && fastga { "-" } else { "+" }
        );
        let mut target_seq = sequence_index
            .fetch_sequence(target_name, target_start, target_end)
            .unwrap_or_else(|msg| {
                error!("{}", message_with_truncate_paf_file(&msg, line));
                std::process::exit(1);
            });
        if fastga && strand == "-" {
            target_seq = reverse_complement(&target_seq);
        }

        // Use specified tracepoint type
        let distance = distance_mode.clone();

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
                if strand == "+" {
                    target_len - target_end
                } else {
                    target_start
                },
                strand == "-",
            )
        } else {
            match tp_type {
                TracepointType::Mixed => {
                    let mixed_tracepoints = parse_mixed_tracepoints(tracepoints_str);
                    match *complexity_metric {
                        ComplexityMetric::EditDistance => {
                            if heuristic {
                                let max_value = heuristic_max_complexity
                                    .expect("missing max-complexity with heuristic");
                                let mut aligner = distance.create_aligner(None);
                                mixed_tracepoints_to_cigar_with_aligner(
                                    &mixed_tracepoints,
                                    &query_seq,
                                    &target_seq,
                                    0,
                                    0,
                                    ComplexityMetric::EditDistance,
                                    &mut aligner,
                                    true,
                                    max_value,
                                )
                            } else {
                                mixed_tracepoints_to_cigar(
                                    &mixed_tracepoints,
                                    &query_seq,
                                    &target_seq,
                                    0,
                                    0,
                                    ComplexityMetric::EditDistance,
                                    &distance,
                                )
                            }
                        }
                        ComplexityMetric::DiagonalDistance => mixed_tracepoints_to_cigar(
                            &mixed_tracepoints,
                            &query_seq,
                            &target_seq,
                            0,
                            0,
                            ComplexityMetric::DiagonalDistance,
                            &distance,
                        ),
                    }
                }
                TracepointType::Variable => {
                    let variable_tracepoints = parse_variable_tracepoints(tracepoints_str);
                    match *complexity_metric {
                        ComplexityMetric::EditDistance => {
                            if heuristic {
                                let max_value = heuristic_max_complexity
                                    .expect("missing max-complexity with heuristic");
                                let mut aligner = distance.create_aligner(None);
                                variable_tracepoints_to_cigar_with_aligner(
                                    &variable_tracepoints,
                                    &query_seq,
                                    &target_seq,
                                    0,
                                    0,
                                    ComplexityMetric::EditDistance,
                                    &mut aligner,
                                    true,
                                    max_value,
                                )
                            } else {
                                variable_tracepoints_to_cigar(
                                    &variable_tracepoints,
                                    &query_seq,
                                    &target_seq,
                                    0,
                                    0,
                                    ComplexityMetric::EditDistance,
                                    &distance,
                                )
                            }
                        }
                        ComplexityMetric::DiagonalDistance => variable_tracepoints_to_cigar(
                            &variable_tracepoints,
                            &query_seq,
                            &target_seq,
                            0,
                            0,
                            ComplexityMetric::DiagonalDistance,
                            &distance,
                        ),
                    }
                }
                TracepointType::Standard => {
                    let tracepoints = parse_tracepoints(tracepoints_str);
                    match *complexity_metric {
                        ComplexityMetric::EditDistance => {
                            if heuristic {
                                let max_value = heuristic_max_complexity
                                    .expect("missing max-complexity with heuristic");
                                let mut aligner = distance.create_aligner(None);
                                tracepoints_to_cigar_with_aligner(
                                    &tracepoints,
                                    &query_seq,
                                    &target_seq,
                                    0,
                                    0,
                                    ComplexityMetric::EditDistance,
                                    &mut aligner,
                                    true,
                                    max_value,
                                )
                            } else {
                                tracepoints_to_cigar(
                                    &tracepoints,
                                    &query_seq,
                                    &target_seq,
                                    0,
                                    0,
                                    ComplexityMetric::EditDistance,
                                    &distance,
                                )
                            }
                        }
                        ComplexityMetric::DiagonalDistance => tracepoints_to_cigar(
                            &tracepoints,
                            &query_seq,
                            &target_seq,
                            0,
                            0,
                            ComplexityMetric::DiagonalDistance,
                            &distance,
                        ),
                    }
                }
                TracepointType::Fastga => {
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

        let (residue_matches, alignment_block_length) = calculate_paf_fields_from_cigar(&cigar);

        // Build the new line
        let mut new_fields: Vec<String> = Vec::new();

        for (i, field) in fields.iter().enumerate() {
            match i {
                9 => {
                    // Update residue_matches from reconstructed CIGAR
                    new_fields.push(residue_matches.to_string());
                }
                10 => {
                    // Update alignment_block_length from reconstructed CIGAR
                    new_fields.push(alignment_block_length.to_string());
                }
                _ => {
                    if field.starts_with("tp:Z:") {
                        // Optionally keep old values
                        if keep_old_stats {
                            if let Some(old_gi) = existing_gi {
                                new_fields.push(format!("giold:f:{}", &old_gi[5..]));
                            }
                            if let Some(old_bi) = existing_bi {
                                new_fields.push(format!("biold:f:{}", &old_bi[5..]));
                            }
                            if let Some(old_sc) = existing_sc {
                                new_fields.push(format!("scold:i:{}", &old_sc[5..]));
                            }
                        }

                        // Always add new identity stats and alignment score before the CIGAR
                        new_fields.push(format!("gi:f:{:.12}", gap_compressed_identity));
                        new_fields.push(format!("bi:f:{:.12}", block_identity));
                        new_fields.push(format!("sc:i:{}", alignment_score));

                        // Replace tracepoints with CIGAR
                        new_fields.push(format!("cg:Z:{}", cigar));
                    } else if field.starts_with("gi:f:")
                        || field.starts_with("bi:f:")
                        || field.starts_with("sc:i:")
                    {
                        // Skip existing gi, bi, and sc fields - they will be replaced
                        continue;
                    } else {
                        new_fields.push(field.to_string());
                    }
                }
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
                let chars = s.chars();
                let mut len_str = String::new();

                // Read digits
                for c in chars {
                    if c.is_ascii_digit() {
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

fn parse_distance(distance: DistanceChoice, penalties: Option<&str>) -> Result<Distance, String> {
    match distance {
        DistanceChoice::Edit => Ok(Distance::Edit {}),
        DistanceChoice::GapAffine => {
            let values = parse_penalty_values(penalties, 3, "mismatch,gap_open,gap_ext")?;
            Ok(Distance::GapAffine {
                mismatch: values[0],
                gap_opening: values[1],
                gap_extension: values[2],
            })
        }
        DistanceChoice::GapAffine2p => {
            let values = parse_penalty_values(
                penalties,
                5,
                "mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2",
            )?;
            Ok(Distance::GapAffine2p {
                mismatch: values[0],
                gap_opening1: values[1],
                gap_extension1: values[2],
                gap_opening2: values[3],
                gap_extension2: values[4],
            })
        }
    }
}

fn parse_penalty_values(
    penalties: Option<&str>,
    expected: usize,
    description: &str,
) -> Result<Vec<i32>, String> {
    let Some(penalties) = penalties else {
        return Err(format!("missing --penalties (expected {})", description));
    };

    let values: Vec<i32> = penalties
        .split(',')
        .map(|token| token.trim())
        .filter(|token| !token.is_empty())
        .map(|token| {
            token
                .parse::<i32>()
                .map_err(|e| format!("invalid penalty '{}': {}", token, e))
        })
        .collect::<Result<Vec<_>, _>>()?;

    if values.len() != expected {
        return Err(format!(
            "expected {} values ({}), found {}",
            expected,
            description,
            values.len()
        ));
    }

    Ok(values)
}

#[cfg(debug_assertions)]
fn parse_penalties(
    penalties: &str,
) -> Result<(i32, i32, i32, i32, i32), Box<dyn std::error::Error>> {
    let values = parse_penalty_values(
        Some(penalties),
        5,
        "mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2",
    )
    .map_err(|msg| -> Box<dyn std::error::Error> { msg.into() })?;
    Ok((values[0], values[1], values[2], values[3], values[4]))
}

#[cfg(debug_assertions)]
fn get_cigar_diagonal_bounds(cigar: &str) -> (i64, i64) {
    let mut current_diagonal = 0; // Current diagonal position
    let mut min_diagonal = 0; // Lowest diagonal reached
    let mut max_diagonal = 0; // Highest diagonal reached

    // Parse CIGAR string with numerical counts
    let mut num_buffer = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
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
        if c.is_ascii_digit() {
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
