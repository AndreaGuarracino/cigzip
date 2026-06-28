mod sequence;

use crate::sequence::{collect_sequence_paths, SequenceIndex};
use clap::{Parser, ValueEnum};
use flate2::read::MultiGzDecoder;
#[cfg(debug_assertions)]
use indicatif::ProgressBar;
#[cfg(debug_assertions)]
use indicatif::ProgressStyle;
use lib_wfa2::affine_wavefront::Distance;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};
// tpa: CIGAR stats, formatting, and tracepoints → CIGAR reconstruction
#[cfg(debug_assertions)]
use tpa::calculate_alignment_score;
use tpa::{format_ab_tracepoints, format_tracepoints, reconstruct_cigar_with_aligner, CigarStats};
// tracepoints: CIGAR → tracepoints encoding
#[cfg(debug_assertions)]
use tracepoints::{
    align_sequences_wfa, cigar_ops_to_cigar_string, cigar_to_tracepoints_raw,
    cigar_to_variable_tracepoints_raw, mixed_tracepoints_to_cigar, tracepoints_to_cigar,
    variable_tracepoints_to_cigar,
};
use std::collections::HashMap;
use tracepoints::{
    cigar_to_mixed_tracepoints, cigar_to_tracepoints_fastga,
    cigar_to_tracepoints_fastga_with_contigs, cigar_to_tracepoints_from_ops,
    cigar_to_variable_tracepoints, ComplexityMetric, TracepointData, TracepointType,
};

/// Load a FASTGA contig table: lines "scaffold<TAB>sbeg<TAB>send" -> scaffold -> sorted (sbeg,send).
/// Zero-length contigs (leading-gap artifacts) are dropped.
fn load_contig_table(path: &str) -> io::Result<HashMap<String, Vec<(usize, usize)>>> {
    let mut map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let (sbeg, send) = (
            f[1].parse::<usize>().unwrap_or(0),
            f[2].parse::<usize>().unwrap_or(0),
        );
        if send > sbeg {
            map.entry(f[0].to_string()).or_default().push((sbeg, send));
        }
    }
    for v in map.values_mut() {
        v.sort_unstable();
    }
    Ok(map)
}

/// Start of the contig (ACGT run) containing scaffold position `pos`.
/// Returns 0 when `scaffold` has no table entry (contig-aware decoding off) or `pos` is in no
/// contig; both keep query_offset = query_start (absolute grid).
fn contig_start_in_table(
    table: &HashMap<String, Vec<(usize, usize)>>,
    scaffold: &str,
    pos: usize,
) -> usize {
    match table.get(scaffold) {
        Some(contigs) => contigs
            .iter()
            .find(|&&(s, e)| s <= pos && pos < e)
            .map(|&(s, _)| s)
            .unwrap_or(0),
        None => 0,
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

/// Memory mode for WFA aligner
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, Default)]
enum MemoryModeChoice {
    #[default]
    High,
    Medium,
    Low,
    Ultralow,
    Adaptive,
}

impl MemoryModeChoice {
    fn to_lib_wfa2(self) -> lib_wfa2::affine_wavefront::MemoryMode {
        match self {
            Self::High | Self::Adaptive => lib_wfa2::affine_wavefront::MemoryMode::High,
            Self::Medium => lib_wfa2::affine_wavefront::MemoryMode::Medium,
            Self::Low => lib_wfa2::affine_wavefront::MemoryMode::Low,
            Self::Ultralow => lib_wfa2::affine_wavefront::MemoryMode::Ultralow,
        }
    }
}

impl fmt::Display for MemoryModeChoice {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MemoryModeChoice::High => write!(f, "high"),
            MemoryModeChoice::Medium => write!(f, "medium"),
            MemoryModeChoice::Low => write!(f, "low"),
            MemoryModeChoice::Ultralow => write!(f, "ultralow"),
            MemoryModeChoice::Adaptive => write!(f, "adaptive"),
        }
    }
}

/// Largest single-segment span in a tracepoint record. Used by adaptive memory-mode selection.
fn max_segment_span(tp: &TracepointData) -> usize {
    use tracepoints::MixedRepresentation;
    match tp {
        TracepointData::Standard(v) | TracepointData::Fastga(v) => {
            v.iter().map(|(a, b)| (*a).max(*b)).max().unwrap_or(0)
        }
        TracepointData::FastgaNoDiff(v) => v.iter().copied().max().unwrap_or(0),
        TracepointData::Variable(v) => {
            v.iter().map(|(a, b)| (*a).max(b.unwrap_or(0))).max().unwrap_or(0)
        }
        TracepointData::Mixed(v) => v
            .iter()
            .map(|it| match it {
                MixedRepresentation::Tracepoint(a, b) => (*a).max(*b),
                MixedRepresentation::CigarOp(len, _) => *len,
            })
            .max()
            .unwrap_or(0),
    }
}

/// Per-thread aligner(s) for decode.
struct Aligners {
    primary: lib_wfa2::affine_wavefront::AffineWavefronts,
    /// Present only in adaptive mode.
    medium: Option<lib_wfa2::affine_wavefront::AffineWavefronts>,
    /// Adaptive applies only where medium actually helps: DB-TP + gap-affine2p + mc>=64.
    gate: bool,
    /// Route a record to the medium aligner when `max_segment_span * divergence > threshold`.
    threshold: f64,
}

impl Aligners {
    /// Aligner for one record, given its largest segment span and input divergence.
    #[inline]
    fn select(
        &mut self,
        max_seg: usize,
        divergence: f64,
    ) -> &mut lib_wfa2::affine_wavefront::AffineWavefronts {
        match self.medium.as_mut() {
            Some(med) if self.gate && (max_seg as f64) * divergence > self.threshold => med,
            _ => &mut self.primary,
        }
    }
}

/// Build a per-thread `Aligners`: the primary aligner in 
// `primary_mode`, plus a medium aligner only when `adaptive`.
fn make_aligners(
    distance_mode: &Distance,
    primary_mode: &lib_wfa2::affine_wavefront::MemoryMode,
    adaptive: bool,
    gate: bool,
    threshold: u64,
) -> Aligners {
    let primary = distance_mode.create_aligner(None, Some(primary_mode));
    let medium = if adaptive {
        Some(distance_mode.create_aligner(None, Some(&lib_wfa2::affine_wavefront::MemoryMode::Medium)))
    } else {
        None
    };
    Aligners {
        primary,
        medium,
        gate,
        threshold: threshold as f64,
    }
}

/// Common options shared between all commands
#[derive(Parser, Debug)]
struct CommonOpts {
    /// PAF file for alignments (use "-" to read from standard input)
    #[arg(short = 'p', long = "paf")]
    paf: String,

    /// Number of threads to use (default: 4)
    #[arg(short, long = "threads", default_value_t = 4)]
    threads: usize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[arg(short, long, default_value = "1")]
    verbose: u8,
}

/// Tracepoint conversion options shared by encode/decode/compress.
#[derive(Parser, Debug)]
#[command(next_help_heading = "Tracepoint options")]
struct TracepointOpts {
    /// Tracepoint type (standard, mixed, variable, fastga, fastga-no-diff)
    #[arg(long = "type", default_value = "standard", value_parser = TracepointType::from_str, value_name = "TYPE")]
    tp_type: TracepointType,

    /// Maximum complexity (max_diff for standard/mixed/variable; trace spacing for fastga; default 32, 100 for fastga)
    #[arg(long = "max-complexity")]
    max_complexity: Option<u32>,

    /// Complexity metric (edit-distance, diagonal-distance). Not used with fastga
    #[arg(long = "complexity-metric", value_parser = ComplexityMetric::from_str, value_name = "METRIC")]
    complexity_metric: Option<ComplexityMetric>,

    /// Distance metric (edit, gap-affine, gap-affine2p)
    #[arg(long = "distance", default_value_t = DistanceChoice::Edit)]
    distance: DistanceChoice,

    /// Gap penalties for gap-affine distances (mismatch,go1,ge1[,go2,ge2]; ignored with edit)
    #[arg(long)]
    penalties: Option<String>,

    /// FASTGA contig table (TSV: scaffold<TAB>sbeg<TAB>send per ACGT run, from the GDB). Enables
    /// contig-aware fastga: encode/compress split at N-gaps (byte-identical to PAFtoALN); decode
    /// places the trace grid contig-locally.
    #[arg(long = "fastga-contigs")]
    fastga_contigs: Option<String>,
}

#[derive(Parser, Debug)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Encode CIGAR to tracepoints (text PAF cg:Z: → text PAF tp:Z:)
    Encode {
        #[clap(flatten)]
        common: CommonOpts,

        #[clap(flatten)]
        tp: TracepointOpts,

        /// Output file (use "-" for stdout)
        #[arg(short = 'o', long = "output")]
        output: Option<String>,

        /// Skip adding optional fields (gi/bi/sc fields)
        #[arg(long = "minimal")]
        minimal: bool,
    },
    /// Decode tracepoints to CIGAR (text PAF tp:Z: → text PAF cg:Z:)
    Decode {
        #[clap(flatten)]
        common: CommonOpts,

        #[clap(flatten)]
        tp: TracepointOpts,

        /// Sequence files in FASTA or AGC format (repeatable)
        #[arg(long = "sequence-files", value_name = "FILE", num_args = 1..)]
        sequence_files: Vec<String>,

        /// File listing sequence file paths (one per line)
        #[arg(long = "sequence-list", value_name = "FILE")]
        sequence_list: Option<String>,

        /// Keep original gi/bi/sc/df fields as go/bo/so/do
        #[arg(long = "keep-old-stats")]
        keep_old_stats: bool,

        /// Trace spacing for fastga (default: 100)
        #[arg(long)]
        trace_spacing: Option<u32>,

        /// Disable banded alignment (banded is on by default)
        #[arg(long = "no-banded")]
        no_banded: bool,

        /// Memory mode for WFA aligner (high, medium, low, ultralow, adaptive)
        #[arg(long = "memory-mode", default_value_t = MemoryModeChoice::High)]
        memory_mode: MemoryModeChoice,

        /// Adaptive threshold: route a record to medium when max_segment_span * divergence
        /// exceeds this (only with --memory-mode adaptive, DB-TP + gap-affine2p + mc>=64)
        #[arg(long = "adaptive-threshold", default_value_t = 5000)]
        adaptive_threshold: u64,
    },
    /// Compress text PAF to binary TPA (text PAF → binary TPA)
    Compress {
        /// Input text PAF file (accepts cg:Z: or tp:Z: tags)
        #[arg(short = 'i', long = "input", visible_alias = "paf")]
        input: String,

        /// Output binary TPA file
        #[arg(short = 'o', long = "output")]
        output: String,

        #[clap(flatten)]
        tp: TracepointOpts,

        /// Compression strategy (automatic, benchmark, raw, zigzag-delta, 2d-delta, etc.)
        #[arg(
            long = "strategy",
            default_value = "automatic",
            value_name = "STRATEGY"
        )]
        strategy_str: String,

        /// Compress all records together with BGZIP
        #[arg(long = "all-records")]
        all_records: bool,

        /// Number of threads to use (default: 4)
        #[arg(short, long = "threads", default_value_t = 4)]
        threads: usize,

        /// Verbosity level (0 = error, 1 = info, 2 = debug)
        #[arg(short, long, default_value = "1")]
        verbose: u8,
    },
    /// Decompress binary TPA to text PAF (binary TPA → text PAF)
    Decompress {
        /// Input binary TPA file
        #[arg(short = 'i', long = "input")]
        input: String,

        /// Output text PAF file (use "-" for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Also decode tracepoints to CIGAR (requires --sequence-files)
        #[arg(long = "decode")]
        decode: bool,

        /// Sequence files in FASTA or AGC format (required with --decode)
        #[arg(long = "sequence-files", value_name = "FILE", num_args = 1.., required_if_eq("decode", "true"))]
        sequence_files: Vec<String>,

        /// File listing sequence file paths (one per line)
        #[arg(long = "sequence-list", value_name = "FILE")]
        sequence_list: Option<String>,

        /// FASTGA contig table (required with --decode to correctly reconstruct contig-aware FASTGA TPA)
        #[arg(long = "fastga-contigs", value_name = "FILE")]
        fastga_contigs: Option<String>,

        /// Keep original gi/bi/sc/df fields as go/bo/so/do
        #[arg(long = "keep-old-stats")]
        keep_old_stats: bool,

        /// Memory mode for WFA aligner (high, medium, low, ultralow, adaptive)
        #[arg(long = "memory-mode", default_value_t = MemoryModeChoice::High)]
        memory_mode: MemoryModeChoice,

        /// Adaptive threshold: route a record to medium when max_segment_span * divergence
        /// exceeds this (only with --memory-mode adaptive, DB-TP + gap-affine2p + mc>=64)
        #[arg(long = "adaptive-threshold", default_value_t = 5000)]
        adaptive_threshold: u64,

        /// Verbosity level (0 = error, 1 = info, 2 = debug)
        #[arg(short, long, default_value = "1")]
        verbose: u8,
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

        /// Sequence files in FASTA or AGC format (repeatable)
        #[arg(long = "sequence-files", value_name = "FILE", num_args = 1..)]
        sequence_files: Vec<String>,

        /// File listing sequence file paths (one per line)
        #[arg(long = "sequence-list", value_name = "FILE")]
        sequence_list: Option<String>,

        /// Gap penalties in the format mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2
        #[arg(long, default_value = "5,8,2,24,1")]
        penalties: String,

        /// Maximum complexity value for tracepoint segmentation
        #[arg(long = "max-complexity", default_value = "32")]
        max_complexity: u32,

        /// Verbosity level (0 = error, 1 = info, 2 = debug)
        #[arg(short, long, default_value = "0")]
        verbose: u8,
    },
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments.
    let args = Args::parse();

    match args {
        Args::Encode {
            common,
            tp,
            output,
            minimal,
        } => {
            let TracepointOpts {
                tp_type,
                max_complexity,
                complexity_metric,
                distance,
                penalties,
                fastga_contigs,
            } = tp;
            setup_logger(common.verbose);

            // Set the thread pool size before any rayon use
            rayon::ThreadPoolBuilder::new()
                .num_threads(common.threads)
                .build_global()?;

            let is_fastga = matches!(tp_type, TracepointType::Fastga | TracepointType::FastgaNoDiff);
            let max_complexity = max_complexity.unwrap_or(if is_fastga { 100 } else { 32 });

            if is_fastga && complexity_metric.is_some() {
                error!("--complexity-metric cannot be used with --type fastga");
                error!("FastGA uses its own segmentation algorithm based on trace spacing");
                std::process::exit(1);
            }

            // Default to EditDistance for non-fastga if not specified
            let complexity_metric = complexity_metric.unwrap_or(ComplexityMetric::EditDistance);

            // Parse distance model for score calculation
            let distance_mode = match parse_distance(distance, penalties.as_deref()) {
                Ok(dist) => dist,
                Err(msg) => {
                    error!("{}", msg);
                    std::process::exit(1);
                }
            };

            if is_fastga {
                info!(
                    "Encoding CIGAR to fastga tracepoints (trace_spacing={})",
                    max_complexity
                );
            } else {
                info!(
                    "Encoding CIGAR to {} tracepoints (max_complexity={}, complexity-metric={})",
                    tp_type.as_str(),
                    max_complexity,
                    complexity_metric
                );
            }

            // Create output writer
            let writer: Arc<Mutex<Box<dyn Write + Send>>> = if let Some(output_path) = output {
                let file = File::create(&output_path)?;
                Arc::new(Mutex::new(Box::new(BufWriter::new(file))))
            } else {
                Arc::new(Mutex::new(Box::new(BufWriter::new(std::io::stdout()))))
            };

            // Optional FASTGA contig table: scaffold -> sorted (sbeg,send) ACGT runs.
            // Enables PAFtoALN-identical contig-aware fastga encoding (split at query N-gaps).
            let contig_table: HashMap<String, Vec<(usize, usize)>> = match &fastga_contigs {
                Some(path) => match load_contig_table(path) {
                    Ok(t) => {
                        info!("Loaded fastga contig table: {} scaffolds from {}", t.len(), path);
                        t
                    }
                    Err(e) => {
                        error!("Failed to read --fastga-contigs {}: {}", path, e);
                        std::process::exit(1);
                    }
                },
                None => HashMap::new(),
            };

            let paf_reader = get_paf_reader(&common.paf)?;

            paf_reader
                .lines()
                .map_while(Result::ok)
                .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
                .par_bridge()
                .for_each(|line| {
                    process_compress_line(
                        &line,
                        &tp_type,
                        max_complexity,
                        &complexity_metric,
                        &distance_mode,
                        &writer,
                        minimal,
                        &contig_table,
                    );
                });

            // Flush the writer
            writer.lock().unwrap().flush()?;
        }
        Args::Decode {
            common,
            tp,
            sequence_files,
            sequence_list,
            trace_spacing,
            no_banded,
            keep_old_stats,
            memory_mode,
            adaptive_threshold,
        } => {
            let TracepointOpts {
                tp_type,
                max_complexity,
                complexity_metric,
                distance,
                penalties,
                fastga_contigs,
            } = tp;
            setup_logger(common.verbose);

            let banded = !no_banded;

            // Set the thread pool size before any rayon use
            rayon::ThreadPoolBuilder::new()
                .num_threads(common.threads)
                .build_global()?;

            // Optional FASTGA contig table: enables decoding contig-aware fastga tracepoints by
            // placing the trace grid contig-locally (query_start - contig_start), matching the encoder.
            let contig_table: HashMap<String, Vec<(usize, usize)>> = match &fastga_contigs {
                Some(path) => match load_contig_table(path) {
                    Ok(t) => {
                        info!("Loaded fastga contig table: {} scaffolds from {}", t.len(), path);
                        t
                    }
                    Err(e) => {
                        error!("Failed to read --fastga-contigs {}: {}", path, e);
                        std::process::exit(1);
                    }
                },
                None => HashMap::new(),
            };

            // Detect input format
            let is_binary = tpa::is_tpa_file(&common.paf).unwrap_or(false);

            // Determine if we're using fastga based on tp_type
            let is_fastga = matches!(tp_type, TracepointType::Fastga | TracepointType::FastgaNoDiff);

            // Validate that complexity-metric is not used with fastga
            if is_fastga && complexity_metric.is_some() {
                error!("--complexity-metric cannot be used with --type fastga");
                error!("FastGA uses its own segmentation algorithm based on trace spacing");
                std::process::exit(1);
            }

            // Default to EditDistance for non-fastga if not specified
            let complexity_metric = complexity_metric.unwrap_or(ComplexityMetric::EditDistance);

            let banded_max_complexity = if banded {
                if is_fastga {
                    // FastGA uses per-segment edit distance for banding; no global max-complexity needed
                    Some(0)
                } else {
                    Some(max_complexity.unwrap_or_else(|| {
                        error!("Banded alignment requires specifying --max-complexity for non-fastga types");
                        std::process::exit(1);
                    }))
                }
            } else {
                if let Some(value) = max_complexity {
                    warn!(
                        "Ignoring --max-complexity={} because --no-banded was requested",
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
                    Some(penalties.unwrap_or_else(|| default.to_string()))
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

            // Adaptive memory mode: only with DB-TP (diagonal-distance) + gap-affine2p + mc>=64
            let adaptive = memory_mode == MemoryModeChoice::Adaptive;
            let adaptive_gate = adaptive
                && matches!(complexity_metric, ComplexityMetric::DiagonalDistance)
                && matches!(distance_mode, Distance::GapAffine2p { .. })
                && max_complexity.unwrap_or(32) >= 64;
            if adaptive {
                info!(
                    "memory-mode=adaptive: per-record routing {} (needs DB-TP + gap-affine2p + mc>=64); threshold max_seg*div > {}",
                    if adaptive_gate { "ENABLED" } else { "OFF (regime not eligible -> all high)" },
                    adaptive_threshold
                );
            }

            let penalties_summary = penalties_value.as_deref().unwrap_or("n/a");

            info!(
                "Converting {} tracepoints to CIGAR (complexity-metric={}, distance={}, penalties={}, banded={}{})",
                tp_type.as_str(),
                complexity_metric,
                distance,
                penalties_summary,
                if banded { "enabled" } else { "disabled" },
                banded_max_complexity.map(|mc| format!(", max-complexity={}", mc)).unwrap_or_default()
            );

            // Validate and apply conditional defaults
            let trace_spacing = if is_fastga {
                trace_spacing.or(max_complexity).unwrap_or(100)
            } else {
                if trace_spacing.is_some() {
                    error!("--trace-spacing should only be used with --type fastga");
                    std::process::exit(1);
                }
                0 // Not used when not fastga
            };

            if is_binary {
                // Read directly from binary format using streaming parallelization
                info!("Reading from binary PAF format...");
                let mut reader = match tpa::TpaReader::new(&common.paf) {
                    Ok(r) => r,
                    Err(e) => {
                        error!("Failed to open binary PAF: {}", e);
                        std::process::exit(1);
                    }
                };

                let memory_mode_wfa2 = memory_mode.to_lib_wfa2();

                // Decode chunks in parallel (per-thread aligner) but emit in input order.
                let stdout = io::stdout();
                let mut writer = BufWriter::new(stdout.lock());
                let mut it = reader.iter_records()?;
                loop {
                    // collect::<io::Result<Vec<_>>>()? propagates a mid-stream read error instead
                    // of silently truncating (no more map_while(Result::ok)).
                    let chunk: Vec<_> = it.by_ref().take(2048).collect::<io::Result<Vec<_>>>()?;
                    if chunk.is_empty() {
                        break;
                    }
                    let lines: Vec<String> = chunk
                        .par_iter()
                        .map_init(
                            || {
                                (
                                    make_aligners(
                                        &distance_mode,
                                        &memory_mode_wfa2,
                                        adaptive,
                                        adaptive_gate,
                                        adaptive_threshold,
                                    ),
                                    Vec::<u8>::new(),
                                    Vec::<u8>::new(),
                                )
                            },
                            |(aligners, query_buf, target_buf), record| {
                                process_decompress_record(
                                    record,
                                    sequence_index.as_ref(),
                                    &distance_mode,
                                    is_fastga,
                                    trace_spacing,
                                    &complexity_metric,
                                    banded,
                                    banded_max_complexity,
                                    keep_old_stats,
                                    aligners,
                                    query_buf,
                                    target_buf,
                                    &contig_table,
                                )
                            },
                        )
                        .collect();
                    for l in &lines {
                        writeln!(writer, "{}", l)?;
                    }
                }
                writer.flush()?;
            } else {
                // Read from text PAF; decode chunks in parallel but emit in input order.
                let paf_reader = get_paf_reader(&common.paf)?;
                let memory_mode_wfa2 = memory_mode.to_lib_wfa2();
                let stdout = io::stdout();
                let mut writer = BufWriter::new(stdout.lock());
                // Keep Err lines so a read failure propagates instead of silently truncating.
                let mut lines_it = paf_reader.lines().filter(|r| match r {
                    Ok(line) => !line.trim().is_empty() && !line.starts_with('#'),
                    Err(_) => true,
                });
                loop {
                    let chunk: Vec<String> =
                        lines_it.by_ref().take(2048).collect::<io::Result<Vec<_>>>()?;
                    if chunk.is_empty() {
                        break;
                    }
                    let out: Vec<Option<String>> = chunk
                        .par_iter()
                        .map_init(
                            || {
                                (
                                    make_aligners(
                                        &distance_mode,
                                        &memory_mode_wfa2,
                                        adaptive,
                                        adaptive_gate,
                                        adaptive_threshold,
                                    ),
                                    Vec::<u8>::new(),
                                    Vec::<u8>::new(),
                                )
                            },
                            |(aligners, query_buf, target_buf), line| {
                                process_decompress_line(
                                    line,
                                    &tp_type,
                                    sequence_index.as_ref(),
                                    &distance_mode,
                                    is_fastga,
                                    trace_spacing,
                                    &complexity_metric,
                                    banded,
                                    banded_max_complexity,
                                    keep_old_stats,
                                    aligners,
                                    query_buf,
                                    target_buf,
                                    &contig_table,
                                )
                            },
                        )
                        .collect();
                    for l in out.into_iter().flatten() {
                        writeln!(writer, "{}", l)?;
                    }
                }
                writer.flush()?;
            }
        }
        Args::Compress {
            input,
            output,
            tp,
            strategy_str,
            all_records,
            threads,
            verbose,
        } => {
            let TracepointOpts {
                tp_type,
                max_complexity,
                complexity_metric,
                distance,
                penalties,
                fastga_contigs,
            } = tp;
            setup_logger(verbose);

            // Set the thread pool size before any rayon use
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()?;

            let is_fastga = matches!(tp_type, TracepointType::Fastga | TracepointType::FastgaNoDiff);
            let max_complexity = max_complexity.unwrap_or(if is_fastga { 100 } else { 32 });

            // Validate complexity_metric usage
            if is_fastga && complexity_metric.is_some() {
                error!("--complexity-metric cannot be used with --type fastga");
                error!("FastGA uses its own segmentation algorithm based on trace spacing");
                std::process::exit(1);
            }

            // Default to EditDistance for non-fastga if not specified
            let complexity_metric = if is_fastga {
                ComplexityMetric::EditDistance // Dummy value, not used
            } else {
                complexity_metric.unwrap_or(ComplexityMetric::EditDistance)
            };

            // Optional FASTGA contig table: enables contig-aware per-segment fastga compression
            // (each cg:Z record splits into per-segment records, like `encode --fastga-contigs`).
            let contig_table = match &fastga_contigs {
                Some(path) => {
                    if !matches!(tp_type, TracepointType::Fastga) {
                        warn!("--fastga-contigs is only used with --type fastga; ignoring");
                        None
                    } else {
                        match load_contig_table(path) {
                            Ok(t) => {
                                info!("Loaded fastga contig table: {} scaffolds from {}", t.len(), path);
                                Some(t)
                            }
                            Err(e) => {
                                error!("Failed to load --fastga-contigs '{}': {}", path, e);
                                std::process::exit(1);
                            }
                        }
                    }
                }
                None => None,
            };

            if all_records {
                info!("All-records mode enabled: header/string table plain, records in BGZIP");
            }

            let tpa_distance = match parse_distance(distance, penalties.as_deref()) {
                Ok(d) => d,
                Err(e) => {
                    error!("Invalid distance parameters: {}", e);
                    std::process::exit(1);
                }
            };

            // Auto-detect input format by checking first non-empty line
            let mut has_cigar = false;
            let mut has_tracepoints = false;

            let reader = get_paf_reader(&input)?;
            for line in reader.lines().take(10).flatten() {
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }
                has_cigar = line.contains("\tcg:Z:");
                has_tracepoints = line.contains("\ttp:Z:");
                break;
            }

            // Log what was detected (informational only - TPA auto-detects per row)
            if has_cigar && has_tracepoints {
                info!("Input PAF has both cg:Z: and tp:Z: tags; will prefer tracepoints");
            } else if has_cigar {
                info!("Detected CIGAR tags (cg:Z:); will convert to tracepoints");
            } else if has_tracepoints {
                info!("Detected tracepoint tags (tp:Z:); compressing directly");
            } else {
                warn!(
                    "First line has neither cg:Z: nor tp:Z: tags; TPA will error per row if needed"
                );
            }

            // Parse strategy string: "single" or "first;second" for dual mode
            let config = if strategy_str.contains(';') {
                // Dual strategy mode: "first;second"
                let parts: Vec<&str> = strategy_str.splitn(2, ';').collect();
                let first_str = parts[0];
                let second_str = parts[1];

                let (first_strategy, first_layer) =
                    match tpa::CompressionStrategy::from_str_with_layer(first_str) {
                        Ok((s, l)) => (s, l),
                        Err(e) => {
                            error!("Invalid first strategy '{}': {}", first_str, e);
                            std::process::exit(1);
                        }
                    };

                let (second_strategy, second_layer) =
                    match tpa::CompressionStrategy::from_str_with_layer(second_str) {
                        Ok((s, l)) => (s, l),
                        Err(e) => {
                            error!("Invalid second strategy '{}': {}", second_str, e);
                            std::process::exit(1);
                        }
                    };

                info!(
                    "Dual strategy mode: {} ({:?}) / {} ({:?})",
                    first_str, first_layer, second_str, second_layer
                );

                let mut cfg = tpa::CompressionConfig::new()
                    .dual_strategy(first_strategy, second_strategy)
                    .dual_layer(first_layer, second_layer)
                    .tp_type(tp_type)
                    .max_complexity(max_complexity)
                    .complexity_metric(complexity_metric)
                    .distance(tpa_distance)
                    .contig_table(contig_table);

                if all_records {
                    cfg = cfg.all_records();
                }
                cfg
            } else {
                // Single strategy mode (used for both streams)
                let (strategy, layer) =
                    match tpa::CompressionStrategy::from_str_with_layer(&strategy_str) {
                        Ok((s, l)) => (s, l),
                        Err(e) => {
                            error!("Invalid strategy '{}': {}", strategy_str, e);
                            std::process::exit(1);
                        }
                    };

                let mut cfg = tpa::CompressionConfig::new()
                    .strategy(strategy)
                    .layer(layer)
                    .tp_type(tp_type)
                    .max_complexity(max_complexity)
                    .complexity_metric(complexity_metric)
                    .distance(tpa_distance)
                    .contig_table(contig_table);

                if all_records {
                    cfg = cfg.all_records();
                }
                cfg
            };

            if let Err(e) = tpa::paf_to_tpa(&input, &output, config) {
                error!("Compression failed: {}", e);
                std::process::exit(1);
            }

            info!("Compressed {} to {}", input, output);
        }
        Args::Decompress {
            input,
            output,
            decode,
            sequence_files,
            sequence_list,
            fastga_contigs,
            keep_old_stats,
            memory_mode,
            adaptive_threshold,
            verbose,
        } => {
            setup_logger(verbose);

            if decode {
                // Decompress and decode: binary → tracepoints → CIGAR
                info!("Decompressing and decoding binary PAF to CIGAR...");

                // Validate sequence files
                let sequence_paths = match collect_sequence_paths(sequence_files, sequence_list) {
                    Ok(paths) if !paths.is_empty() => paths,
                    Ok(_) => {
                        error!("--decode requires sequence files via --sequence-files or --sequence-list");
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

                // Read binary PAF
                let mut reader = match tpa::TpaReader::new(&input) {
                    Ok(r) => r,
                    Err(e) => {
                        error!("Failed to open binary PAF: {}", e);
                        std::process::exit(1);
                    }
                };

                // Get header metadata
                let (distance_mode, tp_type, complexity_metric, max_complexity) = {
                    let header = reader.header();
                    (
                        header.distance(),
                        header.tp_type(),
                        header.complexity_metric(),
                        header.max_complexity(),
                    )
                };

                // Process records using streaming parallelization
                let is_fastga = matches!(tp_type, TracepointType::Fastga | TracepointType::FastgaNoDiff);
                let trace_spacing = if is_fastga { max_complexity } else { 0 };
                let memory_mode_wfa2 = memory_mode.to_lib_wfa2();

                // Adaptive memory mode: gate to DB-TP + gap-affine2p + mc>=64.
                let adaptive = memory_mode == MemoryModeChoice::Adaptive;
                let adaptive_gate = adaptive
                    && matches!(complexity_metric, ComplexityMetric::DiagonalDistance)
                    && matches!(distance_mode, Distance::GapAffine2p { .. })
                    && max_complexity >= 64;
                if adaptive {
                    info!(
                        "memory-mode=adaptive: per-record routing {} (needs DB-TP + gap-affine2p + mc>=64); threshold max_seg*div > {}",
                        if adaptive_gate { "ENABLED" } else { "OFF (regime not eligible -> all high)" },
                        adaptive_threshold
                    );
                }

                // Contig table is required to reconstruct contig-aware FASTGA TPA correctly (the
                // trace grid is contig-local). Without it, N-gap-split records decode at the wrong
                // offset.
                let contig_table: HashMap<String, Vec<(usize, usize)>> = match &fastga_contigs {
                    Some(path) => load_contig_table(path).unwrap_or_else(|e| {
                        error!("Failed to read --fastga-contigs {}: {}", path, e);
                        std::process::exit(1);
                    }),
                    None => HashMap::new(),
                };
                if is_fastga && contig_table.is_empty() {
                    warn!("Decoding contig-aware FASTGA TPA without --fastga-contigs; reconstruction may be incorrect at N-gap boundaries");
                }

                // Decode chunks in parallel (per-thread aligner) but emit in input order, to the
                // selected output (file or stdout) — not unordered stdout from worker threads.
                let mut writer: Box<dyn Write> = if output == "-" {
                    Box::new(BufWriter::new(io::stdout()))
                } else {
                    Box::new(BufWriter::new(File::create(&output)?))
                };
                let mut it = reader.iter_records()?;
                loop {
                    // collect::<io::Result<Vec<_>>>()? propagates a mid-stream read error instead
                    // of silently truncating (no more map_while(Result::ok)).
                    let chunk: Vec<_> = it.by_ref().take(2048).collect::<io::Result<Vec<_>>>()?;
                    if chunk.is_empty() {
                        break;
                    }
                    let lines: Vec<String> = chunk
                        .par_iter()
                        .map_init(
                            || {
                                (
                                    make_aligners(
                                        &distance_mode,
                                        &memory_mode_wfa2,
                                        adaptive,
                                        adaptive_gate,
                                        adaptive_threshold,
                                    ),
                                    Vec::<u8>::new(),
                                    Vec::<u8>::new(),
                                )
                            },
                            |(aligners, query_buf, target_buf), record| {
                                process_decompress_record(
                                    record,
                                    &sequence_index,
                                    &distance_mode,
                                    is_fastga,
                                    trace_spacing,
                                    &complexity_metric,
                                    false, // banded
                                    None,  // banded_max_complexity
                                    keep_old_stats,
                                    aligners,
                                    query_buf,
                                    target_buf,
                                    &contig_table,
                                )
                            },
                        )
                        .collect();
                    for l in &lines {
                        writeln!(writer, "{}", l)?;
                    }
                }
                writer.flush()?;

                if output == "-" {
                    info!("Decompressed and decoded {} to stdout", input);
                } else {
                    info!("Decompressed and decoded {} to {}", input, output);
                }
            } else {
                // Decompress only: binary → tracepoints
                info!("Decompressing binary PAF to text format with tracepoints...");

                if let Err(e) = tpa::tpa_to_paf(&input, &output) {
                    error!("Decompression failed: {}", e);
                    std::process::exit(1);
                }

                if output == "-" {
                    info!("Decompressed {} to stdout", input);
                } else {
                    info!("Decompressed {} to {}", input, output);
                }
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

            // Set the thread pool size before any rayon use
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()?;

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

                // Process using streaming parallelization
                paf_reader
                    .lines()
                    .map_while(Result::ok)
                    .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
                    .par_bridge()
                    .for_each(|line| {
                        process_debug_line(
                            &line,
                            sequence_index.as_ref(),
                            mismatch,
                            gap_open1,
                            gap_ext1,
                            gap_open2,
                            gap_ext2,
                            max_complexity,
                        );
                    });

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
                let aligner = distance.create_aligner(None, None);
                let mut paf_cigar = Vec::new();
                align_sequences_wfa(
                    &query_seq[a_start..a_end],
                    &target_seq[b_start..b_end],
                    &aligner,
                    &mut paf_cigar,
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

/// Process a single PAF line for debugging (called from par_bridge())
#[cfg(debug_assertions)]
fn process_debug_line(
    line: &str,
    sequence_index: &SequenceIndex,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
    max_complexity: u32,
) {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        warn!(
            "{}",
            message_with_truncate_paf_file("Skipping malformed PAF line", line)
        );
        return;
    }

    let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
        warn!(
            "{}",
            message_with_truncate_paf_file("Skipping CIGAR-less PAF line", line)
        );
        return;
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
    let aligner = distance.create_aligner(None, None);
    let mut realn_cigar = Vec::new();
    align_sequences_wfa(&query_seq, &target_seq, &aligner, &mut realn_cigar);
    let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);
    // let paf_cigar = &realn_cigar;

    // Convert CIGAR to tracepoints using query (A) and target (B) coordinates.
    let tracepoints =
        cigar_to_tracepoints(paf_cigar, max_complexity, ComplexityMetric::EditDistance);
    let variable_tracepoints =
        cigar_to_variable_tracepoints(paf_cigar, max_complexity, ComplexityMetric::EditDistance);

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
        .any(|(a, b)| a.0 != b.0 || (b.1.is_some() && Some(a.1) != b.1))
    {
        println!("Tracepoints mismatch! {}", line);
        println!("\t         tracepoints: {:?}", tracepoints);
        println!("\tvariable_tracepoints: {:?}", variable_tracepoints);
        println!("\t     tracepoints_raw: {:?}", tracepoints_raw);
        println!("\tvariable_tracepoints_raw: {:?}", variable_tracepoints_raw);
        println!("\t tracepoints_diagonal: {:?}", tracepoints_diagonal);
        println!(
            "\t mixed_tracepoints_diagonal: {:?}",
            mixed_tracepoints_diagonal
        );
        println!(
            "\t variable_tracepoints_diagonal: {:?}",
            variable_tracepoints_diagonal
        );
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

    let (
        matches,
        mismatches,
        insertions,
        inserted_bp,
        deletions,
        deleted_bp,
        paf_gap_compressed_id,
        paf_block_id,
    ) = calculate_cigar_stats(paf_cigar);
    let (
        tracepoints_matches,
        tracepoints_mismatches,
        tracepoints_insertions,
        tracepoints_inserted_bp,
        tracepoints_deletions,
        tracepoints_deleted_bp,
        tracepoints_gap_compressed_id,
        tracepoints_block_id,
    ) = calculate_cigar_stats(&cigar_from_tracepoints);
    let (
        _variable_tracepoints_matches,
        _variable_tracepoints_mismatches,
        _variable_tracepoints_insertions,
        _variable_tracepoints_inserted_bp,
        _variable_tracepoints_deletions,
        _variable_tracepoints_deleted_bp,
        _variable_tracepoints_gap_compressed_id,
        _variable_tracepoints_block_id,
    ) = calculate_cigar_stats(&cigar_from_variable_tracepoints);
    let (
        tracepoints_raw_matches,
        tracepoints_raw_mismatches,
        tracepoints_raw_insertions,
        tracepoints_raw_inserted_bp,
        tracepoints_raw_deletions,
        tracepoints_raw_deleted_bp,
        tracepoints_raw_gap_compressed_id,
        tracepoints_raw_block_id,
    ) = calculate_cigar_stats(&cigar_from_tracepoints_raw);
    let (
        _variable_tracepoints_raw_matches,
        _variable_tracepoints_raw_mismatches,
        _variable_tracepoints_raw_insertions,
        _variable_tracepoints_raw_inserted_bp,
        _variable_tracepoints_raw_deletions,
        _variable_tracepoints_raw_deleted_bp,
        _variable_tracepoints_raw_gap_compressed_id,
        _variable_tracepoints_raw_block_id,
    ) = calculate_cigar_stats(&cigar_from_variable_tracepoints_raw);
    // Calculate stats for diagonal CIGAR reconstructions
    let (
        tracepoints_diagonal_matches,
        tracepoints_diagonal_mismatches,
        tracepoints_diagonal_insertions,
        tracepoints_diagonal_inserted_bp,
        tracepoints_diagonal_deletions,
        tracepoints_diagonal_deleted_bp,
        tracepoints_diagonal_gap_compressed_id,
        tracepoints_diagonal_block_id,
    ) = calculate_cigar_stats(&cigar_from_tracepoints_diagonal);
    let (
        _mixed_tracepoints_diagonal_matches,
        _mixed_tracepoints_diagonal_mismatches,
        _mixed_tracepoints_diagonal_insertions,
        _mixed_tracepoints_diagonal_inserted_bp,
        _mixed_tracepoints_diagonal_deletions,
        _mixed_tracepoints_diagonal_deleted_bp,
        _mixed_tracepoints_diagonal_gap_compressed_id,
        _mixed_tracepoints_diagonal_block_id,
    ) = calculate_cigar_stats(&cigar_from_mixed_tracepoints_diagonal);
    let (
        _variable_tracepoints_diagonal_matches,
        _variable_tracepoints_diagonal_mismatches,
        _variable_tracepoints_diagonal_insertions,
        _variable_tracepoints_diagonal_inserted_bp,
        _variable_tracepoints_diagonal_deletions,
        _variable_tracepoints_diagonal_deleted_bp,
        _variable_tracepoints_diagonal_gap_compressed_id,
        _variable_tracepoints_diagonal_block_id,
    ) = calculate_cigar_stats(&cigar_from_variable_tracepoints_diagonal);
    let (
        realign_matches,
        realign_mismatches,
        realign_insertions,
        realign_inserted_bp,
        realign_deletions,
        realign_deleted_bp,
        realign_gap_compressed_id,
        realign_block_id,
    ) = calculate_cigar_stats(&realn_cigar);

    let score_from_realign = compute_alignment_score_from_cigar(
        &realn_cigar,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );
    let score_from_paf = compute_alignment_score_from_cigar(
        paf_cigar, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2,
    );
    let score_from_tracepoints = compute_alignment_score_from_cigar(
        &cigar_from_tracepoints,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );
    let _score_from_variable_tracepoints = compute_alignment_score_from_cigar(
        &cigar_from_variable_tracepoints,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );
    let score_from_tracepoints_raw = compute_alignment_score_from_cigar(
        &cigar_from_tracepoints_raw,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );
    let _score_from_variable_tracepoints_raw = compute_alignment_score_from_cigar(
        &cigar_from_variable_tracepoints_raw,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );
    // Calculate scores for diagonal CIGAR reconstructions
    let score_from_tracepoints_diagonal = compute_alignment_score_from_cigar(
        &cigar_from_tracepoints_diagonal,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );
    let _score_from_mixed_tracepoints_diagonal = compute_alignment_score_from_cigar(
        &cigar_from_mixed_tracepoints_diagonal,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );
    let _score_from_variable_tracepoints_diagonal = compute_alignment_score_from_cigar(
        &cigar_from_variable_tracepoints_diagonal,
        mismatch,
        gap_open1,
        gap_ext1,
        gap_open2,
        gap_ext2,
    );

    if tracepoints.len() > tracepoints_diagonal.len()
    //if variable_tracepoints.len() != variable_tracepoints_diagonal.len()
    //if cigar_from_tracepoints != cigar_from_variable_tracepoints || (paf_cigar != cigar_from_tracepoints && score_from_paf != score_from_tracepoints) //&& paf_gap_compressed_id != tracepoints_gap_compressed_id
    {
        println!("CIGAR mismatch! {}", line);
        println!("\t seqa: {}", String::from_utf8(query_seq.clone()).unwrap());
        println!(
            "\t seqb: {}",
            String::from_utf8(target_seq.clone()).unwrap()
        );
        println!(
            "\t                      CIGAR from realign: {}",
            realn_cigar
        );
        println!("\t                          CIGAR from PAF: {}", paf_cigar);
        println!(
            "\t                  CIGAR from tracepoints: {}",
            cigar_from_tracepoints
        );
        println!(
            "\t              CIGAR from tracepoints_raw: {}",
            cigar_from_tracepoints_raw
        );
        //println!("\t         CIGAR from variable_tracepoints: {}", cigar_from_variable_tracepoints);
        println!(
            "\t     CIGAR from variable_tracepoints_raw: {}",
            cigar_from_variable_tracepoints_raw
        );
        println!(
            "\t         CIGAR from tracepoints_diagonal: {}",
            cigar_from_tracepoints_diagonal
        );
        //println!("\t   CIGAR from mixed_tracepoints_diagonal: {}", cigar_from_mixed_tracepoints_diagonal);
        //println!("\tCIGAR from variable_tracepoints_diagonal: {}", cigar_from_variable_tracepoints_diagonal);
        println!(
            "\t                      CIGAR score from realign: {}",
            score_from_realign
        );
        println!(
            "\t                          CIGAR score from PAF: {}",
            score_from_paf
        );
        println!(
            "\t                  CIGAR score from tracepoints: {}",
            score_from_tracepoints
        );
        println!(
            "\t              CIGAR score from tracepoints_raw: {}",
            score_from_tracepoints_raw
        );
        //println!("\t         CIGAR score from variable tracepoints: {}", score_from_variable_tracepoints);
        //println!("\t     CIGAR score from variable_tracepoints_raw: {}", score_from_variable_tracepoints_raw);
        println!(
            "\t         CIGAR score from tracepoints_diagonal: {}",
            score_from_tracepoints_diagonal
        );
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
        println!(
            "\t          tracepoints_diagonal: {:?}",
            tracepoints_diagonal
        );
        //println!("\t    mixed_tracepoints_diagonal: {:?}", mixed_tracepoints_diagonal);
        //println!("\t variable_tracepoints_diagonal: {:?}", variable_tracepoints_diagonal);
        println!(
            "\t                           bounds CIGAR from PAF: {:?}",
            get_cigar_diagonal_bounds(paf_cigar)
        );
        println!(
            "\t                   bounds CIGAR from tracepoints: {:?}",
            get_cigar_diagonal_bounds(&cigar_from_tracepoints)
        );
        println!(
            "\t               bounds CIGAR from tracepoints_raw: {:?}",
            get_cigar_diagonal_bounds(&cigar_from_tracepoints_raw)
        );
        println!(
            "\t          bounds CIGAR from tracepoints_diagonal: {:?}",
            get_cigar_diagonal_bounds(&cigar_from_tracepoints_diagonal)
        );
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
}

/// Calculate all CIGAR stats including gap compressed identity and block identity
#[cfg(debug_assertions)]
fn calculate_cigar_stats(cigar: &str) -> (usize, usize, usize, usize, usize, usize, f32, f32) {
    let stats = CigarStats::from_cigar(cigar);
    (
        stats.matches,
        stats.mismatches,
        stats.insertions(),
        stats.inserted_bp(),
        stats.deletions(),
        stats.deleted_bp(),
        stats.gap_compressed_identity(),
        stats.block_identity(),
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
    let distance = Distance::GapAffine2p {
        mismatch,
        gap_opening1: gap_open1,
        gap_extension1: gap_ext1,
        gap_opening2: gap_open2,
        gap_extension2: gap_ext2,
    };
    calculate_alignment_score(cigar, &distance)
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

fn get_paf_reader(paf: &str) -> io::Result<Box<dyn BufRead + Send>> {
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

/// Process a single line for compression
fn process_compress_line(
    line: &str,
    tp_type: &TracepointType,
    max_complexity: u32,
    complexity_metric: &ComplexityMetric,
    distance_mode: &Distance,
    writer: &Arc<Mutex<Box<dyn Write + Send>>>,
    minimal: bool,
    contig_table: &HashMap<String, Vec<(usize, usize)>>,
) {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        warn!(
            "{}",
            message_with_truncate_paf_file("Skipping malformed PAF line", line)
        );
        return;
    }

    let Some(cg_field) = fields.iter().find(|&&s| s.starts_with("cg:Z:")) else {
        warn!(
            "{}",
            message_with_truncate_paf_file("Skipping CIGAR-less PAF line", line)
        );
        return;
    };
    let cigar = &cg_field[5..];

    // Check for existing df:i: field
    let existing_df = fields
        .iter()
        .find(|&&s| s.starts_with("df:i:"))
        .and_then(|&s| s[5..].parse::<usize>().ok());

    // Handle FastGA tracepoints with potential overflow splitting
    if matches!(tp_type, TracepointType::Fastga | TracepointType::FastgaNoDiff) {
        // FastGA path unchanged: it does its own CIGAR handling (N-gap splitting).
        let stats = CigarStats::from_cigar(cigar);
        process_fastga_with_overflow(
            &fields,
            cigar,
            *tp_type,
            max_complexity,
            stats.gap_compressed_identity(),
            stats.block_identity(),
            stats.alignment_score(distance_mode),
            existing_df,
            complexity_metric,
            writer,
            minimal,
            contig_table,
        );
    } else {
        // Standard/Mixed/Variable: parse the CIGAR once, share the ops between stats and the
        // tracepoint conversion instead of parsing the string twice.
        let ops = tracepoints::cigar_str_to_cigar_ops(cigar);
        let stats = CigarStats::from_ops(&ops);
        process_single_record(
            &fields,
            cigar,
            ops,
            tp_type,
            max_complexity,
            stats.gap_compressed_identity(),
            stats.block_identity(),
            stats.alignment_score(distance_mode),
            existing_df,
            complexity_metric,
            writer,
            minimal,
        );
    }
}

fn process_fastga_with_overflow(
    fields: &[&str],
    cigar: &str,
    tp_type: TracepointType,
    max_complexity: u32,
    gap_compressed_identity: f32,
    block_identity: f32,
    alignment_score: i32,
    existing_df: Option<usize>,
    _complexity_metric: &ComplexityMetric,
    writer: &Arc<Mutex<Box<dyn Write + Send>>>,
    minimal: bool,
    contig_table: &HashMap<String, Vec<(usize, usize)>>,
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

    // If either scaffold has a contig table, take the contig-aware path (split at query AND target
    // N-gaps, contig-local grid) to match PAFtoALN exactly. Pass an empty slice for a scaffold with
    // no table entry.
    let empty_contigs: Vec<(usize, usize)> = Vec::new();
    let qcontigs = contig_table.get(fields[0]);
    let tcontigs = contig_table.get(fields[5]);
    let segments = if qcontigs.is_some() || tcontigs.is_some() {
        cigar_to_tracepoints_fastga_with_contigs(
            cigar,
            max_complexity,
            query_start,
            query_end,
            query_len,
            target_start,
            target_end,
            target_len,
            complement,
            qcontigs.unwrap_or(&empty_contigs),
            tcontigs.unwrap_or(&empty_contigs),
        )
    } else {
        cigar_to_tracepoints_fastga(
            cigar,
            max_complexity,
            query_start,
            query_end,
            query_len,
            target_start,
            target_end,
            target_len,
            complement,
        )
    };

    // Debug: log when an alignment is split into multiple segments
    if segments.len() > 1 {
        debug!(
            "SPLIT: alignment {}:{}-{} -> {}:{}-{} ({}) split into {} segments (cause: indel overflow >200bp in tracepoint):",
            fields[0], query_start, query_end,
            fields[5], target_start, target_end,
            if complement { "-" } else { "+" },
            segments.len()
        );
        for (i, (tps, (qs, qe, ts, te))) in segments.iter().enumerate() {
            // Convert target coordinates back to original space for display
            let (display_ts, display_te) = if complement {
                (target_len - te, target_len - ts)
            } else {
                (*ts, *te)
            };
            debug!(
                "  Segment {}/{}: query {}..{} ({} bp), target {}..{} ({} bp), {} tracepoints",
                i + 1,
                segments.len(),
                qs,
                qe,
                qe - qs,
                display_ts,
                display_te,
                display_te - display_ts,
                tps.len()
            );

            // Log the gap (overflow-causing indel) between this segment and the next
            if i + 1 < segments.len() {
                let (_, (next_qs, _, next_ts, _)) = &segments[i + 1];
                let query_gap = next_qs.saturating_sub(*qe);
                let target_gap = if complement {
                    // In reversed space, target coords go in opposite direction
                    te.saturating_sub(*next_ts)
                } else {
                    next_ts.saturating_sub(*te)
                };

                let gap_type = if query_gap > 0 && target_gap == 0 {
                    format!("{}bp insertion", query_gap)
                } else if target_gap > 0 && query_gap == 0 {
                    format!("{}bp deletion", target_gap)
                } else if query_gap > 0 && target_gap > 0 {
                    format!("{}bp query gap + {}bp target gap", query_gap, target_gap)
                } else {
                    "unknown gap".to_string()
                };
                debug!("    -> split caused by: {}", gap_type);
            }
        }
    }

    // For FASTGA mode, output MULTIPLE records (one per segment) like PAFtoALN.c does
    // Each segment has its own coordinates and tracepoints
    // Gaps between segments are implicit in coordinate discontinuities

    let mut valid_segment_count = 0;
    for (tracepoints, (seg_query_start, seg_query_end, seg_target_start, seg_target_end)) in
        segments.iter()
    {
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
        // Empty tracepoints (0-length segments) -> "0,0" (fastga) / "0" (fastga-no-diff), per FASTGA behavior.
        let nodiff = matches!(tp_type, TracepointType::FastgaNoDiff);
        let tracepoints_str = if tracepoints.is_empty() {
            if nodiff { "0".to_string() } else { "0,0".to_string() }
        } else if nodiff {
            let singles: Vec<usize> = tracepoints.iter().map(|(_d, b)| *b).collect();
            format_tracepoints(&TracepointData::FastgaNoDiff(singles))
        } else {
            format_ab_tracepoints(tracepoints)
        };
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
                            new_fields.push(format!("gi:f:{:.6}", gap_compressed_identity));
                            new_fields.push(format!("bi:f:{:.6}", block_identity));

                            // Add df fields (do = df-original, preserves old value)
                            if let Some(old_df) = existing_df {
                                new_fields.push(format!("do:i:{}", old_df));
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

        let mut w = writer.lock().unwrap();
        writeln!(w, "{}", new_fields.join("\t")).unwrap();
    }
}

fn process_single_record(
    fields: &[&str],
    cigar: &str,
    ops: Vec<(usize, char)>,
    tp_type: &TracepointType,
    max_complexity: u32,
    gap_compressed_identity: f32,
    block_identity: f32,
    alignment_score: i32,
    existing_df: Option<usize>,
    complexity_metric: &ComplexityMetric,
    writer: &Arc<Mutex<Box<dyn Write + Send>>>,
    minimal: bool,
) {
    // Convert CIGAR based on tracepoint type and complexity metric
    let (tracepoints_str, df_value) = match tp_type {
        TracepointType::Standard => {
            // Reuse the ops already parsed for stats instead of re-parsing the CIGAR string.
            let tp = cigar_to_tracepoints_from_ops(&ops, max_complexity, *complexity_metric);
            (format_tracepoints(&TracepointData::Standard(tp)), None)
        }
        TracepointType::Mixed => {
            let tp = cigar_to_mixed_tracepoints(cigar, max_complexity, *complexity_metric);
            (
                format_tracepoints(&TracepointData::Mixed(tp)),
                None::<usize>,
            )
        }
        TracepointType::Variable => {
            let tp = cigar_to_variable_tracepoints(cigar, max_complexity, *complexity_metric);
            (format_tracepoints(&TracepointData::Variable(tp)), None)
        }
        TracepointType::Fastga | TracepointType::FastgaNoDiff => {
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
                        new_fields.push(format!("do:i:{}", old_df));
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

    let mut w = writer.lock().unwrap();
    writeln!(w, "{}", new_fields.join("\t")).unwrap();
}

/// Process a single PAF line for decompression (called from par_bridge())
fn process_decompress_line(
    line: &str,
    tp_type: &TracepointType,
    sequence_index: &SequenceIndex,
    distance_mode: &Distance,
    fastga: bool,
    trace_spacing: u32,
    complexity_metric: &ComplexityMetric,
    banded: bool,
    banded_max_complexity: Option<u32>,
    keep_old_stats: bool,
    aligners: &mut Aligners,
    query_buf: &mut Vec<u8>,
    target_buf: &mut Vec<u8>,
    contig_table: &HashMap<String, Vec<(usize, usize)>>,
) -> Option<String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        warn!(
            "{}",
            message_with_truncate_paf_file("Skipping malformed PAF line", line)
        );
        return None;
    }

    let Some(tp_field) = fields.iter().find(|&&s| s.starts_with("tp:Z:")) else {
        warn!(
            "{}",
            message_with_truncate_paf_file("Skipping tracepoints-less PAF line", line)
        );
        return None;
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
    if let Err(msg) =
        sequence_index.fetch_sequence_into(query_name, query_start, query_end, query_buf)
    {
        error!("{}", message_with_truncate_paf_file(&msg, line));
        std::process::exit(1);
    }
    if strand == "-" && !fastga {
        reverse_complement_in_place(query_buf);
    }

    debug!(
        "Fetching target sequence {}:{}-{} on {} strand",
        target_name,
        target_start,
        target_end,
        if strand == "-" && fastga { "-" } else { "+" }
    );
    if let Err(msg) =
        sequence_index.fetch_sequence_into(target_name, target_start, target_end, target_buf)
    {
        error!("{}", message_with_truncate_paf_file(&msg, line));
        std::process::exit(1);
    }
    if fastga && strand == "-" {
        reverse_complement_in_place(target_buf);
    }

    // Parse tracepoints based on type (preserve Fastga vs FastgaNoDiff)
    let tracepoints = tpa::parse_tracepoints(tracepoints_str, *tp_type)
        .expect("Failed to parse tracepoints");

    // Compute offsets (only used for FastGA)
    let (query_offset, target_offset, complement) = if fastga {
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
        let target_offset = if strand == "+" {
            target_len - target_end
        } else {
            target_start
        };
        // Query grid is the absolute query_start, made contig-local when a contig table is given.
        let query_offset =
            query_start - contig_start_in_table(contig_table, query_name, query_start);
        // fastga/1aln tracepoints are query-forward vs revcomp-target, so '-' strand needs reversal.
        let complement = strand == "-";
        (query_offset, target_offset, complement)
    } else {
        (0, 0, false)
    };

    let max_value = if banded {
        Some(banded_max_complexity.expect("missing max-complexity with banded"))
    } else {
        None
    };
    // Adaptive memory-mode pick (no-op unless adaptive + eligible regime). Divergence from
    // the input PAF: col10 = residue matches, col11 = alignment block length.
    let divergence = {
        let matches = fields.get(9).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let blocklen = fields.get(10).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        if blocklen > 0.0 {
            1.0 - matches / blocklen
        } else {
            0.0
        }
    };
    let aligner = aligners.select(max_segment_span(&tracepoints), divergence);
    let cigar = reconstruct_cigar_with_aligner(
        &tracepoints,
        &query_buf[..],
        &target_buf[..],
        query_offset,
        target_offset,
        *complexity_metric,
        trace_spacing,
        complement,
        aligner,
        max_value,
    );

    // Parse CIGAR once and compute all stats
    let stats = CigarStats::from_cigar(&cigar);
    let gap_compressed_identity = stats.gap_compressed_identity();
    let block_identity = stats.block_identity();
    let alignment_score = stats.alignment_score(distance_mode);
    let residue_matches = stats.matches;
    let alignment_block_length =
        stats.matches + stats.mismatches + stats.inserted_bp() + stats.deleted_bp();

    // Check for existing gi:f:, bi:f:, and sc:i: fields
    let existing_gi = fields.iter().find(|&&s| s.starts_with("gi:f:"));
    let existing_bi = fields.iter().find(|&&s| s.starts_with("bi:f:"));
    let existing_sc = fields.iter().find(|&&s| s.starts_with("sc:i:"));

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
                    // Optionally keep old values (go=gi-original, bo=bi-original, so=sc-original)
                    if keep_old_stats {
                        if let Some(old_gi) = existing_gi {
                            new_fields.push(format!("go:f:{}", &old_gi[5..]));
                        }
                        if let Some(old_bi) = existing_bi {
                            new_fields.push(format!("bo:f:{}", &old_bi[5..]));
                        }
                        if let Some(old_sc) = existing_sc {
                            new_fields.push(format!("so:i:{}", &old_sc[5..]));
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

    Some(new_fields.join("\t"))
}

/// Process a single binary record for decompression (called from par_bridge())
fn process_decompress_record(
    record: &tpa::AlignmentRecord,
    sequence_index: &SequenceIndex,
    distance_mode: &Distance,
    fastga: bool,
    trace_spacing: u32,
    complexity_metric: &ComplexityMetric,
    banded: bool,
    banded_max_complexity: Option<u32>,
    keep_old_stats: bool,
    aligners: &mut Aligners,
    query_buf: &mut Vec<u8>,
    target_buf: &mut Vec<u8>,
    contig_table: &HashMap<String, Vec<(usize, usize)>>,
) -> String {
    // Get sequence names directly from record
    let query_name = &record.query_name;
    let target_name = &record.target_name;

    debug!(
        "Fetching query sequence {}:{}-{} on {} strand",
        query_name,
        record.query_start,
        record.query_end,
        if record.strand == '-' && !fastga {
            "-"
        } else {
            "+"
        }
    );
    if let Err(msg) = sequence_index.fetch_sequence_into(
        query_name,
        record.query_start as usize,
        record.query_end as usize,
        query_buf,
    ) {
        error!("{}", msg);
        std::process::exit(1);
    }
    if record.strand == '-' && !fastga {
        reverse_complement_in_place(query_buf);
    }

    debug!(
        "Fetching target sequence {}:{}-{} on {} strand",
        target_name,
        record.target_start,
        record.target_end,
        if record.strand == '-' && fastga {
            "-"
        } else {
            "+"
        }
    );
    if let Err(msg) = sequence_index.fetch_sequence_into(
        target_name,
        record.target_start as usize,
        record.target_end as usize,
        target_buf,
    ) {
        error!("{}", msg);
        std::process::exit(1);
    }
    if fastga && record.strand == '-' {
        reverse_complement_in_place(target_buf);
    }

    // Compute offsets (only used for FastGA)
    let (query_offset, target_offset, complement) = if fastga {
        let target_len = record.target_len as usize;
        let target_offset = if record.strand == '+' {
            target_len - (record.target_end as usize)
        } else {
            record.target_start as usize
        };
        let qs = record.query_start as usize;
        // Contig-local trace grid when a contig table is provided (matches the contig-aware encoder).
        let query_offset = qs - contig_start_in_table(contig_table, query_name, qs);
        (query_offset, target_offset, record.strand == '-')
    } else {
        (0, 0, false)
    };

    let max_value = if banded {
        Some(banded_max_complexity.expect("missing max-complexity with banded"))
    } else {
        None
    };
    // Adaptive memory-mode pick: medium for divergent, big-segment records (no-op unless
    // adaptive mode + eligible regime). divergence from the input PAF identity.
    let divergence = if record.alignment_block_len > 0 {
        1.0 - (record.residue_matches as f64) / (record.alignment_block_len as f64)
    } else {
        0.0
    };
    let aligner = aligners.select(max_segment_span(&record.tracepoints), divergence);
    let cigar = reconstruct_cigar_with_aligner(
        &record.tracepoints,
        &query_buf[..],
        &target_buf[..],
        query_offset,
        target_offset,
        *complexity_metric,
        trace_spacing,
        complement,
        aligner,
        max_value,
    );

    // Parse CIGAR once and compute all stats
    let stats = CigarStats::from_cigar(&cigar);
    let gap_compressed_identity = stats.gap_compressed_identity();
    let block_identity = stats.block_identity();
    let alignment_score = stats.alignment_score(distance_mode);

    // Build output line - start with core PAF fields
    let mut new_fields = vec![
        query_name.to_string(),
        record.query_len.to_string(),
        record.query_start.to_string(),
        record.query_end.to_string(),
        record.strand.to_string(),
        target_name.to_string(),
        record.target_len.to_string(),
        record.target_start.to_string(),
        record.target_end.to_string(),
        record.residue_matches.to_string(),
        record.alignment_block_len.to_string(),
        record.mapping_quality.to_string(),
    ];

    // Add optional tags from binary record (all tags except tp:Z: which we replace with cg:Z:)
    for tag in &record.tags {
        let key_str = std::str::from_utf8(&tag.key).unwrap_or("??");

        // Skip tp tag (we're replacing it with cigar)
        if key_str == "tp" {
            // Optionally keep old gi/bi/sc as go/bo/so
            if keep_old_stats {
                // Look for these in other tags... but for now skip
            }
            continue;
        }

        // Skip gi, bi, sc tags - they will be replaced
        if key_str == "gi" || key_str == "bi" || key_str == "sc" {
            if keep_old_stats {
                let old_prefix = format!("{}old", key_str);
                let type_char = tag.tag_type as char;
                let value_str = match &tag.value {
                    tpa::TagValue::Int(v) => v.to_string(),
                    tpa::TagValue::Float(v) => v.to_string(),
                    tpa::TagValue::String(v) => v.clone(),
                };
                new_fields.push(format!("{}:{}:{}", old_prefix, type_char, value_str));
            }
            continue;
        }

        // Add other tags as-is
        let type_char = tag.tag_type as char;
        let value_str = match &tag.value {
            tpa::TagValue::Int(v) => v.to_string(),
            tpa::TagValue::Float(v) => v.to_string(),
            tpa::TagValue::String(v) => v.clone(),
        };
        new_fields.push(format!("{}:{}:{}", key_str, type_char, value_str));
    }

    // Add new identity stats and alignment score
    new_fields.push(format!("gi:f:{:.12}", gap_compressed_identity));
    new_fields.push(format!("bi:f:{:.12}", block_identity));
    new_fields.push(format!("sc:i:{}", alignment_score));

    // Replace tracepoints with CIGAR
    new_fields.push(format!("cg:Z:{}", cigar));

    new_fields.join("\t")
}

/// Combines a message with the first 9 columns of a PAF line.
fn message_with_truncate_paf_file(message: &str, line: &str) -> String {
    let truncated_line = line.split('\t').take(9).collect::<Vec<&str>>().join("\t");
    format!("{}: {} ...", message, truncated_line)
}

/// Parse comma-separated penalty values from a string
fn parse_penalty_values(
    penalties: Option<&str>,
    expected: usize,
    description: &str,
) -> Result<Vec<i32>, String> {
    let Some(penalties) = penalties else {
        return Err(format!("missing penalties (expected {})", description));
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

/// Convert DistanceChoice + penalties to lib_wfa2 Distance
fn parse_distance(choice: DistanceChoice, penalties: Option<&str>) -> Result<Distance, String> {
    match choice {
        DistanceChoice::Edit => Ok(Distance::Edit),
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

/// In-place reverse-complement (same base mapping as `reverse_complement`), so a
/// reused fetch buffer can be complemented without a fresh allocation.
fn reverse_complement_in_place(seq: &mut [u8]) {
    let comp = |c: u8| match c {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _ => b'N',
    };
    let n = seq.len();
    let mut i = 0;
    let mut j = n;
    while i < j {
        j -= 1;
        let a = comp(seq[i]);
        let b = comp(seq[j]);
        seq[i] = b;
        seq[j] = a;
        i += 1;
    }
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
