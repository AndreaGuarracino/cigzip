use std::cmp::min;
use clap::Parser;
use std::io::{BufRead, BufReader};
use std::fs::File;
use rust_htslib::faidx::Reader as FastaReader;
use libwfa2::affine_wavefront::{AffineWavefronts, AlignmentSpan, AlignmentStatus, AlignmentScope, MemoryMode};

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

    /// Delta value for tracepoints
    #[arg(long, default_value = "100")]
    delta: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments.
    let args = Args::parse();
    let tokens: Vec<&str> = args.penalties.split(',').collect();
    if tokens.len() != 5 {
        eprintln!("Error: penalties must be provided as mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2");
        std::process::exit(1);
    }
    let mismatch: i32 = tokens[0].parse()?;
    let gap_open1: i32 = tokens[1].parse()?;
    let gap_ext1: i32 = tokens[2].parse()?;
    let gap_open2: i32 = tokens[3].parse()?;
    let gap_ext2: i32 = tokens[4].parse()?;
    eprintln!("Penalties: {},{},{},{},{}", mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);

    let delta = args.delta;

    // If both a PAF file and FASTA file are provided, process each PAF record.
    if let (Some(paf_path), Some(fasta_path)) = (args.paf, args.fasta) {
        eprintln!("Using FASTA file: {}", fasta_path);
        eprintln!("Processing PAF file: {}", paf_path);

        // Open the FASTA file using rust-htslib.
        let fasta_reader = FastaReader::from_path(&fasta_path)
            .expect("Error reading FASTA file");

        // Open the PAF file (or use stdin if "-" is provided).
        let reader: Box<dyn BufRead> = if paf_path == "-" {
            Box::new(BufReader::new(std::io::stdin()))
        } else {
            Box::new(BufReader::new(File::open(&paf_path)?))
        };

        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                eprintln!("Skipping malformed PAF line {}: {}", i + 1, line);
                continue;
            }
            // Parse mandatory PAF fields.
            let query_name = fields[0];
            let _query_len: usize = fields[1].parse().unwrap_or(0);
            let query_start: usize = fields[2].parse()?;
            let query_end: usize = fields[3].parse()?;
            let strand = fields[4];
            if strand == "-" {
                continue; // Still not supported
            }
            let target_name = fields[5];
            let _target_len: usize = fields[6].parse().unwrap_or(0);
            let target_start: usize = fields[7].parse()?;
            let target_end: usize = fields[8].parse()?;

            // Find the cg:Z: field (the CIGAR string).
            let cg_field = fields.iter().find(|&&s| s.starts_with("cg:Z:"));
            if cg_field.is_none() {
                eprintln!("Skipping line {}: no cg field found", i + 1);
                continue;
            }
            let paf_cigar = cg_field.unwrap().strip_prefix("cg:Z:").unwrap();

            // Print PAF row
            // eprintln!("{}", line);
            // eprintln!("Line {}: Query: {}:{}-{}", i + 1, query_name, query_start, query_end);
            // eprintln!("Line {}: Target: {}:{}-{}", i + 1, target_name, target_start, target_end);

            // Convert CIGAR to tracepoints using query (A) and target (B) coordinates.
            let tracepoints = cigar_to_tracepoints(paf_cigar, query_start, query_end, target_start, target_end, delta);

            // Fetch query sequence from FASTA.
            // Note: rust-htslib uses 1-based coordinates in region strings.
            let query_seq_bytes = fasta_reader.fetch_seq(query_name, query_start, query_end - 1)?;
            let query_seq = String::from_utf8(query_seq_bytes)?.to_ascii_uppercase();

            // Fetch target sequence from FASTA.
            let target_seq_bytes = fasta_reader.fetch_seq(target_name, target_start, target_end - 1)?;
            let target_seq = String::from_utf8(target_seq_bytes)?.to_ascii_uppercase();

            // eprintln!("Line {}: Query sequence: {} (len: {})", i + 1, query_seq, query_seq.len());
            // eprintln!("Line {}: Target sequence: {} (len: {})", i + 1, target_seq, target_seq.len());

            // Reconstruct the CIGAR from tracepoints.
            let recon_cigar = tracepoints_to_cigar(
                &tracepoints,
                &query_seq,
                &target_seq,
                0,
                query_seq.len(),
                0,
                target_seq.len(),
                delta,
                mismatch,
                gap_open1,
                gap_ext1,
                gap_open2,
                gap_ext2,
            );

            //let realn_cigar = align_segment_dual_gap_affine_wfa(&query_seq, &target_seq, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
            //let realn_cigar = cigar_ops_to_cigar_string(&realn_cigar);

            if paf_cigar == recon_cigar {
                //eprintln!("Line {}: Conversion successful.", i + 1);
            } else {
                eprintln!("Line {}: Conversion mismatch!", i + 1);
                //eprintln!("\t           Query seq: {:?}", query_seq);
                //eprintln!("\t          Target seq: {:?}", target_seq);
                eprintln!("\t            Tracepoints: {:?}", tracepoints);
                eprintln!("\t CIGAR from tracepoints: {}", recon_cigar);
                eprintln!("\t     CIGAR from the PAF: {}", paf_cigar);
                //eprintln!("\t CIGAR from realignment: {}", realn_cigar);
            }
        }
    } else {
        // Fallback: run default example if no PAF/FASTA provided.
        eprintln!("No PAF and FASTA provided, running default example.");

        let a_seq: String = "GAACAGAGAAATGGTGGAATTCAAATACAAAAAAACCGCAAAATTAAAAATCTTGCGGCTCTCTGAACTCATTTTCATGAGTGAATTTGGCGGAACGGACGGGACTCGAACCCGCGACCCCCTGCGTGACAGGCAGGTATTCTAACCGACTGAACTACCGCTCCGCCGTTGTGTTCCGTTGGGAACGGGCGAATATTACGGATTTGCCTCACCCTTCGTCAACGGTTTTTCTCATCTTTTGAATCGTTTGCTGCAAAAATCGCCCAAGCCGCTATTTTTAGCGCCTTTTACAGGTATTTATGCCCGCCAGAGGCAGCTTCCGCCCTTCTTCTCCACCAGATCAAGACGGGCTTCCTGAGCTGCAAGCTCTTCATCTGTCGCAAAAACAACGCGTAACTTACTTGCCTGACGTACAATGCGCTGAATTGTTGCTTCACCTTGTTGCTGCTGTGTCTCTCCTTCCATCGCAAAAGCCATCGACGTTTGACCACCGGTCATCG".to_owned();
        let b_seq: String = "ACCGAATAAATCTTTGTCTGTCAGTTGCTGGGGAGTGTTTTTATTAATATTACTTTGTGTGACTTGAGCTTTTTTTTCTTCTTTTCATCAATAGAAGATTCAGAGGCACATCCTGCAAGCAGCGAGAATAACACACATATGCAACAAAGCTTTTTTGAATTCATAATTGGACACTCCCTCGCCTGTCATAATGATTAGAATATTTAGCATTGATAACGGACTCTAATTATATACCCCACTTAAATAATAACAAACAAAATCAATTTGAAATATACATTCATTTAATCAATATAAATTTTATTCTCTGGAAATATATTTAAGCTGAATGGTTCGATATCATTATTGAATTTTAATGAATGAGACAACAAGCGAGAATCACGCACTTACGCCAGAATAAAACATTGAACAGAGAAATGGTGGAATTCAAATACAAAAAAACCGCAAAATTAAAAATCTTGCGGGCTCTCTGAACTCATTTTCATGAGTGAATTTGGCGGAAC".to_owned();

        let a_start = 0;
        let a_end = a_seq.len();
        let b_start = 0;
        let b_end = b_seq.len();

        let cigar = align_segment_dual_gap_affine_wfa(&a_seq[a_start..a_end], &b_seq[b_start..b_end],
            mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        let cigar = cigar_ops_to_cigar_string(&cigar);

        let tracepoints = cigar_to_tracepoints(&cigar, a_start, a_end, b_start, b_end, delta);
        eprintln!("Tracepoints (d, b) = {:?}", tracepoints);

        let recon_cigar = tracepoints_to_cigar(&tracepoints, &a_seq, &b_seq,
            a_start, a_end, b_start, b_end, delta, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2);
        eprintln!("     Original CIGAR: {}", cigar);
        eprintln!("Reconstructed CIGAR: {}", recon_cigar);
        assert!(cigar == recon_cigar);
    }

    Ok(())
}

/// With the inverted logic, a is consumed by insertions:
/// so =, X, and I (and M) consume A.
fn consumes_a(op: char) -> bool {
    op == '=' || op == 'X' || op == 'I' || op == 'M'
}

/// With the inverted logic, b is consumed by deletions:
/// so =, X, and D (and M) consume B.
fn consumes_b(op: char) -> bool {
    op == '=' || op == 'X' || op == 'D' || op == 'M'
}

/// Returns true if the op counts as an edit (difference).
fn is_edit(op: char) -> bool {
    op == 'X' || op == 'I' || op == 'D'
}

/// Convert a CIGAR string into tracepoints.
/// Given:
/// - `cigar`: the alignment CIGAR string in extended format.
/// - `a_start`, `a_end`: the positions on sequence A that the alignment spans.
/// - `b_start`, `b_end`: similarly for sequence B.
/// - `delta`: the spacing for tracepoints.
/// Returns a vector of (diff_count, b_bases) pairs for each A–segment.
fn cigar_to_tracepoints(
    cigar: &str,
    a_start: usize,
    a_end: usize,
    _b_start: usize,
    _b_end: usize,
    delta: usize,
) -> Vec<(usize, usize)> {
    let ops = cigar_str_to_cigar_ops(cigar);

    // Calculate the next A boundary (tracepoint) after a_start.
    let mut next_thresh = if a_start % delta == 0 {
        a_start + delta
    } else {
        ((a_start / delta) + 1) * delta
    };

    let mut a_pos = a_start;
    let mut current_diffs = 0;
    let mut current_b_bases = 0;

    let mut tracepoints = Vec::new();

    for (mut len, op) in ops {
        // For op that does not consume A (e.g. deletion now) we process it wholly.
        if !consumes_a(op) {
            if consumes_b(op) {
                current_b_bases += len;
            }
            if is_edit(op) {
                current_diffs += len;
            }
            continue;
        }
        // For op that consumes A, we may need to “split” it if it spans a tracepoint boundary.
        while len > 0 {
            let remaining_to_thresh = next_thresh.saturating_sub(a_pos);
            let step = min(len, remaining_to_thresh);

            // Determine how many bases on B this step consumes.
            let b_step = if consumes_b(op) { step } else { 0 };

            // Update counters.
            a_pos += step;
            current_b_bases += b_step;
            if is_edit(op) {
                current_diffs += step;
            }

            len -= step;

            // If we’ve reached a tracepoint boundary (i.e. a_pos == next_thresh)
            if a_pos == next_thresh {
                // Record the tracepoint for this segment.
                tracepoints.push((current_diffs, current_b_bases));
                // Reset accumulators for the next segment.
                current_diffs = 0;
                current_b_bases = 0;
                next_thresh += delta;
            }
        }
    }
    // Always record a final segment if we haven't reached a_end or if there is leftover.
    if a_pos < a_end || current_diffs != 0 || current_b_bases != 0 {
        tracepoints.push((current_diffs, current_b_bases));
    }
    tracepoints
}

/// Given tracepoints (the per-segment (d,b) pairs), re‐construct the full CIGAR string.
/// To “fill in” the alignment for each segment we use a global alignment with affine gap penalties
/// on the corresponding sub–sequences. The A boundaries are determined by a_start, a_end and delta;
/// for each segment the B–interval length is taken from the tracepoint b value.
///
/// Note: here we ignore the d value (edit count) from the tracepoint.
fn tracepoints_to_cigar(
    tracepoints: &[(usize, usize)],
    a_seq: &str,
    b_seq: &str,
    a_start: usize,
    a_end: usize,
    b_start: usize,
    _b_end: usize,
    delta: usize,
    mismatch: i32,
    gap_open_i: i32,
    gap_extend_i: i32,
    gap_open_d: i32,
    gap_extend_d: i32,
) -> String {
    // First, compute the number of intervals based solely on A consumption.
    // (This does not count any trailing B-only operations.)
    let consumed_intervals = ((a_end - a_start) + delta - 1) / delta;
    // If we recorded more tracepoints than that, we have a trailing segment.
    let extra = if tracepoints.len() > consumed_intervals { 1 } else { 0 };
    //let total_intervals = consumed_intervals + extra;

    // Now compute the boundaries.
    // For the consumed intervals, we use the standard boundaries.
    let mut boundaries = Vec::new();
    boundaries.push(a_start);
    let mut next = if a_start % delta == 0 {
        a_start + delta
    } else {
        ((a_start / delta) + 1) * delta
    };
    while boundaries.len() < consumed_intervals {
        boundaries.push(next);
        next += delta;
    }
    boundaries.push(a_end);
    // If there's an extra tracepoint (trailing B-only segment),
    // add an extra boundary equal to a_end so that the final interval is zero-length in A.
    if extra == 1 {
        boundaries.push(a_end);
    }
    // Now, boundaries.len() - 1 should equal total_intervals (which equals tracepoints.len()).
    if boundaries.len() - 1 != tracepoints.len() {
        panic!(
            "Mismatch: {} intervals vs {} tracepoint pairs",
            boundaries.len() - 1,
            tracepoints.len()
        );
    }

    let mut cigar_ops = Vec::new();
    let mut current_b = b_start;
    for (i, &(_d, b_len)) in tracepoints.iter().enumerate() {
        let a_left = boundaries[i];
        let a_right = boundaries[i + 1];

        // Extract the sub–sequences.
        let a_sub = &a_seq[a_left..a_right];
        let b_sub = &b_seq[current_b..current_b + (b_len as usize)];

        // Add debug output
        // eprintln!("Segment {}: A[{}..{}] (len {}), B[{}..{}] (len {})",
        // i, a_left, a_right, a_right - a_left,
        // current_b, current_b + b_len, b_len);

        // Align the two segments using affine gap penalties.
        let mut seg_cigar = align_segment_dual_gap_affine_wfa(
            a_sub,
            b_sub,
            mismatch,
            gap_open_i,
            gap_extend_i,
            gap_open_d,
            gap_extend_d,
        );

        // Append to our overall CIGAR operations.
        cigar_ops.append(&mut seg_cigar);
        current_b += b_len as usize;
    }

    // Merge adjacent operations if needed.
    let merged = merge_cigar_ops(cigar_ops);
    cigar_ops_to_cigar_string(&merged)
}

/// Helper: Merge two vectors of CIGAR operations (merging adjacent ops of the same kind).
fn merge_cigar_ops(ops: Vec<(usize, char)>) -> Vec<(usize, char)> {
    if ops.is_empty() {
        return ops;
    }
    let mut merged = Vec::new();
    let (mut count, mut op) = ops[0];
    for &(c, o) in ops.iter().skip(1) {
        if o == op {
            count += c;
        } else {
            merged.push((count, op));
            op = o;
            count = c;
        }
    }
    merged.push((count, op));
    merged
}

/// Parse a CIGAR string into a vector of (length, op) pairs.
fn cigar_str_to_cigar_ops(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_digit(10) {
            num.push(ch);
        } else {
            if let Ok(n) = num.parse::<usize>() {
                ops.push((n, ch));
            }
            num.clear();
        }
    }
    ops
}

fn cigar_u8_to_cigar_ops(ops: &[u8]) -> Vec<(usize, char)> {
    let mut result = Vec::new();
    let mut count = 1;
    let mut current_op = if ops[0] == 77 { '=' } else { ops[0] as char };
    
    for &byte in ops.iter().skip(1) {
        let op = if byte == 77 { '=' } else { byte as char };
        if op == current_op {
            count += 1;
        } else {
            result.push((count, current_op));
            current_op = op;
            count = 1;
        }
    }
    
    // Push the final operation
    result.push((count, current_op));
    
    result
}

/// Convert a vector of (length, op) pairs to a CIGAR string.
fn cigar_ops_to_cigar_string(ops: &[(usize, char)]) -> String {
    ops.iter()
        .map(|(len, op)| format!("{}{}", len, op))
        .collect::<Vec<_>>()
        .join("")
}

/// A simple Needleman–Wunsch alignment that returns a CIGAR string in extended format.
/// Here we use unit cost for a mismatch and gap. With inverted gap logic:
/// - Insertion ('I') now consumes an A–base (vertical move).
/// - Deletion ('D') now consumes a B–base (horizontal move).
fn align_segment(a: &str, b: &str) -> Vec<(usize, char)> {
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let n = a_bytes.len();
    let m = b_bytes.len();

    // Create DP table.
    let mut dp = vec![vec![0; m + 1]; n + 1];
    for i in 0..=n {
        dp[i][0] = i;
    }
    for j in 0..=m {
        dp[0][j] = j;
    }
    for i in 1..=n {
        for j in 1..=m {
            let cost = if a_bytes[i - 1] == b_bytes[j - 1] { 0 } else { 1 };
            dp[i][j] = min(
                min(dp[i - 1][j - 1] + cost, dp[i - 1][j] + 1),
                dp[i][j - 1] + 1,
            );
        }
    }

    // Traceback to recover operations.
    let mut i = n;
    let mut j = m;
    let mut ops_rev = Vec::new();
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = if a_bytes[i - 1] == b_bytes[j - 1] { 0 } else { 1 };
            if dp[i][j] == dp[i - 1][j - 1] + cost {
                // Use '=' for a match and 'X' for a mismatch.
                ops_rev.push(if cost == 0 { '=' } else { 'X' });
                i -= 1;
                j -= 1;
                continue;
            }
        }
        // Inverted gap logic:
        // Vertical move (i-1, j) is now an insertion (consuming A).
        if i > 0 && dp[i][j] == dp[i - 1][j] + 1 {
            ops_rev.push('I');
            i -= 1;
        }
        // Horizontal move (i, j-1) is now a deletion (consuming B).
        else if j > 0 && dp[i][j] == dp[i][j - 1] + 1 {
            ops_rev.push('D');
            j -= 1;
        }
    }
    ops_rev.reverse();

    // Compress consecutive identical operations into (length, op) pairs.
    let mut cigar = Vec::new();
    if !ops_rev.is_empty() {
        let mut count = 1;
        let mut last = ops_rev[0];
        for &op in ops_rev.iter().skip(1) {
            if op == last {
                count += 1;
            } else {
                cigar.push((count, last));
                last = op;
                count = 1;
            }
        }
        cigar.push((count, last));
    }
    cigar
}

/// align_segment_gap_affine performs global alignment using affine gap penalties.
/// It returns a vector of (length, op) pairs representing the CIGAR string.
///  
/// The scoring is:
///   - Match: 0
///   - Mismatch: mismatch
///   - Gap open: gap_open
///   - Gap extend: gap_extend
///
/// With the inverted logic:
/// - Ix (gap in B) now produces an 'I' (insertion; consumes A).
/// - Iy (gap in A) now produces a 'D' (deletion; consumes B).
fn align_segment_gap_affine(
    a: &str,
    b: &str,
    mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<(usize, char)> {
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let n = a_bytes.len();
    let m = b_bytes.len();
    let inf = i32::MAX / 2;

    // Three matrices: M (match/mismatch), Ix (gap in B), Iy (gap in A)
    let mut dp_m = vec![vec![inf; m + 1]; n + 1];
    let mut dp_ix = vec![vec![inf; m + 1]; n + 1];
    let mut dp_iy = vec![vec![inf; m + 1]; n + 1];

    #[derive(Clone, Copy, Debug)]
    enum State {
        M,
        Ix,
        Iy,
    }
    let mut ptr_m = vec![vec![None; m + 1]; n + 1];
    let mut ptr_ix = vec![vec![None; m + 1]; n + 1];
    let mut ptr_iy = vec![vec![None; m + 1]; n + 1];

    // Initialization.
    dp_m[0][0] = 0;
    dp_ix[0][0] = inf;
    dp_iy[0][0] = inf;
    ptr_m[0][0] = None;

    for i in 0..=n {
        dp_m[i][0] = i as i32;
    }
    for j in 0..=m {
        dp_m[0][j] = j as i32;
    }
    for i in 1..=n {
        if i == 1 {
            dp_ix[i][0] = gap_open;
            ptr_ix[i][0] = Some(State::M);
        } else {
            dp_ix[i][0] = dp_ix[i - 1][0] + gap_extend;
            ptr_ix[i][0] = Some(State::Ix);
        }
        dp_m[i][0] = dp_ix[i][0];
        ptr_m[i][0] = Some(State::Ix);
        dp_iy[i][0] = inf;
    }
    for j in 1..=m {
        if j == 1 {
            dp_iy[0][j] = gap_open;
            ptr_iy[0][j] = Some(State::M);
        } else {
            dp_iy[0][j] = dp_iy[0][j - 1] + gap_extend;
            ptr_iy[0][j] = Some(State::Iy);
        }
        dp_m[0][j] = dp_iy[0][j];
        ptr_m[0][j] = Some(State::Iy);
        dp_ix[0][j] = inf;
    }

    // Fill in the DP matrices.
    for i in 1..=n {
        for j in 1..=m {
            let sub_cost = if a_bytes[i - 1] == b_bytes[j - 1] { 0 } else { mismatch };

            // M: coming diagonally from M, Ix, or Iy.
            let cand_m = dp_m[i - 1][j - 1] + sub_cost;
            let cand_ix = dp_ix[i - 1][j - 1] + sub_cost;
            let cand_iy = dp_iy[i - 1][j - 1] + sub_cost;
            if cand_m <= cand_ix && cand_m <= cand_iy {
                dp_m[i][j] = cand_m;
                ptr_m[i][j] = Some(State::M);
            } else if cand_ix <= cand_iy {
                dp_m[i][j] = cand_ix;
                ptr_m[i][j] = Some(State::Ix);
            } else {
                dp_m[i][j] = cand_iy;
                ptr_m[i][j] = Some(State::Iy);
            }

            // Ix: gap in B (vertical move, now yielding an 'I').
            let cand_from_m = dp_m[i - 1][j] + gap_open + gap_extend;
            let cand_from_ix = dp_ix[i - 1][j] + gap_extend;
            if cand_from_m <= cand_from_ix {
                dp_ix[i][j] = cand_from_m;
                ptr_ix[i][j] = Some(State::M);
            } else {
                dp_ix[i][j] = cand_from_ix;
                ptr_ix[i][j] = Some(State::Ix);
            }

            // Iy: gap in A (horizontal move, now yielding a 'D').
            let cand_from_m = dp_m[i][j - 1] + gap_open + gap_extend;
            let cand_from_iy = dp_iy[i][j - 1] + gap_extend;
            if cand_from_m <= cand_from_iy {
                dp_iy[i][j] = cand_from_m;
                ptr_iy[i][j] = Some(State::M);
            } else {
                dp_iy[i][j] = cand_from_iy;
                ptr_iy[i][j] = Some(State::Iy);
            }
        }
    }

    // Traceback: choose the best ending state.
    let mut i = n;
    let mut j = m;
    let (mut current_state, _final_score) = if dp_m[n][m] <= dp_ix[n][m] && dp_m[n][m] <= dp_iy[n][m] {
        (State::M, dp_m[n][m])
    } else if dp_ix[n][m] <= dp_iy[n][m] {
        (State::Ix, dp_ix[n][m])
    } else {
        (State::Iy, dp_iy[n][m])
    };

    let mut ops_rev = Vec::new();
    while i > 0 || j > 0 {
        if i == 0 {
            // Only deletion possible.
            ops_rev.push('D');
            j -= 1;
            continue;
        } else if j == 0 {
            // Only insertion possible.
            ops_rev.push('I');
            i -= 1;
            continue;
        }
        match current_state {
            State::M => {
                let prev = ptr_m[i][j].unwrap();
                let op = if a_bytes[i - 1] == b_bytes[j - 1] { '=' } else { 'X' };
                ops_rev.push(op);
                i -= 1;
                j -= 1;
                current_state = prev;
            }
            State::Ix => {
                let prev = ptr_ix[i][j].unwrap();
                ops_rev.push('I'); // now 'I' (vertical move, consuming A)
                i -= 1;
                current_state = prev;
            }
            State::Iy => {
                let prev = ptr_iy[i][j].unwrap();
                ops_rev.push('D'); // now 'D' (horizontal move, consuming B)
                j -= 1;
                current_state = prev;
            }
        }
    }
    ops_rev.reverse();

    // Compress consecutive operations.
    let mut cigar = Vec::new();
    if !ops_rev.is_empty() {
        let mut count = 1;
        let mut last = ops_rev[0];
        for &op in ops_rev.iter().skip(1) {
            if op == last {
                count += 1;
            } else {
                cigar.push((count, last));
                last = op;
                count = 1;
            }
        }
        cigar.push((count, last));
    }
    cigar
}

/// Enum for traceback state.
#[derive(Clone, Copy, Debug)]
enum State {
    M,  // Match/mismatch
    Ix, // Insertion (gap in B; vertical move, consumes A)
    Iy, // Deletion (gap in A; horizontal move, consumes B)
}

/// align_segment_dual_gap_affine performs global alignment with affine gap penalties,
/// where the gap penalties for insertions (Ix) and deletions (Iy) are specified separately.
/// 
/// Parameters:
/// - `a`: first sequence
/// - `b`: second sequence
/// - `mismatch`: cost for a mismatch (substitution)
/// - `gap_open_i`: gap opening penalty for an insertion (vertical gap, consuming A)
/// - `gap_extend_i`: gap extension penalty for an insertion
/// - `gap_open_d`: gap opening penalty for a deletion (horizontal gap, consuming B)
/// - `gap_extend_d`: gap extension penalty for a deletion
///
/// Returns a vector of (length, op) pairs (e.g. [(10, '='), (2, 'X'), (5, 'I'), …])
/// representing the CIGAR string for the alignment. In our inverted logic, an 'I'
/// denotes an insertion (consuming A) and a 'D' denotes a deletion (consuming B).
fn align_segment_dual_gap_affine(
    a: &str,
    b: &str,
    mismatch: i32,
    gap_open_i: i32,
    gap_extend_i: i32,
    gap_open_d: i32,
    gap_extend_d: i32,
) -> Vec<(usize, char)> {
    let a_bytes = a.as_bytes();
    let b_bytes = b.as_bytes();
    let n = a_bytes.len();
    let m = b_bytes.len();
    let inf = i32::MAX / 2;

    // Create DP matrices for three states: M (match/mismatch), Ix (insertion), Iy (deletion)
    let mut dp_m = vec![vec![inf; m + 1]; n + 1];
    let mut dp_ix = vec![vec![inf; m + 1]; n + 1];
    let mut dp_iy = vec![vec![inf; m + 1]; n + 1];

    // Pointer matrices to record traceback decisions.
    let mut ptr_m = vec![vec![None; m + 1]; n + 1];
    let mut ptr_ix = vec![vec![None; m + 1]; n + 1];
    let mut ptr_iy = vec![vec![None; m + 1]; n + 1];

    // Initialization.
    dp_m[0][0] = 0;
    dp_ix[0][0] = inf;
    dp_iy[0][0] = inf;
    ptr_m[0][0] = None;

    // Initialize first column: only vertical moves (insertions, i.e. gap in B)
    for i in 1..=n {
        if i == 1 {
            dp_ix[i][0] = gap_open_i + gap_extend_i;
            ptr_ix[i][0] = Some(State::M);
        } else {
            dp_ix[i][0] = dp_ix[i - 1][0] + gap_extend_i;
            ptr_ix[i][0] = Some(State::Ix);
        }
        dp_m[i][0] = dp_ix[i][0];
        ptr_m[i][0] = Some(State::Ix);
        dp_iy[i][0] = inf;
    }

    // Initialize first row: only horizontal moves (deletions, i.e. gap in A)
    for j in 1..=m {
        if j == 1 {
            dp_iy[0][j] = gap_open_d + gap_extend_d;
            ptr_iy[0][j] = Some(State::M);
        } else {
            dp_iy[0][j] = dp_iy[0][j - 1] + gap_extend_d;
            ptr_iy[0][j] = Some(State::Iy);
        }
        dp_m[0][j] = dp_iy[0][j];
        ptr_m[0][j] = Some(State::Iy);
        dp_ix[0][j] = inf;
    }

    // Fill in the DP matrices.
    for i in 1..=n {
        for j in 1..=m {
            let sub_cost = if a_bytes[i - 1] == b_bytes[j - 1] { 0 } else { mismatch };

            // Compute M: coming diagonally from any state.
            let cand_m = dp_m[i - 1][j - 1] + sub_cost;
            let cand_ix = dp_ix[i - 1][j - 1] + sub_cost;
            let cand_iy = dp_iy[i - 1][j - 1] + sub_cost;
            if cand_m <= cand_ix && cand_m <= cand_iy {
                dp_m[i][j] = cand_m;
                ptr_m[i][j] = Some(State::M);
            } else if cand_ix <= cand_iy {
                dp_m[i][j] = cand_ix;
                ptr_m[i][j] = Some(State::Ix);
            } else {
                dp_m[i][j] = cand_iy;
                ptr_m[i][j] = Some(State::Iy);
            }

            // Compute Ix: vertical gap (insertion, consuming A).
            let cand_from_m = dp_m[i - 1][j] + gap_open_i + gap_extend_i;
            let cand_from_ix = dp_ix[i - 1][j] + gap_extend_i;
            if cand_from_m <= cand_from_ix {
                dp_ix[i][j] = cand_from_m;
                ptr_ix[i][j] = Some(State::M);
            } else {
                dp_ix[i][j] = cand_from_ix;
                ptr_ix[i][j] = Some(State::Ix);
            }

            // Compute Iy: horizontal gap (deletion, consuming B).
            let cand_from_m = dp_m[i][j - 1] + gap_open_d + gap_extend_d;
            let cand_from_iy = dp_iy[i][j - 1] + gap_extend_d;
            if cand_from_m <= cand_from_iy {
                dp_iy[i][j] = cand_from_m;
                ptr_iy[i][j] = Some(State::M);
            } else {
                dp_iy[i][j] = cand_from_iy;
                ptr_iy[i][j] = Some(State::Iy);
            }
        }
    }

    // Traceback: choose the best ending state.
    let mut i = n;
    let mut j = m;
    let (mut current_state, _final_score) = if dp_m[n][m] <= dp_ix[n][m] && dp_m[n][m] <= dp_iy[n][m] {
        (State::M, dp_m[n][m])
    } else if dp_ix[n][m] <= dp_iy[n][m] {
        (State::Ix, dp_ix[n][m])
    } else {
        (State::Iy, dp_iy[n][m])
    };

    let mut ops_rev = Vec::new();
    while i > 0 || j > 0 {
        if i == 0 {
            // Must be a deletion.
            ops_rev.push('D');
            j -= 1;
            continue;
        } else if j == 0 {
            // Must be an insertion.
            ops_rev.push('I');
            i -= 1;
            continue;
        }
        match current_state {
            State::M => {
                let prev = ptr_m[i][j].unwrap();
                // Diagonal move: match if bases are equal, else mismatch.
                let op = if a_bytes[i - 1] == b_bytes[j - 1] { '=' } else { 'X' };
                ops_rev.push(op);
                i -= 1;
                j -= 1;
                current_state = prev;
            }
            State::Ix => {
                let prev = ptr_ix[i][j].unwrap();
                ops_rev.push('I'); // Insertion (vertical move, consuming A)
                i -= 1;
                current_state = prev;
            }
            State::Iy => {
                let prev = ptr_iy[i][j].unwrap();
                ops_rev.push('D'); // Deletion (horizontal move, consuming B)
                j -= 1;
                current_state = prev;
            }
        }
    }
    ops_rev.reverse();

    // Compress consecutive operations.
    let mut cigar = Vec::new();
    if !ops_rev.is_empty() {
        let mut count = 1;
        let mut last = ops_rev[0];
        for &op in ops_rev.iter().skip(1) {
            if op == last {
                count += 1;
            } else {
                cigar.push((count, last));
                last = op;
                count = 1;
            }
        }
        cigar.push((count, last));
    }
    cigar
}

fn align_segment_dual_gap_affine_wfa(
    a: &str, 
    b: &str,
    mismatch: i32,
    gap_open_i: i32,
    gap_extend_i: i32, 
    gap_open_d: i32,
    gap_extend_d: i32,
) -> Vec<(usize, char)> {
    // Create aligner and configure settings
    let mut aligner = AffineWavefronts::default();
    
    // Set alignment scope to compute full alignment
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    
    // Set end-to-end alignment mode
    aligner.set_alignment_span(AlignmentSpan::End2End);
    
    // Set memory mode to High for best accuracy
    aligner.set_memory_mode(MemoryMode::High);
    
    // Set the penalties
    aligner.set_penalties(0, mismatch, gap_open_i, gap_extend_i);//, gap_open_d, gap_extend_d);
    // unsafe {
    //     let wf_aligner = aligner.aligner_mut();
    //     (*wf_aligner).penalties.match_ = 0;  // Match score
    //     (*wf_aligner).penalties.mismatch = mismatch;
    //     (*wf_aligner).penalties.gap_opening1 = gap_open_i;
    //     (*wf_aligner).penalties.gap_extension1 = gap_extend_i;
    //     (*wf_aligner).penalties.gap_opening2 = gap_open_d;
    //     (*wf_aligner).penalties.gap_extension2 = gap_extend_d;
    // }

    // Do the alignment
    let status = aligner.align(b.as_bytes(), a.as_bytes());

    match status {
        AlignmentStatus::Completed => {
            cigar_u8_to_cigar_ops(aligner.cigar())
        },
        s => {
            eprintln!("Alignment failed with status: {:?}", s);
            Vec::new()
        }
    }
}