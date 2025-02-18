use std::cmp::{min};

/// Parse a CIGAR string (e.g. "43=5X10I2D") into a vector of (length, op) pairs.
fn parse_cigar(cigar: &str) -> Vec<(usize, char)> {
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

/// Returns true if the given op consumes an A–base.
fn consumes_a(op: char) -> bool {
    // In our extended CIGAR:
    // =, X, and D consume A. (I does not.)
    op == '=' || op == 'X' || op == 'D' || op == 'M'
}

/// Returns true if the given op consumes a B–base.
fn consumes_b(op: char) -> bool {
    // In our extended CIGAR:
    // =, X, and I consume B. (D does not.)
    op == '=' || op == 'X' || op == 'I' || op == 'M'
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
    b_start: usize,
    _b_end: usize, // not used directly here
    delta: usize,
) -> Vec<(u8, u8)> {
    let ops = parse_cigar(cigar);

    // Calculate the next A boundary (tracepoint) after a_start.
    let mut next_thresh = if a_start % delta == 0 {
        a_start + delta
    } else {
        ((a_start / delta) + 1) * delta
    };

    let mut a_pos = a_start;
    //let mut b_pos = b_start;
    let mut current_diffs = 0;
    let mut current_b_bases = 0;

    let mut tracepoints = Vec::new();

    for (mut len, op) in ops {
        // For op that does not consume A (e.g. insertion) we can process it wholly.
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
            //b_pos += b_step;
            current_b_bases += b_step;
            if is_edit(op) {
                current_diffs += step;
            }

            len -= step;

            // If we’ve reached a tracepoint boundary (i.e. a_pos == next_thresh)
            if a_pos == next_thresh {
                // Push the tracepoint for this segment.
                // We assume that d and b fit in u8 (as in the original design).
                tracepoints.push((current_diffs as u8, current_b_bases as u8));
                // Reset accumulators for the next segment.
                current_diffs = 0;
                current_b_bases = 0;
                next_thresh += delta;
            }
        }
    }
    // If there is a “final” segment (when a_end is not an exact multiple of delta)
    if a_pos < a_end || !tracepoints.is_empty() && a_pos != next_thresh - delta {
        tracepoints.push((current_diffs as u8, current_b_bases as u8));
    }
    tracepoints
}

/// A simple Needleman–Wunsch alignment that returns a CIGAR string in extended format.
/// Here we use unit cost for a mismatch and gap.
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
        if i > 0 && dp[i][j] == dp[i - 1][j] + 1 {
            ops_rev.push('D'); // deletion (gap in b)
            i -= 1;
        } else if j > 0 && dp[i][j] == dp[i][j - 1] + 1 {
            ops_rev.push('I'); // insertion (extra base in b)
            j -= 1;
        }
    }
    ops_rev.reverse();

    // Now compress consecutive identical operations into (length, op) pairs.
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

/// Helper: Merge two vectors of CIGAR operations (merging adjacent ops of the same kind).
fn merge_cigar_ops(mut ops: Vec<(usize, char)>) -> Vec<(usize, char)> {
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

/// Convert a vector of (length, op) pairs to a CIGAR string.
fn cigar_vec_to_string(ops: &[(usize, char)]) -> String {
    ops.iter().map(|(len, op)| format!("{}{}", len, op)).collect::<Vec<_>>().join("")
}

/// Given tracepoints (the per-segment (d,b) pairs), re‐construct the full CIGAR string.
/// To “fill in” the alignment for each segment we use a simple global alignment (NW)
/// on the corresponding sub–sequences. The A boundaries are determined by a_start, a_end and delta;
/// for each segment the B–interval length is taken from the tracepoint b value.
///
/// Note: here we ignore the d value (edit count) from the tracepoint.
fn tracepoints_to_cigar(
    tracepoints: &[(u8, u8)],
    a_seq: &str,
    b_seq: &str,
    a_start: usize,
    a_end: usize,
    b_start: usize,
    _b_end: usize, // not used directly
    delta: usize,
) -> String {
    // Determine the A boundaries.
    let mut boundaries = Vec::new();
    boundaries.push(a_start);
    let mut next = if a_start % delta == 0 {
        a_start + delta
    } else {
        ((a_start / delta) + 1) * delta
    };
    while next < a_end {
        boundaries.push(next);
        next += delta;
    }
    boundaries.push(a_end);

    // We expect the number of intervals to equal the number of tracepoint pairs.
    if boundaries.len() - 1 != tracepoints.len() {
        panic!("Mismatch: {} intervals vs {} tracepoint pairs", boundaries.len()-1, tracepoints.len());
    }

    let mut cigar_ops = Vec::new();
    let mut current_b = b_start;
    for (i, &(d, b_len)) in tracepoints.iter().enumerate() {
        let a_left = boundaries[i];
        let a_right = boundaries[i + 1];

        // Extract the sub–sequences.
        // (Assume that a_seq and b_seq are indexed by the absolute positions.)
        let a_sub = &a_seq[a_left..a_right];
        let b_sub = &b_seq[current_b..current_b + (b_len as usize)];

        // Align the two segments.
        let mut seg_cigar = align_segment(a_sub, b_sub);
        // (Optionally one could check that the edit distance equals d.)

        // Append to our overall CIGAR operations.
        cigar_ops.append(&mut seg_cigar);
        current_b += b_len as usize;
    }

    // Merge adjacent operations if needed.
    let merged = merge_cigar_ops(cigar_ops);
    cigar_vec_to_string(&merged)
}

/// ---
///
/// The code below is an example “driver” that shows how one might call these functions.
///
/// (In a real application you would have proper sequence data and alignment coordinates.)
fn main() {
    let a_seq: String = "ACCGTGAAAGGTGTTGCTGCAGGCAAGGTCAACATTCCGGTTGTATCCGGTAATGGTGAACTTGCTGTGGTTGCAGAAATCACCGTCACCGAAGTTAATCCGGAGAGTCAGCGATGTTCCTGAAAACCGATCATTTGAATATAACGGTGTGAGCGTCACGCTTTCTGAACTGTCAGCCCTCAGCGAATTGAGCATCTCGCCCAGCTGAAACGACAGGCAGAACAGGCGGGATCCAGTCTCAATCGACAGGTGAGCGTGGAAGATCTCGTCAAACCGGTGCTTTTCTGGTGGCGATGTCCCTGTGGCATAGCCATCCGCTAGAAGACAAAGATGCCGTCCATGAATGAAGCCGTTAAACAAATTGAGCAGGAAGTGCTTACCACCTGGCCCACAGAGGCAATTGCTCAGGCTGAAAATGTGGTAATGCGTCTGTCCGGTATGTCTGAGTTTGTTGTGAATGATGCCACCTGAACAGGCAGATGACGCCGGGCCAGCAGAGC".to_owned();
    let b_seq: String = "ACCGTGAAAGGTGTTTGCAGGCAAGGTCAACATTCCGGTTGTATCCGGTAATGGTGAACTTGCTGTGGTTGCAGAAATCACCGTCACCGAAGTTAATCCGGAGAGTCAGCGATGTTCCTGAAAACCGATCATTTGAATATAACGGTGTGAGCGTCACGCTTTCTGAACTGTCAGCCCTCAGCGAATTGAGCATCTCGCCCAGCTGAAACGACAGGCAGAACAGGCGGGATCCAGTCTCAATCGACAGGTGAGCGTGGAAGATCTCGTCAAACCGGTGCTTTTCTGGTGGCGATGTCCCTGTGGCATAGCCATCCGCTAGAAGACAAAGATGCCGTCCATGAATGAAGCCGTTAAACAAATTGAGCAGGAAGTGCTTACCACCTGGCCCACAGAGGCAATTGCTCAGGCTGAAAATGTGGTAATGCGTCTGTCCGGTATGTCTGAGTTTGTTGTGAATGATGCCACCTGAACAGGCAGATGACGCCGGGCCAGCAGAGC".to_owned();

    // Example: suppose we have an alignment spanning A[57,432] and B[1002,1391]
    // and the daligner computed the following CIGAR (extended) string.
    let cigar = "15=2D483=";
    let a_start = 0;
    let a_end = a_seq.len();
    let b_start = 0;
    let b_end = b_seq.len();
    let delta = 100;

    // Convert CIGAR -> tracepoints.
    let tracepoints = cigar_to_tracepoints(cigar, a_start, a_end, b_start, b_end, delta);
    println!("Tracepoints (d, b) = {:?}", tracepoints);
    // (For example, one might see: [(10,48), (25,105), (18,95), (22,101), (7,40)])

    // For the reverse conversion, we need the sequences.
    // (For this example we simply create dummy sequences of the proper lengths.)
    // For simplicity we fill with dummy letters.

    let recon_cigar = tracepoints_to_cigar(&tracepoints, &a_seq, &b_seq, a_start, a_end, b_start, b_end, delta);
    println!("Reconstructed CIGAR: {}", recon_cigar);
}

