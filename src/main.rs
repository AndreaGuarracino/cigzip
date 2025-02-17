use std::error::Error;

#[derive(Debug, Clone)]
struct Path {
    abpos: i32,    // Source start
    aepos: i32,    // Source end
    bbpos: i32,    // Target start
    bepos: i32,    // Target end
    diffs: i32,    // Total differences
}

#[derive(Debug)]
struct Alignment {
    path: Path,
    tspace: u32,           // Trace spacing
    trace: Vec<u8>,        // Raw trace points array like C
    flags: u32,            // Includes strand info
}

const COMP_FLAG: u32 = 0x1;


impl Alignment {
    // Get target position at trace index
    fn trace_target_pos(&self, idx: usize) -> u32 {
        let pos_idx = idx * 2;
        self.trace[pos_idx] as u32
    }

    // Get differences at trace index 
    fn trace_diffs(&self, idx: usize) -> u8 {
        let diff_idx = idx * 2 + 1;
        self.trace[diff_idx]
    }

    // Add a trace point
    fn add_trace(&mut self, target_pos: u32, diffs: u8) {
        self.trace.push(target_pos as u8);
        self.trace.push(diffs);
    }

    // Number of trace points
    fn trace_count(&self) -> usize {
        self.trace.len() / 2
    }
}

// Parse CIGAR operations
#[derive(Debug, PartialEq)]
struct CigarOp {
    val: u32,
}

impl CigarOp {
    pub fn new(len: i32, op: char) -> Self {
        let val = match op {
            '=' => 0,
            'X' => 1,
            'I' => 2,
            'D' => 3,
            'M' => 4,
            _ => panic!("Invalid CIGAR operation: {}", op),
        };
        Self { val: (val << 29) | (len as u32) }
    }

    pub fn op(&self) -> char {
        match self.val >> 29 {
            0 => '=',
            1 => 'X',
            2 => 'I',
            3 => 'D',
            4 => 'M',
            _ => panic!("Invalid CIGAR operation: {}", self.val >> 29),
        }
    }

    pub fn len(&self) -> i32 {
        (self.val & ((1 << 29) - 1)) as i32
    }
}

fn cigar_to_trace_points(
    cigar: &str,
    query_start: i32,
    target_start: i32,
    target_end: i32,
    tspace: u32,
    is_reverse: bool,
) -> Result<Alignment, Box<dyn Error>> {
    let mut trace = Vec::new();
    let mut current_pos = target_start;
    let mut current_diffs = 0u32;
    let mut next_trace = target_start + tspace as i32;
    let mut total_diffs = 0;
    let mut curr_query_pos = query_start;

    let mut path = Path {
        abpos: query_start,
        aepos: query_start,
        bbpos: target_start,
        bepos: target_end,
        diffs: 0,
    };

    let mut len = 0i32;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            len = len * 10 + (c as i32 - '0' as i32);
            continue;
        }

        match c {
            '=' | 'M' => {
                let mut remaining = len;
                while remaining > 0 {
                    let to_next = next_trace - current_pos;
                    if to_next <= remaining {
                        trace.push(next_trace as u8);
                        trace.push(current_diffs.min(255) as u8);
                        remaining -= to_next;
                        current_pos = next_trace;
                        curr_query_pos += to_next;
                        next_trace += tspace as i32;
                        current_diffs = 0;
                    } else {
                        current_pos += remaining;
                        curr_query_pos += remaining;
                        remaining = 0;
                    }
                }
            }
            'X' => {
                total_diffs += len;
                let mut remaining = len;
                while remaining > 0 {
                    let to_next = next_trace - current_pos;
                    if to_next <= remaining {
                        current_diffs += to_next as u32;
                        trace.push(next_trace as u8);
                        trace.push(current_diffs.min(255) as u8);
                        remaining -= to_next;
                        current_pos = next_trace;
                        curr_query_pos += to_next;
                        next_trace += tspace as i32;
                        current_diffs = 0;
                    } else {
                        current_diffs += remaining as u32;
                        current_pos += remaining;
                        curr_query_pos += remaining;
                        remaining = 0;
                    }
                }
            }
            'I' => {
                total_diffs += len;
                current_diffs += len as u32;
                curr_query_pos += len;
                
                if current_pos == next_trace {
                    trace.push(next_trace as u8);
                    trace.push(current_diffs.min(255) as u8);
                    next_trace += tspace as i32;
                    current_diffs = 0;
                }
            }
            'D' => {
                total_diffs += len;
                current_diffs += len as u32;
                
                let mut remaining = len;
                while remaining > 0 {
                    let to_next = next_trace - current_pos;
                    if to_next <= remaining {
                        trace.push(next_trace as u8);
                        trace.push(current_diffs.min(255) as u8);
                        remaining -= to_next;
                        current_pos = next_trace;
                        next_trace += tspace as i32;
                        current_diffs = 0;
                    } else {
                        current_pos += remaining;
                        remaining = 0;
                    }
                }
            }
            _ => return Err("Invalid CIGAR operation".into()),
        }
        len = 0;
    }

    // Add final trace point if needed
    if current_pos > trace.last().map(|&tp| tp as i32).unwrap_or(target_start) {
        trace.push(current_pos as u8);
        trace.push(current_diffs.min(255) as u8);
    }

    path.diffs = total_diffs;
    path.aepos = curr_query_pos;
    path.bepos = current_pos;

    let flags = if is_reverse { COMP_FLAG } else { 0 };

    Ok(Alignment {
        path,
        tspace,
        trace,
        flags,
    })
}

fn trace_points_to_cigar(alignment: &Alignment) -> String {
    let mut cigar = String::new();
    
    // For reverse strand, start at target end and work backwards
    let trace_iter = if alignment.flags & COMP_FLAG != 0 {
        let len = alignment.trace_count();
        Box::new((0..len).rev()) as Box<dyn Iterator<Item = usize>>
    } else {
        Box::new(0..alignment.trace_count()) as Box<dyn Iterator<Item = usize>>
    };
    
    let mut current_pos = if alignment.flags & COMP_FLAG != 0 {
        alignment.path.bepos
    } else {
        alignment.path.bbpos
    };
    
    // Process all trace points to build CIGAR string
    for idx in trace_iter {
        let target_pos = alignment.trace_target_pos(idx) as i32;
        let diffs = alignment.trace_diffs(idx);
        
        let length = (target_pos - current_pos).abs();
        
        if length > 0 {
            if diffs == 0 {
                cigar.push_str(&format!("{}=", length));
            } else {
                cigar.push_str(&format!("{}M", length));
            }
        }
        
        current_pos = target_pos;
    }
    
    // Handle remaining length after last trace point
    let final_pos = if alignment.flags & COMP_FLAG != 0 {
        alignment.path.bbpos
    } else {
        alignment.path.bepos  
    };
    
    let remaining = (final_pos - current_pos).abs();
    if remaining > 0 {
        if alignment.trace_diffs(alignment.trace_count() - 1) == 0 {
            cigar.push_str(&format!("{}=", remaining));
        } else {
            cigar.push_str(&format!("{}M", remaining));
        }
    }

    cigar
}

fn main() -> Result<(), Box<dyn Error>> {
    // Example usage
    let cigar = "100=10X20=5I30=2D40=";
    let query_start = 0;
    let target_start = 0;
    let trace_spacing = 50;

    println!("Original CIGAR: {}", cigar);

    // Calculate target length from CIGAR
    let target_len = get_target_len_from_cigar(cigar)?;

    let alignment = cigar_to_trace_points(
        cigar, 
        query_start, 
        target_start, 
        target_start + target_len, 
        trace_spacing, 
        false
    )?;
    println!("Trace points: {:?}", alignment);

    let reconstructed_cigar = trace_points_to_cigar(&alignment);
    println!("Reconstructed CIGAR: {}", reconstructed_cigar);

    Ok(())
}

fn get_target_len_from_cigar(cigar: &str) -> Result<i32, Box<dyn Error>> {
    let mut len = 0;
    let mut curr_len = 0i32;
    
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            len = len * 10 + (c as i32 - '0' as i32);
        } else {
            match c {
                '=' | 'X' | 'D' | 'M' => curr_len += len as i32,
                'I' => (),
                _ => return Err("Invalid CIGAR operation".into()),
            }
            len = 0;
        }
    }
    Ok(curr_len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "Invalid CIGAR operation: Y")]
    fn test_invalid_cigar_op() {
        let _ = CigarOp::new(10, 'Y');
    }

    #[test]
    fn test_all_mismatches() {
        let cigar = "100X";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        
        // With 50bp trace spacing, we get 3 trace points:
        // - At position 50 with 50 differences
        // - At position 100 with 50 differences
        // - A final point at position 100 with 0 differences
        assert_eq!(alignment.trace_count(), 3);
        assert_eq!(alignment.trace_target_pos(0), 50);
        assert_eq!(alignment.trace_diffs(0), 50);
        assert_eq!(alignment.trace_target_pos(1), 100);
        assert_eq!(alignment.trace_diffs(1), 50);
        assert_eq!(alignment.trace_target_pos(2), 100);
        assert_eq!(alignment.trace_diffs(2), 0);
        
        // Verify total path properties
        assert_eq!(alignment.path.diffs, 100); // All positions are mismatches
        assert_eq!(alignment.path.bepos - alignment.path.bbpos, 100); // Total length
        assert_eq!(alignment.path.aepos - alignment.path.abpos, 100); // Query length
    }

    #[test]
    fn test_complex_pattern() {
        let cigar = "25=5X20=10I40=5D30=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        
        assert_eq!(alignment.trace_target_pos(0), 50);
        assert_eq!(alignment.trace_diffs(0), 5); // From the 5X
        assert_eq!(alignment.path.diffs, 20); // 10 mismatches + 5 inserted bases + 5 deleted bases
    }


    #[test]
    fn test_non_zero_start() {
        let cigar = "50=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 100, 200, 200 + target_len, 25, false).unwrap();
        
        assert_eq!(alignment.path.abpos, 100);
        assert_eq!(alignment.path.bbpos, 200);
        assert_eq!(alignment.trace_target_pos(0), 225);
        assert_eq!(alignment.trace_diffs(0), 0);
    }
    
    #[test]
    fn test_short_alignment() {
        let cigar = "10=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        assert_eq!(alignment.trace_count(), 1);
        assert_eq!(alignment.trace_target_pos(0), 10);
    }

    #[test]
    fn test_reverse_strand() {
        let cigar = "100=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, true).unwrap();
        
        assert!(alignment.flags & COMP_FLAG != 0);
        let reconstructed = trace_points_to_cigar(&alignment);
        assert_eq!(reconstructed, "50=50=");
    }

    #[test]
    fn test_only_indels() {
        let cigar = "5I5D";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        assert_eq!(alignment.trace_count(), 1);
        assert_eq!(alignment.trace_target_pos(0), 5);
        assert_eq!(alignment.path.diffs, 10); // 5I + 5D = 10 base differences
    }

    #[test]
    fn test_perfect_matches() {
        let cigar = "100=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        let reconstructed = trace_points_to_cigar(&alignment);
        assert_eq!(reconstructed, "50=50=");
    }

    #[test]
    fn test_cigar_op_properties() {
        let op = CigarOp::new(100, '=');
        assert_eq!(op.len(), 100);
        assert_eq!(op.op(), '=');

        let op = CigarOp::new(50, 'X');
        assert_eq!(op.len(), 50);
        assert_eq!(op.op(), 'X');
    }

    #[test]
    fn test_boundary_conditions() {
        let cases = vec![
            ("1=", 1),
            ("1000=", 1000),
            ("49=2I49=", 98),
            ("50=2I48=", 98),
            ("51=2I47=", 98),
        ];
    
        for (cigar, expected_target_len) in cases {
            let target_len = get_target_len_from_cigar(cigar).unwrap();
            let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
            
            assert_eq!(
                alignment.path.bepos - alignment.path.bbpos,
                expected_target_len,
                "Target length mismatch for CIGAR: {}", cigar
            );
    
            let original_target_consuming = get_target_len_from_cigar(cigar).unwrap();
                    
            assert_eq!(
                original_target_consuming,
                expected_target_len,
                "Original CIGAR target length doesn't match expected for: {}", cigar
            );
        }
    }
    
    #[test]
    fn test_trace_spacing_variations() {
        let cigar = "100=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let spacings = vec![10, 25, 50, 75, 100];
        
        for spacing in spacings {
            let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, spacing, false).unwrap();
            let last_idx = alignment.trace_count() - 1;
            assert_eq!(alignment.trace_target_pos(last_idx), 100);
            
            let expected_points = (100 + spacing) / spacing;
            assert_eq!(alignment.trace_count(), expected_points as usize);
        }
    }

    #[test]
    fn test_indel_at_boundary() {
        let cigar = "50=2I48=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        assert_eq!(alignment.trace_count(), 2);
        assert_eq!(alignment.trace_target_pos(0), 50);
        assert_eq!(alignment.trace_diffs(0), 0);
        assert_eq!(alignment.trace_target_pos(1), 98);
        assert_eq!(alignment.trace_diffs(1), 2);
    }

    #[test]
    fn test_mixed_operations() {
        let cigar = "25=10X15=5I20=5D25=";
        let target_len = get_target_len_from_cigar(cigar).unwrap();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 30, false).unwrap();
        
        // 10X + 5I + 5D = 20 differences total
        assert_eq!(alignment.path.diffs, 20);
        
        // First point at 30 contains 10 mismatches (X)
        assert_eq!(alignment.trace_target_pos(0), 30);
        assert_eq!(alignment.trace_diffs(0), 5);
        
        // Second point contains the insertion (5I)
        assert_eq!(alignment.trace_target_pos(1), 60);
        assert_eq!(alignment.trace_diffs(1), 10);
        
        // Third point contains the deletion (5D)
        assert_eq!(alignment.trace_target_pos(2), 90);
        assert_eq!(alignment.trace_diffs(2), 5);
    }
}