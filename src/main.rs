use std::error::Error;

#[derive(Debug, Clone)]
struct TracePoint {
    target_pos: i64,    // Position in target sequence
    differences: u32,   // Number of differences in this interval
}

#[derive(Debug)]
struct Alignment {
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    trace_spacing: u32,
    trace_points: Vec<TracePoint>,
}

// Parse CIGAR operations
#[derive(Debug, PartialEq)]
enum CigarOp {
    Match(u32),      // M
    Equal(u32),      // =
    Mismatch(u32),   // X
    Insertion(u32),  // I
    Deletion(u32),   // D
}

impl CigarOp {
    fn length(&self) -> u32 {
        match self {
            CigarOp::Match(len) => *len,
            CigarOp::Equal(len) => *len,
            CigarOp::Mismatch(len) => *len,
            CigarOp::Insertion(len) => *len,
            CigarOp::Deletion(len) => *len,
        }
    }
}

fn parse_cigar(cigar: &str) -> Result<Vec<CigarOp>, Box<dyn Error>> {
    let mut ops = Vec::new();
    let mut current_len = 0;

    for c in cigar.chars() {
        if c.is_digit(10) {
            current_len = current_len * 10 + c.to_digit(10).unwrap();
        } else {
            let op = match c {
                'M' => CigarOp::Match(current_len),
                '=' => CigarOp::Equal(current_len),
                'X' => CigarOp::Mismatch(current_len),
                'I' => CigarOp::Insertion(current_len),
                'D' => CigarOp::Deletion(current_len),
                _ => return Err("Invalid CIGAR operation".into()),
            };
            ops.push(op);
            current_len = 0;
        }
    }
    Ok(ops)
}

fn cigar_to_trace_points(
    cigar: &str,
    query_start: i64,
    target_start: i64,
    trace_spacing: u32,
) -> Result<Alignment, Box<dyn Error>> {
    let cigar_ops = parse_cigar(cigar)?;
    let mut trace_points = Vec::new();
    
    let mut query_pos = query_start;
    let mut target_pos = target_start;
    let mut current_differences = 0;
    
    // Track next trace point position
    let mut next_trace_point = target_start + trace_spacing as i64;

    for op in &cigar_ops {
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) => {
                let len = *len as i64;
                let mut remaining_len = len;
                
                while remaining_len > 0 {
                    let distance_to_next_trace = next_trace_point - target_pos;
                    if distance_to_next_trace <= remaining_len {
                        // Add trace point
                        trace_points.push(TracePoint {
                            target_pos: next_trace_point,
                            differences: current_differences,
                        });
                        
                        // Update positions
                        target_pos += distance_to_next_trace;
                        query_pos += distance_to_next_trace;
                        remaining_len -= distance_to_next_trace;
                        
                        // Reset for next trace point
                        next_trace_point += trace_spacing as i64;
                        current_differences = 0;
                    } else {
                        // Not enough length to reach next trace point
                        target_pos += remaining_len;
                        query_pos += remaining_len;
                        remaining_len = 0;
                    }
                }
            }
            CigarOp::Mismatch(len) => {
                let len = *len as i64;
                
                // Handle mismatches that cross trace point boundaries
                let mut remaining_len = len;
                while remaining_len > 0 {
                    let distance_to_next_trace = next_trace_point - target_pos;
                    if distance_to_next_trace <= remaining_len {
                        current_differences += distance_to_next_trace as u32;
                        trace_points.push(TracePoint {
                            target_pos: next_trace_point,
                            differences: current_differences,
                        });
                        
                        target_pos += distance_to_next_trace;
                        query_pos += distance_to_next_trace;
                        remaining_len -= distance_to_next_trace;
                        
                        next_trace_point += trace_spacing as i64;
                        current_differences = 0;
                    } else {
                        current_differences += remaining_len as u32;
                        target_pos += remaining_len;
                        query_pos += remaining_len;
                        remaining_len = 0;
                    }
                }
            }
            CigarOp::Insertion(len) => {
                // Insertions affect query position but not target position
                query_pos += *len as i64;
                // Count as one difference regardless of length
                current_differences += 1;
                
                // If we're exactly at a trace point, include this difference
                if target_pos == next_trace_point {
                    trace_points.push(TracePoint {
                        target_pos,
                        differences: current_differences,
                    });
                    next_trace_point += trace_spacing as i64;
                    current_differences = 0;
                }
            }
            CigarOp::Deletion(len) => {
                let del_len = *len as i64;
                
                // Handle deletions that cross trace point boundaries
                let mut remaining_len = del_len;
                while remaining_len > 0 {
                    let distance_to_next_trace = next_trace_point - target_pos;
                    if distance_to_next_trace <= remaining_len {
                        // Count deletion as one difference in current interval
                        current_differences += 1;
                        trace_points.push(TracePoint {
                            target_pos: next_trace_point,
                            differences: current_differences,
                        });
                        
                        target_pos += distance_to_next_trace;
                        remaining_len -= distance_to_next_trace;
                        
                        next_trace_point += trace_spacing as i64;
                        current_differences = 0;
                    } else {
                        current_differences += 1;
                        target_pos += remaining_len;
                        remaining_len = 0;
                    }
                }
            }
        }
    }

    // Add final trace point if we have pending differences
    if target_pos > (trace_points.last().map(|tp| tp.target_pos).unwrap_or(target_start)) {
        trace_points.push(TracePoint {
            target_pos,
            differences: current_differences,
        });
    }

    Ok(Alignment {
        query_start,
        query_end: query_pos,
        target_start,
        target_end: target_pos,
        trace_spacing,
        trace_points,
    })
}

fn trace_points_to_cigar(alignment: &Alignment) -> String {
    let mut cigar = String::new();
    let mut current_pos = alignment.target_start;
    
    for (i, trace) in alignment.trace_points.iter().enumerate() {
        let length = trace.target_pos - current_pos;
        
        if length > 0 {
            if trace.differences == 0 {
                // All matches
                cigar.push_str(&format!("{}=", length));
            } else {
                // Mix of matches and mismatches
                // This is a simplification - in reality you'd need the sequences 
                // to know exact positions of mismatches
                cigar.push_str(&format!("{}M", length));
            }
        }
        
        current_pos = trace.target_pos;
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

    let alignment = cigar_to_trace_points(cigar, query_start, target_start, trace_spacing)?;
    println!("Trace points: {:?}", alignment);

    let reconstructed_cigar = trace_points_to_cigar(&alignment);
    println!("Reconstructed CIGAR: {}", reconstructed_cigar);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // CIGAR Parsing Tests
    #[test]
    fn test_basic_cigar_parsing() {
        let cigar = "100=10X20=5I30=2D40=";
        let ops = parse_cigar(cigar).unwrap();
        assert_eq!(ops, vec![
            CigarOp::Equal(100),
            CigarOp::Mismatch(10),
            CigarOp::Equal(20),
            CigarOp::Insertion(5),
            CigarOp::Equal(30),
            CigarOp::Deletion(2),
            CigarOp::Equal(40),
        ]);
    }

    #[test]
    fn test_invalid_cigar_op() {
        let result = parse_cigar("10Y"); // Invalid operation Y
        assert!(result.is_err());
    }

    #[test]
    fn test_all_mismatches() {
        let cigar = "100X";
        let alignment = cigar_to_trace_points(cigar, 0, 0, 50).unwrap();
        assert_eq!(alignment.trace_points.len(), 2);
        assert_eq!(alignment.trace_points[0].target_pos, 50);
        assert_eq!(alignment.trace_points[0].differences, 50);
        assert_eq!(alignment.trace_points[1].target_pos, 100);
        assert_eq!(alignment.trace_points[1].differences, 50);
    }

    #[test]
    fn test_complex_pattern() {
        let cigar = "25=5X20=10I40=5D30="; // Multiple operations crossing trace boundaries
        let alignment = cigar_to_trace_points(cigar, 0, 0, 50).unwrap();
        // Verify each trace point carefully
        let points = &alignment.trace_points;
        
        // First point should be at 50, counting differences up to that point
        assert_eq!(points[0].target_pos, 50);
        assert_eq!(points[0].differences, 5); // From the 5X
        
        // Subsequent points should continue tracking position and differences
        println!("{:?}", points);
    }

    #[test]
    fn test_non_zero_start() {
        let cigar = "50=";
        let alignment = cigar_to_trace_points(cigar, 100, 200, 25).unwrap();
        assert_eq!(alignment.query_start, 100);
        assert_eq!(alignment.target_start, 200);
        assert_eq!(alignment.trace_points[0].target_pos, 225);
    }

    #[test]
    fn test_short_alignment() {
        let cigar = "10="; // Shorter than trace spacing
        let alignment = cigar_to_trace_points(cigar, 0, 0, 50).unwrap();
        assert_eq!(alignment.trace_points.len(), 1);
        assert_eq!(alignment.trace_points[0].target_pos, 10);
    }

    #[test]
    fn test_only_indels() {
        let cigar = "5I5D";
        let alignment = cigar_to_trace_points(cigar, 0, 0, 50).unwrap();
        assert_eq!(alignment.trace_points.len(), 1);
        assert_eq!(alignment.trace_points[0].target_pos, 5); // Only the deletion affects target position
        assert_eq!(alignment.trace_points[0].differences, 2); // Both indels count as differences
    }

    // Reconstruction Tests
    #[test]
    fn test_perfect_reconstruction() {
        let original_cigar = "100=";
        let alignment = cigar_to_trace_points(original_cigar, 0, 0, 50).unwrap();
        let reconstructed = trace_points_to_cigar(&alignment);
        assert_eq!(reconstructed, "50=50=");
    }

    #[test]
    fn test_approximate_reconstruction() {
        let original_cigar = "30=10X10=";
        let alignment = cigar_to_trace_points(original_cigar, 0, 0, 25).unwrap();
        let reconstructed = trace_points_to_cigar(&alignment);
        
        // Verify total length matches
        let original_len: u32 = parse_cigar(original_cigar).unwrap().iter().map(|op| op.length()).sum();
        let reconstructed_len: u32 = parse_cigar(&reconstructed).unwrap().iter().map(|op| op.length()).sum();
        assert_eq!(original_len, reconstructed_len);
    }

    #[test]
    fn test_boundary_conditions() {
        // Test various boundary conditions
        let cases = vec![
            ("1=", 1),           // Minimal alignment
            ("1000=", 1000),     // Much longer than trace spacing
            ("49=2I49=", 98),    // Indel just before trace point
            ("50=2I48=", 98),    // Indel exactly at trace point
            ("51=2I47=", 98),    // Indel just after trace point
        ];
    
        for (cigar, expected_target_len) in cases {
            let alignment = cigar_to_trace_points(cigar, 0, 0, 50).unwrap();
            let reconstructed = trace_points_to_cigar(&alignment);
            
            // Verify target length
            assert_eq!(alignment.target_end - alignment.target_start, expected_target_len,
                      "Target length mismatch for CIGAR: {}", cigar);
            
            // Verify total alignment length preservation
            let original_target_consuming: u32 = parse_cigar(cigar)
                .unwrap()
                .iter()
                .filter_map(|op| match op {
                    CigarOp::Match(len) | CigarOp::Equal(len) | 
                    CigarOp::Mismatch(len) | CigarOp::Deletion(len) => Some(*len),
                    _ => None
                })
                .sum();
                
            let reconstructed_target_consuming: u32 = parse_cigar(&reconstructed)
                .unwrap()
                .iter()
                .filter_map(|op| match op {
                    CigarOp::Match(len) | CigarOp::Equal(len) | 
                    CigarOp::Mismatch(len) | CigarOp::Deletion(len) => Some(*len),
                    _ => None
                })
                .sum();
                
            assert_eq!(original_target_consuming as i64, expected_target_len,
                      "Original CIGAR target length doesn't match expected for: {}", cigar);
            assert_eq!(reconstructed_target_consuming as i64, expected_target_len,
                      "Reconstructed CIGAR target length doesn't match expected for: {}", cigar);
        }
    }

    #[test]
    fn test_trace_spacing_variations() {
        let cigar = "100=";
        let spacings = vec![10, 25, 50, 75, 100];
        
        for spacing in spacings {
            let alignment = cigar_to_trace_points(cigar, 0, 0, spacing).unwrap();
            assert_eq!(alignment.trace_points.last().unwrap().target_pos, 100);
            
            // Number of trace points should be ceil(100/spacing)
            let expected_points = (100 + spacing - 1) / spacing;
            assert_eq!(alignment.trace_points.len(), expected_points as usize);
        }
    }

    #[test]
    fn test_cigar_parsing() {
        let cigar = "100=10X20=5I30=2D40=";
        let ops = parse_cigar(cigar).unwrap();
        assert_eq!(ops[0], CigarOp::Equal(100));
        assert_eq!(ops[1], CigarOp::Mismatch(10));
    }

    #[test]
    fn test_trace_point_conversion() {
        let cigar = "100=10X20=";
        let alignment = cigar_to_trace_points(cigar, 0, 0, 50).unwrap();

        // For "100=10X20=" with trace_spacing=50:
        // - First trace point at pos 50 (after 50=)
        // - Second trace point at pos 100 (after another 50=)
        // - Third trace point at pos 130 (after remaining 10X20=)
        assert_eq!(alignment.trace_points.len(), 3);
        
        // Verify positions and differences
        assert_eq!(alignment.trace_points[0].target_pos, 50);
        assert_eq!(alignment.trace_points[0].differences, 0);
        
        assert_eq!(alignment.trace_points[1].target_pos, 100);
        assert_eq!(alignment.trace_points[1].differences, 0);
        
        assert_eq!(alignment.trace_points[2].target_pos, 130);
        assert_eq!(alignment.trace_points[2].differences, 10); // from the 10X
    }

    #[test]
    fn test_indel_at_boundary() {
        let cigar = "50=2I48=";
        let alignment = cigar_to_trace_points(cigar, 0, 0, 50).unwrap();
        assert_eq!(alignment.trace_points.len(), 2);
        assert_eq!(alignment.trace_points[0].target_pos, 50);
        assert_eq!(alignment.trace_points[0].differences, 0);  // First trace point has no differences
        assert_eq!(alignment.trace_points[1].target_pos, 98);
        assert_eq!(alignment.trace_points[1].differences, 1);  // Second trace point includes the insertion
    }

    #[test]
    fn test_reconstruction() {
        let original_cigar = "100=10X20=";
        let alignment = cigar_to_trace_points(original_cigar, 0, 0, 50).unwrap();
        let reconstructed = trace_points_to_cigar(&alignment);
        
        // The reconstruction won't be exactly the same due to trace point compression
        // but should have the same total length
        let expected_length: u32 = parse_cigar(original_cigar)
            .unwrap()
            .iter()
            .map(|op| op.length())
            .sum();
        
        let reconstructed_length: u32 = parse_cigar(&reconstructed)
            .unwrap()
            .iter()
            .map(|op| op.length())
            .sum();
            
        assert_eq!(expected_length, reconstructed_length);
    }
}