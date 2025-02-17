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
    is_reverse: bool,
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

fn parse_cigar(cigar: &str) -> Result<Vec<CigarOp>, Box<dyn Error>> {
    let mut ops = Vec::new();
    let mut len: i32 = 0;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            len = len*10 + (c as i32 - '0' as i32);
        } else {
            let op = CigarOp::new(len, c);
            ops.push(op);
            len = 0;
        }
    }

    Ok(ops)
}

fn cigar_to_trace_points(
    cigar: &str,
    query_start: i64,
    target_start: i64,
    target_end: i64,
    trace_spacing: u32,
    is_reverse: bool,
) -> Result<Alignment, Box<dyn Error>> {
    let cigar_ops = parse_cigar(cigar)?;
    let mut trace_points = Vec::new();
    
    // For reverse strand, adjust coordinates
    let (actual_target_start, actual_target_end) = if is_reverse {
        (target_end, target_start)  // Swap coordinates for reverse strand
    } else {
        (target_start, target_end)
    };

    let mut query_pos = query_start;
    let mut target_pos = target_start;
    let mut current_differences = 0;
    let mut next_trace_point = target_start + trace_spacing as i64;

    for op in &cigar_ops {
        match op.op() {
            '=' => {
                let len = op.len() as i64;
                let mut remaining_len = len;
                
                while remaining_len > 0 {
                    let distance_to_next_trace = next_trace_point - target_pos;
                    if distance_to_next_trace <= remaining_len {
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
                        target_pos += remaining_len;
                        query_pos += remaining_len;
                        remaining_len = 0;
                    }
                }
            }
            'X' | 'M' => {
                let len = op.len() as i64;
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
            'I' => {
                query_pos += op.len() as i64;
                current_differences += 1;
                
                if target_pos == next_trace_point {
                    trace_points.push(TracePoint {
                        target_pos,
                        differences: current_differences,
                    });
                    next_trace_point += trace_spacing as i64;
                    current_differences = 0;
                }
            }
            'D' => {
                let del_len = op.len() as i64;
                let mut remaining_len = del_len;
                
                while remaining_len > 0 {
                    let distance_to_next_trace = next_trace_point - target_pos;
                    if distance_to_next_trace <= remaining_len {
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
            _ => return Err("Invalid CIGAR operation".into()),
        }
    }

    // Add final trace point if needed
    if target_pos > (trace_points.last().map(|tp| tp.target_pos).unwrap_or(target_start)) {
        trace_points.push(TracePoint {
            target_pos,
            differences: current_differences,
        });
    }

    // For reverse strand, trace points should decrease
    if is_reverse {
        trace_points.reverse();
    }

    Ok(Alignment {
        query_start,
        query_end: query_pos,
        target_start: actual_target_start,
        target_end: actual_target_end,
        trace_spacing,
        trace_points,
        is_reverse,
    })
}

fn trace_points_to_cigar(alignment: &Alignment) -> String {
    let mut cigar = String::new();
    let mut current_pos = if alignment.is_reverse {
        alignment.target_end
    } else {
        alignment.target_start
    };

    let trace_points = if alignment.is_reverse {
        // Process trace points in reverse order
        alignment.trace_points.iter().rev().collect::<Vec<_>>()
    } else {
        alignment.trace_points.iter().collect()
    };
    
    for trace in trace_points {
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

    let alignment = cigar_to_trace_points(cigar, query_start, target_start, target_start + parse_cigar(cigar)?.iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>(), trace_spacing, false)?;
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
        assert_eq!(ops[0], CigarOp::new(100, '='));
        assert_eq!(ops[1], CigarOp::new(10, 'X'));
        assert_eq!(ops[2], CigarOp::new(20, '='));
        assert_eq!(ops[3], CigarOp::new(5, 'I'));
        assert_eq!(ops[4], CigarOp::new(30, '='));
        assert_eq!(ops[5], CigarOp::new(2, 'D'));
        assert_eq!(ops[6], CigarOp::new(40, '='));
    }

    #[test]
    #[should_panic(expected = "Invalid CIGAR operation: Y")]
    fn test_invalid_cigar_op() {
        let _ = CigarOp::new(10, 'Y');
    }

    #[test]
    fn test_all_mismatches() {
        let cigar = "100X";
        let target_len = parse_cigar(cigar).unwrap().iter()
            .map(|op| match op.op() {
                '=' | 'X' | 'D' | 'M' => op.len() as i64,
                'I' => 0,
                _ => panic!("Invalid op"),
            })
            .sum::<i64>();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        assert_eq!(alignment.trace_points.len(), 2);
        assert_eq!(alignment.trace_points[0].target_pos, 50);
        assert_eq!(alignment.trace_points[0].differences, 50);
        assert_eq!(alignment.trace_points[1].target_pos, 100);
        assert_eq!(alignment.trace_points[1].differences, 50);
    }

    #[test]
    fn test_complex_pattern() {
        let cigar = "25=5X20=10I40=5D30=";
        let target_len = parse_cigar(cigar).unwrap().iter()
            .map(|op| match op.op() {
                '=' | 'X' | 'D' | 'M' => op.len() as i64,
                'I' => 0,
                _ => panic!("Invalid op"),
            })
            .sum::<i64>();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        let points = &alignment.trace_points;
        
        assert_eq!(points[0].target_pos, 50);
        assert_eq!(points[0].differences, 5); // From the 5X
    }

    #[test]
    fn test_non_zero_start() {
        let cigar = "50=";
        let target_len = parse_cigar(cigar).unwrap().iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>();
        let alignment = cigar_to_trace_points(cigar, 100, 200, 200 + target_len, 25, false).unwrap();
        assert_eq!(alignment.query_start, 100);
        assert_eq!(alignment.target_start, 200);
        assert_eq!(alignment.trace_points[0].target_pos, 225);
    }

    #[test]
    fn test_short_alignment() {
        let cigar = "10=";
        let target_len = parse_cigar(cigar).unwrap().iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        assert_eq!(alignment.trace_points.len(), 1);
        assert_eq!(alignment.trace_points[0].target_pos, 10);
    }

    #[test]
    fn test_only_indels() {
        let cigar = "5I5D";
        let target_len = parse_cigar(cigar).unwrap().iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        assert_eq!(alignment.trace_points.len(), 1);
        assert_eq!(alignment.trace_points[0].target_pos, 5);
        assert_eq!(alignment.trace_points[0].differences, 2);
    }

    #[test]
    fn test_perfect_matches() {
        let cigar = "100=";
        let target_len = parse_cigar(cigar).unwrap().iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>();
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
            let target_len = parse_cigar(cigar).unwrap().iter()
            .map(|op| match op.op() {
                '=' | 'X' | 'D' | 'M' => op.len() as i64,
                'I' => 0,
                _ => panic!("Invalid op"),
            })
            .sum::<i64>();
            let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
            
            assert_eq!(
                alignment.target_end - alignment.target_start,
                expected_target_len,
                "Target length mismatch for CIGAR: {}", cigar
            );

            let original_target_consuming: i32 = parse_cigar(cigar)
                .unwrap()
                .iter()
                .map(|op| match op.op() {
                    '=' | 'X' | 'D' | 'M' => op.len(),
                    'I' => 0,
                    _ => panic!("Invalid op"),
                })
                .sum();
                
            assert_eq!(
                original_target_consuming as i64,
                expected_target_len,
                "Original CIGAR target length doesn't match expected for: {}", cigar
            );
        }
    }

    #[test]
    fn test_trace_spacing_variations() {
        let cigar = "100=";
        let target_len = parse_cigar(cigar).unwrap().iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>();
        let spacings = vec![10, 25, 50, 75, 100];
        
        for spacing in spacings {
            let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, spacing, false).unwrap();
            assert_eq!(alignment.trace_points.last().unwrap().target_pos, 100);
            
            let expected_points = (100 + spacing - 1) / spacing;
            assert_eq!(alignment.trace_points.len(), expected_points as usize);
        }
    }

    #[test]
    fn test_indel_at_boundary() {
        let cigar = "50=2I48=";
        let target_len = parse_cigar(cigar).unwrap().iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 50, false).unwrap();
        assert_eq!(alignment.trace_points.len(), 2);
        assert_eq!(alignment.trace_points[0].target_pos, 50);
        assert_eq!(alignment.trace_points[0].differences, 0);
        assert_eq!(alignment.trace_points[1].target_pos, 98);
        assert_eq!(alignment.trace_points[1].differences, 1);
    }

    #[test]
    fn test_mixed_operations() {
        let cigar = "25=10X15=5I20=5D25=";
        let target_len = parse_cigar(cigar).unwrap().iter()
        .map(|op| match op.op() {
            '=' | 'X' | 'D' | 'M' => op.len() as i64,
            'I' => 0,
            _ => panic!("Invalid op"),
        })
        .sum::<i64>();
        let alignment = cigar_to_trace_points(cigar, 0, 0, target_len, 30, false).unwrap();
        
        // Verify total length and differences are tracked correctly
        let total_target_len: i32 = parse_cigar(cigar)
            .unwrap()
            .iter()
            .map(|op| match op.op() {
                '=' | 'X' | 'D' | 'M' => op.len(),
                'I' => 0,
                _ => panic!("Invalid op"),
            })
            .sum();
            
        assert_eq!(
            alignment.target_end - alignment.target_start,
            total_target_len as i64
        );
    }
}