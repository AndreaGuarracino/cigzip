# Alignment-tracepoints

A Rust implementation of efficient sequence alignment storage and reconstruction using tracepoints, providing a significant reduction in storage requirements while maintaining fast alignment reconstruction capabilities.

## Overview
This project implements an efficient approach to sequence alignment storage using tracepoints. Instead of storing complete CIGAR strings, A-sequence (query) vs B-sequence (target) alignments are recorded as `(difference count, consumed bases in B-sequence)` pairs at regular intervals along the A-sequence. This approach achieves:

Key benefits:
 - **Space Efficiency:** Significant reduction in alignment storage size.
 - **Reconstruction Speed:** Fast recovery of full alignment details.
 - **Flexibility:** Adjustable tracepoint intervals (delta) to balance storage vs. runtime.

This design is inspired by Gene Myers' original tracepoint concept and extended here with multiple alignment algorithms.

## Features

- **CIGAR Operations:** Parse and generate extended format CIGAR strings (including =, X, I, D, M).
- **Tracepoint Conversion:** Convert CIGAR strings to a compact tracepoint representation.
- **Alignment Reconstruction:** Rebuild full CIGAR strings from tracepoint data.
- **Multiple Alignment Algorithms:**
   - Basic Needleman–Wunsch (unit-cost),
   - Affine gap penalties,
   - Dual affine gap penalties (with separate insertion and deletion costs).

## How It Works
Tracepoints record alignment information at regular intervals (controlled by the delta parameter). For each interval:

The number of differences (edits) in that segment
The number of bases consumed in sequence B

This creates a compact representation that can be used to efficiently reconstruct the full alignment when needed.

## Performance
Using tracepoints provides significant space savings compared to storing full CIGAR strings:

For a 10Kbp alignment with 15% error rate:

Full CIGAR: ~6KB
Tracepoints (delta=100): ~200 bytes

Alignment reconstruction time is O(εΔn) where:

ε = error rate
Δ = tracepoint interval
n = sequence length

## Quick Start

To get started, clone the repository and build the project with Cargo:

```bash
git clone trace_points
cd trace_points
cargo build --release
```

Run the example driver to see tracepoints and alignment reconstruction in action:

```bash
cargo run --release
```

For more details on the API and configuration (e.g., setting delta), refer to the source code documentation.

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests. For major changes, open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References
Based on Gene Myers' tracepoint concept described in:
[Recording Alignments with Trace Points](https://dazzlerblog.wordpress.com/2015/11/05/trace-points/).
