# Alignment-tracepoints

A Rust implementation for alignment reconstruction that uses variable‐sized tracepoints to efficiently store and rebuild full CIGAR strings.

## Overview
This project implements an efficient approach to sequence alignment storage using variable-sized tracepoints. Instead of storing complete CIGAR strings or using fixed intervals, alignments are recorded as tracepoints that adapt their size based on the local alignment complexity. This approach achieves:

Key benefits:
 - **Space Efficiency:** Significant reduction in alignment storage size.
 - **Reconstruction Speed:** Fast recovery of full alignment details.
 - **Flexibility:** Adjustable tracepoint intervals (delta) to balance storage vs. runtime.

This design is inspired by Gene Myers' original tracepoint concept and extended here with multiple alignment algorithms.

## Features

- **Variable-sized Tracepoints:** Convert CIGAR strings into variable‐sized tracepoints that adapt based on a configurable diff threshold.
- **CIGAR Operations:** Parse and generate extended format CIGAR strings (including =, X, I, D, M).
- **Efficient Storage:** Compact representation that adapts to local alignment complexity.
- **Accurate Reconstruction:** Rebuild full CIGAR strings from tracepoint data using WFA2-lib.

## How It Works
This implementation uses variable-sized tracepoints that accumulate bases and differences until reaching a configurable diff threshold. Longer indels trigger special handling, ensuring that the reconstructed CIGAR string is accurate while minimizing storage. Each tracepoint stores:

- The number of bases consumed in sequence A
- The number of bases consumed in sequence B  
- The number of differences in that segment

## Quick Start

To get started, clone the repository and build the project with Cargo:

```bash
git clone trace_points
cd trace_points
cargo build --release
```

Run the example driver to see variable-sized tracepoints in action for efficient alignment reconstruction:

```bash
cargo run --release -- --max-diff 128
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
