# Alignment-tracepoints

A Rust implementation of sequence alignment with tracepoints, based on Gene Myers' concept for efficient alignment storage and reconstruction. This implementation provides functionality for converting between CIGAR strings and tracepoint representations of sequence alignments.

## Overview
Tracepoints offer a space-efficient way to store sequence alignments while allowing quick reconstruction. Instead of storing full CIGAR strings or edit operations, alignments are stored as pairs of (differences, bases) at regular intervals. This provides a configurable trade-off between storage space and reconstruction time.

## Features

Parse and generate CIGAR strings in extended format (=,X,I,D,M operations)
Convert CIGAR strings to tracepoint representation
Reconstruct CIGAR strings from tracepoints
Multiple alignment algorithms:
Basic Needleman-Wunsch with unit costs
- Affine gap penalties
- Dual affine gap penalties (separate costs for insertions and deletions)

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

## References
Based on Gene Myers' tracepoint concept described in:
[Recording Alignments with Trace Points](https://dazzlerblog.wordpress.com/2015/11/05/trace-points/).
