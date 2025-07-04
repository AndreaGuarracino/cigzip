# cigzip

A tool for compression and decompression of alignment CIGAR strings using tracepoints.

## Overview

`cigzip` leverages an efficient approach to sequence alignment storage using tracepoints. Instead of storing complete CIGAR strings, alignments are compressed into tracepoints that adapt to the local alignment complexity, reducing storage size while allowing for full reconstruction.

## Features

- **Compression**: Convert CIGAR strings to tracepoints for compact storage
- **Decompression**: Reconstruct full CIGAR strings from tracepoints
- **Parallel Processing**: Multi-threaded operation with chunk-based processing
- **Configurable Parameters**: Adjust max-diff, gap penalties, and other settings
- **Support for Compressed Files**: Automatically handles gzipped and bgzipped PAF files

## Usage

```shell
# Compress alignments in a PAF file (convert CIGAR to tracepoints)
cigzip compress --paf alignments.paf [--max-diff 32] [--threads 4] > alignments.tp.paf

# Decompress alignments (convert tracepoints back to CIGAR)
cigzip decompress --paf alignments.tp.paf --query-fasta query.fa --target-fasta target.fa [--penalties "5,8,2,24,1"] [--threads 4] > alignments.cigar.paf
```

### Command Options

#### Common Options
- `--paf FILE`: PAF file containing alignments (use "-" for stdin)
- `--threads N`: Number of threads to use (default: 4)
- `--verbose N`: Verbosity level (0=error, 1=info, 2=debug)

#### Compress-specific Options
- `--max-diff N`: Max-diff value for tracepoints (default: 32)

#### Decompress-specific Options
- `--query-fasta FILE`: FASTA file containing query sequences
- `--target-fasta FILE`: FASTA file containing target sequences
- `--penalties STRING`: Gap penalties in format "mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2" (default: "5,8,2,24,1")

## Building

```shell
git clone https://github.com/AndreaGuarracino/cigzip
cd cigzip
cargo build --release
```

### For GUIX's slaves

```bash
git clone --recursive https://github.com/AndreaGuarracino/cigzip
cd cigzip/WFA2-lib
env -i bash -c 'PATH=/usr/local/bin:/usr/bin:/bin ~/.cargo/bin/cargo build --release'
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References
Inspired by Gene Myers' tracepoint concept described in [Recording Alignments with Trace Points](https://dazzlerblog.wordpress.com/2015/11/05/trace-points/).