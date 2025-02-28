# cigzip

A tool for compression and decompression of alignment CIGAR strings using tracepoints.

## Overview

`cigzip` implements an efficient approach to sequence alignment storage using variable-sized tracepoints. Instead of storing complete CIGAR strings, alignments are compressed into tracepoints that adapt to the local alignment complexity, reducing storage size while allowing for full reconstruction.

## Features

- **Compression**: Convert CIGAR strings to tracepoints for compact storage
- **Decompression**: Reconstruct full CIGAR strings from tracepoints
- **Banded Tracepoints**: Tracepoint representation for faster decompression
- **Parallel Processing**: Multi-threaded operation
- **Configurable Parameters**: Adjust max-diff, gap penalties, and other settings

## Usage

```shell
# Compress alignments in a PAF file (convert CIGAR to tracepoints)
cigzip compress --paf alignments.paf [--banded] [--max-diff 128] [--threads 4] > alignments.tp.paf

# Decompress alignments (convert tracepoints back to CIGAR)
cigzip decompress --paf alignments.tp.paf --fasta sequences.fa [--penalties "3,4,2,24,1"] [--threads 4] > alignments.cigar.paf
```

## Building

You need to build `WFA2-lib` first, which is a submodule of this repository. To do so, run:

```shell
git clone --recursive https://github.com/AndreaGuarracino/cigzip
cd cigzip/WFA2-lib
make clean all
cd ..
```

Then, you can build the project using Cargo:

```shell
# Point to your pre-built WFA2-lib directory
export WFA2LIB_PATH="./WFA2-lib"

# Build your project
cargo build --release
```

### For GUIX's slaves

```bash
git clone --recursive https://github.com/AndreaGuarracino/cigzip
cd cigzip/WFA2-lib
guix shell -C -D -f guix.scm
export CC=gcc; make clean all
exit
cd ..
env -i bash -c 'WFA2LIB_PATH="./WFA2-lib" PATH=/usr/local/bin:/usr/bin:/bin ~/.cargo/bin/cargo build --release'
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## References
Inspired by Gene Myers' tracepoint concept described in [Recording Alignments with Trace Points](https://dazzlerblog.wordpress.com/2015/11/05/trace-points/).
