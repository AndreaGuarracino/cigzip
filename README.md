# cigzip

Command-line tool to convert CIGAR strings to compact tracepoints and back, with efficient binary storage.

## Overview

`cigzip` is built on the [tracepoints](https://github.com/AndreaGuarracino/tracepoints) library (tracepoint sampling and alignment reconstruction) and the [tpa](https://github.com/AndreaGuarracino/tpa) binary format.

## Installation

```sh
git clone https://github.com/AndreaGuarracino/cigzip
cd cigzip
cargo build --release   # binary: target/release/cigzip
```

Build without AGC archive support with `cargo build --release --no-default-features`.

## Usage

The repository provides a small alignment in [`examples/`](examples/):

```sh
cigzip encode     --paf examples/example.paf --distance edit --max-complexity 32 -o example.tp.paf  # CIGAR -> tracepoints
cigzip compress   -i example.tp.paf -o example.tpa                                                   # text -> binary TPA (+ example.tpa.idx)
cigzip decompress -i example.tpa    -o example.tp.back.paf                                           # binary TPA -> text
cigzip decode     --paf example.tp.back.paf --sequence-files examples/example.fa \
                  --distance edit --max-complexity 32 > restored.paf

# restored.paf carries the same CIGAR as examples/example.paf
```

`decode` reads the sequences from FASTA or AGC (`--sequence-files`). Run `cigzip <command> --help` for the full options of `encode`, `decode`, `compress`, and `decompress`.

## Related repositories

- **[tracepoints](https://github.com/AndreaGuarracino/tracepoints)**: the core library implementing tracepoint sampling and reconstruction, which cigzip is built on.
- **[tpa](https://github.com/AndreaGuarracino/tpa)**: the TracePoint Alignment (TPA) binary format library used for storage and random access.

## License

MIT
