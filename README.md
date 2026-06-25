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
cigzip encode     --paf examples/example.paf -o example.tp.paf                                  # CIGAR -> tracepoints
cigzip compress   -i example.tp.paf -o example.tpa                                              # text -> binary TPA (+ example.tpa.idx)
cigzip decompress -i example.tpa    -o example.tp.back.paf                                      # binary TPA -> text
cigzip decode     --paf example.tp.back.paf --sequence-files examples/example.fa --max-complexity 32 > restored.paf
```

`decode` reads the sequences from FASTA or AGC (`--sequence-files`). Run `cigzip <command> --help` for the full options of `encode`, `decode`, `compress`, and `decompress`.

### FastGA tracepoints

Use `--type fastga` for fixed-spacing tracepoints compatible with [FastGA](https://github.com/thegenemyers/FASTGA) (spacing set by `--max-complexity`, default 100). On gapped assemblies, pass a contig table (`--fastga-contigs`: `scaffold<TAB>sbeg<TAB>send` per ACGT run) so encoding splits at N-gaps and stays identical to FastGA's `PAFtoALN`; `decode` needs the same table.

```sh
cigzip compress -i aln.cg.paf -o aln.tpa --type fastga --fastga-contigs contigs.tsv                 # CIGAR PAF -> TPA, one step
cigzip decode   --paf aln.tp.paf --type fastga --fastga-contigs contigs.tsv --sequence-files genome.fa > restored.paf
```

## Related repositories

- **[tracepoints](https://github.com/AndreaGuarracino/tracepoints)**: the core library implementing tracepoint sampling and reconstruction, which cigzip is built on.
- **[tpa](https://github.com/AndreaGuarracino/tpa)**: the TracePoint Alignment (TPA) binary format library used for storage and random access.

## License

MIT
