# cigzip

Encode and decode alignment CIGARs using tracepoints.

## Overview

`cigzip` converts CIGAR strings into compact tracepoints that adapt to local alignment complexity, and reconstructs the original CIGARs on demand.

## Features

- **Encode**: Convert CIGAR strings to tracepoints (compact storage)
- **Decode**: Reconstruct full CIGAR strings from tracepoints
- **Parallel Processing**: Multi-threaded operation with chunk-based processing
- **Configurable Parameters**: Choose tracepoint type/metric, penalties, and banded heuristics
- **Support for Compressed Files**: Automatically handles gzipped and bgzipped PAF files

## Usage

Encode (CIGAR → tracepoints)

```sh
# Standard, default settings
cigzip encode --paf in.paf > out.tp.paf

# Mixed type, diagonal metric, tighter segments
cigzip encode --paf in.paf --type mixed --complexity-metric diagonal-distance   --max-complexity 24 > out.mixed.tp.paf

# FastGA (spacing=100 default)
cigzip encode --paf in.paf --type fastga > out.fastga.tp.paf
```

Decode (tracepoints → CIGAR)

```sh
# Edit distance (unit costs)
cigzip decode --paf out.tp.paf --sequence-files ref.fa > out.cigar.paf

# Gap-affine 2p penalties
cigzip decode --paf out.tp.paf --sequence-files ref1.fa ref2.fa   --distance gap-affine-2p --penalties 5,8,2,24,1 > out.cigar.paf

# Banded heuristics (band from max-complexity)
cigzip decode --paf out.tp.paf --sequence-files ref.fa   --heuristics --max-complexity 100 > out.cigar.paf

# Many FASTAs via list file
cigzip decode --paf out.tp.paf --sequence-list refs.txt > out.cigar.paf

# Keep old stats while replacing gi/bi/sc
cigzip decode --paf out.tp.paf --sequence-files ref.fa --keep-old-stats > out.cigar.paf
```

### Command Options

#### Common Options
- `--paf FILE`: PAF file containing alignments (use "-" for stdin)
- `--threads N`: Number of threads to use (default: 4)
- `--verbose N`: Verbosity level (0=error, 1=info, 2=debug)

#### Encode
- `--type`: Tracepoint representation (`standard`, `mixed`, `variable`, `fastga`)
- `--complexity-metric`: `edit-distance` (default) or `diagonal-distance`
- `--max-complexity N`: Segmentation limit (default: 32; 100 for `fastga`)

#### Decode
- `--sequence-files`: One or more FASTA files with all referenced sequences
- `--sequence-list`: File containing additional FASTA paths (one per line)
- `--distance`: `edit` (unit costs), `gap-affine`, or `gap-affine-2p`
- `--penalties`: Comma list. For `gap-affine-2p`: `mismatch,gap_open1,gap_ext1,gap_open2,gap_ext2` (default `5,8,2,24,1`)
- `--heuristics --max-complexity N`: Enable banded static heuristics (works with both metrics). For `diagonal-distance`, the band equals `max-complexity`.
- `--keep-old-stats`: When replacing fields, also keep previous `gi/bi/sc` as `giold/biold/scold`

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
