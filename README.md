# cigzip

Convert CIGAR strings to compact tracepoints and back, with efficient binary storage.

## Quick Start

```sh
# Encode CIGAR to tracepoints
cigzip encode --paf input.paf > output.tp.paf

# Decode tracepoints back to CIGAR (supports FASTA and AGC)
cigzip decode --paf output.tp.paf --sequence-files ref.fa > restored.paf
cigzip decode --paf output.tp.paf --sequence-files archive.agc > restored.paf

# Compress to binary format
cigzip compress -i input.tp.paf -o output.bpaf

# Decompress
cigzip decompress -i output.bpaf -o output.paf
```

## Commands

### encode - CIGAR to tracepoints

```sh
# Text output (default)
cigzip encode --paf input.paf > output.tp.paf

# Binary output
cigzip encode --paf input.paf --type standard \
  --output-format binary --output-file output.bpaf --strategy varint
```

**Options:**
- `--paf FILE`: Input PAF file with CIGAR strings
- `--type TYPE`: Tracepoint type - `standard` (default), `mixed`, `variable`, `fastga`
- `--complexity-metric METRIC`: `edit-distance` (default), `diagonal-distance`
- `--max-complexity N`: Segment size limit (default: 32, fastga: 100)
- `--output-format FORMAT`: `text` (stdout, default), `binary` (requires `--output-file`)
- `--output-file FILE`: Output file path (required for binary format)
- `--strategy STRATEGY`: Compression strategy for binary - `varint` (default), `huffman`
- `--distance DIST`: Distance metric for alignment - `edit` (default), `gap-affine`, `gap-affine-2p`
- `--penalties VALUES`: Gap penalty values (format depends on distance)
- `--heuristic`: Enable banded alignment (edit distance only, requires `--max-complexity`)
- `--threads N`: Parallel threads (default: 4)
- `--verbose N`: Verbosity level - 0 (error), 1 (info, default), 2 (debug)

### decode - Tracepoints to CIGAR

```sh
# Auto-detects text or binary format
cigzip decode --paf input.tp.paf --sequence-files ref.fa > output.paf

# AGC format (Assembled Genomes Compressor)
cigzip decode --paf input.tp.paf --sequence-files genomes.agc > output.paf

# FASTA list file
cigzip decode --paf input.tp.paf --sequence-list refs.txt > output.paf

# Custom gap penalties
cigzip decode --paf input.tp.paf --sequence-files ref.fa \
  --distance gap-affine-2p --penalties 5,8,2,24,1 > output.paf

# Banded alignment (faster)
cigzip decode --paf input.tp.paf --sequence-files ref.fa \
  --heuristic --max-complexity 100 > output.paf
```

**Sequence File Options:**
- `--sequence-files FILE...`: Reference sequences (FASTA or AGC format)
  - FASTA: `.fa`, `.fasta`, `.fna` (optionally gzipped)
  - AGC: `.agc` archives (query as `contig` or `contig@sample`)
  - Cannot mix FASTA and AGC files
- `--sequence-list FILE`: File containing sequence paths (one per line)

**Alignment Options:**
- `--type TYPE`: Tracepoint type (must match encoding) - `standard` (default), `mixed`, `variable`, `fastga`
- `--complexity-metric METRIC`: `edit-distance` (default), `diagonal-distance`
- `--distance DIST`: Distance metric - `edit` (default), `gap-affine`, `gap-affine-2p`
- `--penalties VALUES`: Gap penalty values (default: 5,8,2,24,1)
- `--heuristic`: Enable banded alignment (edit distance only, requires `--max-complexity`)
- `--max-complexity N`: Maximum complexity for heuristic banding
- `--trace-spacing N`: Trace spacing for fastga type (default: 100)

**Other Options:**
- `--keep-old-stats`: Preserve original gi/bi/sc as giold/biold/scold
- `--threads N`: Parallel threads (default: 4)
- `--verbose N`: Verbosity level - 0 (error), 1 (info, default), 2 (debug)

### compress - Binary PAF compression

Compress PAF files with tracepoints to binary format.

```sh
# Basic compression
cigzip compress -i input.tp.paf -o output.bpaf

# With specific strategy
cigzip compress -i input.tp.paf -o output.bpaf --strategy huffman

# From stdin
cat input.tp.paf | cigzip compress -i - -o output.bpaf
```

**Options:**
- `-i, --input FILE`: Input PAF file with `tp:Z:` tags (or `-` for stdin)
- `-o, --output FILE`: Output binary file
- `--strategy STRATEGY`: Compression strategy - `varint` (default), `huffman`
- `--verbose N`: Verbosity level - 0 (error), 1 (info, default), 2 (debug)

### decompress - Binary to PAF

Decompress binary PAF files to text format with tracepoints.

```sh
# To file
cigzip decompress -i input.bpaf -o output.tp.paf

# To stdout
cigzip decompress -i input.bpaf -o -
```

**Options:**
- `-i, --input FILE`: Input binary PAF file
- `-o, --output FILE`: Output file (or `-` for stdout)
- `--verbose N`: Verbosity level - 0 (error), 1 (info, default), 2 (debug)

## AGC Format Support

cigzip supports AGC (Assembled Genomes Compressor) archives for sequence files. AGC provides efficient compression for pangenome sequences with multiple samples.

**Query formats:**
- `contig` - Query by contig name (must be unique across samples)
- `contig@sample` - Explicitly specify sample name

**Requirements:**
- All sequence files must be same format (all FASTA or all AGC)
- AGC support enabled by default (see Building section to disable)

**Example:**
```sh
cigzip decode --paf alignments.paf --sequence-files genomes.agc > output.paf
```

## Library Usage

Use cigzip as a Rust library for programmatic access to BPAF files:

```rust
use lib_bpaf::BpafReader;

fn main() -> std::io::Result<()> {
    // Open BPAF file with index
    let mut reader = BpafReader::open("alignments.bpaf")?;
    println!("Total records: {}", reader.len());

    // O(1) random access by record ID
    let record = reader.get_alignment_record(1000)?;
    let (tracepoints, tp_type, _, _) = reader.get_tracepoints(1000)?;

    match &tracepoints {
        TracepointType::Standard(tps) => {
            println!("Standard tracepoints: {} items", tps.len());
        }
        TracepointType::Mixed(items) => {
            println!("Mixed tracepoints: {} items", items.len());
        }
        TracepointType::Variable(tps) => {
            println!("Variable tracepoints: {} items", tps.len());
        }
        TracepointType::Fastga(tps) => {
            println!("FastGA tracepoints: {} items", tps.len());
        }
    }

    // Sequential iteration
    for record in reader.iter_records() {
        let record = record?;
        // Process each alignment...
    }

    Ok(())
}
```

Add to `Cargo.toml`:
```toml
[dependencies]
cigzip = { git = "https://github.com/AndreaGuarracino/cigzip" }
lib_bpaf = { git = "https://github.com/AndreaGuarracino/lib_bpaf" }
```

## Building

### Standard Build (with AGC support)

```sh
git clone https://github.com/AndreaGuarracino/cigzip
cd cigzip
cargo build --release
```

Binary: `target/release/cigzip`

### Minimal Build (without AGC support)

To build without AGC support (e.g., in resource-constrained environments):

```sh
cargo build --release --no-default-features
```

This disables the optional `agc-rs` dependency.

## Format Details

### Binary Format (BPAF)

```
[Header] → [Records] → [StringTable]
```

**Features:**
- **Fast O(1) random access**: Byte-aligned varint encoding enables instant tracepoint extraction
- **Compression**: Delta encoding + varint + zstd level 3 (or Huffman + zstd)
- **Deduplication**: Shared string table for sequence names
- **Random access**: External `.bpaf.idx` index for O(1) record lookup
- **Backwards compatible**: Reads all format versions

### Index Format

```
Magic:     BPAI (4 bytes)
Version:   1 (1 byte)
Count:     varint (number of records)
Offsets:   varint[] (byte positions)
```

Index enables instant access to any alignment without scanning.

## License

MIT License - see [LICENSE](LICENSE) file

## References

Inspired by Gene Myers' tracepoint concept: [Recording Alignments with Trace Points](https://dazzlerblog.wordpress.com/2015/11/05/trace-points/)
