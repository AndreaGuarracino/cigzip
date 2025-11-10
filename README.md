# cigzip

Convert CIGAR strings to compact tracepoints and back, with efficient binary storage.

## Quick Start

```sh
# Encode CIGAR to tracepoints
cigzip encode --paf input.paf > output.tp.paf

# Decode tracepoints back to CIGAR
cigzip decode --paf output.tp.paf --sequence-files ref.fa > restored.paf

# Binary compression (6x compression, ~1min for 4GB)
cigzip compress -i input.paf -o output.bpaf

# Decompress (auto-detects format)
cigzip decompress -i output.bpaf -o output.paf
```

## Commands

### encode - CIGAR to tracepoints

```sh
# Text output (default)
cigzip encode --paf input.paf > output.tp.paf

# Binary output
cigzip encode --paf input.paf --type fastga \
  --output-format binary --output-file output.bpaf
```

**Options:**
- `--type`: `standard` (default), `mixed`, `variable`, `fastga`
- `--complexity-metric`: `edit-distance` (default), `diagonal-distance`
- `--max-complexity N`: Segment size limit (default: 32, fastga: 100)
- `--output-format`: `text` (stdout), `binary` (requires `--output-file`)
- `--threads N`: Parallel threads (default: 4)

### decode - Tracepoints to CIGAR

```sh
# Auto-detects text or binary format
cigzip decode --paf input.paf --sequence-files ref.fa > output.paf

# FASTA list file
cigzip decode --paf input.paf --sequence-list refs.txt > output.paf

# Custom gap penalties
cigzip decode --paf input.paf --sequence-files ref.fa \
  --distance gap-affine-2p --penalties 5,8,2,24,1 > output.paf

# Banded alignment (faster)
cigzip decode --paf input.paf --sequence-files ref.fa \
  --heuristic --max-complexity 100 > output.paf
```

**Options:**
- `--sequence-files`: FASTA file(s) with reference sequences
- `--sequence-list`: File containing FASTA paths (one per line)
- `--distance`: `edit` (default), `gap-affine`, `gap-affine-2p`
- `--penalties`: Gap penalty values (default: 5,8,2,24,1)
- `--heuristic`: Enable banded alignment
- `--keep-old-stats`: Preserve original gi/bi/sc as giold/biold/scold

### compress - PAF with tracepoints to binary

```sh
# Compress with automatic strategy
cigzip compress -i input.paf -o output.bpaf

# Compress with specific strategy (automatic, varint-zstd, delta-varint-zstd)
cigzip compress -i input.paf -o output.bpaf --strategy varint-zstd

# Set compression level (1-22, default 3)
cigzip compress -i input.paf -o output.bpaf --strategy automatic,9

# From stdin
cat input.paf | cigzip compress -i - -o output.bpaf
```

**Options:**
- `-i, --input`: Input PAF file with `tp:Z:` tags (or `-` for stdin)
- `-o, --output`: Output binary file
- `--strategy`: `varint` (default, recommended) or `varint-raw`

### decompress - Binary to PAF with tracepoints

```sh
# To file
cigzip decompress -i input.bpaf -o output.paf

# To stdout
cigzip decompress -i input.bpaf -o -
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
        TracepointData::Standard(tps) => {
            println!("Standard tracepoints: {} items", tps.len());
        }
        TracepointData::Fastga(tps) => {
            println!("FastGA tracepoints: {} items", tps.len());
        }
        _ => {}
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

See examples:
- `examples/seek_demo.rs` - O(1) random access
- `examples/offset_demo.rs` - Offset-based access

## Building

```sh
git clone https://github.com/AndreaGuarracino/cigzip
cd cigzip
cargo build --release
```

Binary: `target/release/cigzip`

## Format Details

### Binary Format (BPAF)

```
[Header] → [Records] → [StringTable]
```

- **Fast O(1) random access**: Byte-aligned varint encoding enables instant tracepoint extraction
- **Compression**: Automatic encoding selection + varint + zstd (configurable level 1-22)
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
