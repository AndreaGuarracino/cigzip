

# cigzip

Herramienta de línea de comandos para convertir cadenas CIGAR a tracepoints compactos y viceversa, con almacenamiento binario eficiente.

## Descripción general

`cigzip` se construye sobre la biblioteca [tracepoints](https://github.com/AndreaGuarracino/tracepoints) (muestreo de tracepoints y reconstrucción de alineamientos) y el formato binario [tpa](https://github.com/AndreaGuarracino/tpa).

## Instalación

```sh
git clone https://github.com/AndreaGuarracino/cigzip
cd cigzip
cargo build --release   # binario: target/release/cigzip
```

Compila sin soporte para archivos AGC con `cargo build --release --no-default-features`.

## Uso

El repositorio proporciona un alineamiento de ejemplo en [`examples/`](examples/):

```sh
cigzip encode     --paf examples/example.paf -o example.tp.paf                                  # CIGAR -> tracepoints
cigzip compress   -i example.tp.paf -o example.tpa                                              # texto -> binario TPA (+ example.tpa.idx)
cigzip decompress -i example.tpa    -o example.tp.back.paf                                      # binario TPA -> texto
cigzip decode     --paf example.tp.back.paf --sequence-files examples/example.fa --max-complexity 32 > restored.paf
```

`decode` lee las secuencias desde archivos FASTA o AGC (`--sequence-files`). Ejecuta `cigzip <comando> --help` para ver todas las opciones de `encode`, `decode`, `compress` y `decompress`.

### Tracepoints FastGA

Usa `--type fastga` para tracepoints con espaciado fijo compatibles con [FastGA](https://github.com/thegenemyers/FASTGA) (el espaciado se establece con `--max-complexity`, por defecto 100). En ensamblajes con gaps, proporciona una tabla de contigs (`--fastga-contigs`: `scaffold<TAB>sbeg<TAB>send` por corrida de ACGT) para que la codificación se divida en las brechas N y se mantenga idéntica a `PAFtoALN` de FastGA; `decode` requiere la misma tabla.

```sh
cigzip compress -i aln.cg.paf -o aln.tpa --type fastga --fastga-contigs contigs.tsv                 # CIGAR PAF -> TPA, en un paso
cigzip decode   --paf aln.tp.paf --type fastga --fastga-contigs contigs.tsv --sequence-files genome.fa > restored.paf
```

## Repositorios relacionados

- **[tracepoints](https://github.com/AndreaGuarracino/tracepoints)**: la biblioteca principal que implementa el muestreo y la reconstrucción de tracepoints, sobre la cual se construye cigzip.
- **[tpa](https://github.com/AndreaGuarracino/tpa)**: la biblioteca del formato binario TracePoint Alignment (TPA) utilizada para almacenamiento y acceso aleatorio.

## Licencia

MIT
