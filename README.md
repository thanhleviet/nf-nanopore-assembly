Simple nextflow script to run de-novo assembly on nanopore data with Flye and polish the assembly with medaka.

## Requirements

- conda and mamba or docker for managing bioinformatics software
- This pipeline requires nextflow edge version to use mamba.
  Installing the edge version:
  ```bash
  export NXF_EDGE=1
  nextflow self-update
  ```
## Usage

```
nextflow run thanhleviet/nf-nanopore-assembly \
--reads "/path/to/folder/of/*.fastq.gz" \
--outdir "/path/to/output" \
--genome_size "1.7m" \
-profile conda
```
