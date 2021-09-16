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
nextflow run thanhleviet/nf-campy-nanopore --reads "/path/to/folder/of/*.fastq.gz" --outdir "/path/to/output" -profile conda
```
