# sash

**scwatts/sash** is the UMCCR post-processing WGS workflow.

## Table of contents

* [Requirements](#requirements)
* [Usage](#usage)

## Requirements

* Java
* Nextflow ≥22.10.6
* Docker

## Usage

```bash
nextflow run scwatts/sash \
    -profile docker \
    --input samplesheet.csv \
    --ref_data_path reference_data/ \
    --outdir output/
```
