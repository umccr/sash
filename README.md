# sash

`sash` is the UMCCR post-processing WGS workflow. The workflow takes DRAGEN small variant calls and oncoanalyser results
as input to perform annotation, prioritisation, rescue and filtering, and reporting for the WGS variant data.
Additionally, `sash` runs several sensors for biomarker assessment and genomic characterisation including HRD status,
mutational signatures, purity/ploidy, MSI, and TMB.

While the `sash` workflow utilises a range of tools and software, it is most closely coupled with
[bolt](https://github.com/scwatts/bolt), a Python package that implements the UMCCR post-processing logic and supporting
functionality.

## Table of contents

* [Summary](#summary)
* [Requirements](#requirements)
* [Usage](#usage)

## Summary

The general processes `sash` runs include:

- gpgr for generating the summary Cancer Report
- PCGR to report processed small somatic variants (annotated, rescued, filtered, prioritised)
- CPSR to report processed small germline variants (filtered, annotated, prioritised)
- linxreport to collate SV annotations and plots from LINX
- MultiQC for reporting various WGS statistics / metrics for QC
- SAGE variant calling to supplement DRAGEN small somatic variants
- PURPLE for TMB, MSI, CNV calling, and purity / ploidy estimation
- HRDetect and CHORD for HRD inference
- MutationalPatterns to fit mutational signatures
- PAVE for somatic variant annotation with MNV filtering (see [discussion](https://github.com/umccr/sash/issues/19))

## Requirements

- Java
- Nextflow ≥22.10.6
- Docker

## Usage

Create a samplesheet

```text
id,subject_name,sample_name,filetype,filepath
subject_a.example,subject_a,sample_germline,dragen_germline_dir,/path/to/dragen_germline/
subject_a.example,subject_a,sample_somatic,dragen_somatic_dir,/path/to/dragen_somatic/
subject_a.example,subject_a,sample_somatic,oncoanalyser_dir,/path/to/oncoanalyser/
```

Execute analysis

```bash
nextflow run scwatts/sash \
  -profile docker \
  --input samplesheet.csv \
  --ref_data_path /path/to/reference_data/ \
  --outdir output/
```
