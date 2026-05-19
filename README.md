# sash

`sash` is the UMCCR post-processing WGS workflow. The workflow takes DRAGEN small variant calls and oncoanalyser results
as input to perform annotation, prioritisation, rescue and filtering, and reporting for the WGS variant data.
Additionally, `sash` runs several sensors for biomarker assessment and genomic characterisation including HRD status,
mutational signatures, purity/ploidy, MSI, and TMB.

While the `sash` workflow utilises a range of tools and software, it is most closely coupled with
[bolt](https://github.com/scwatts/bolt), a Python package that implements the UMCCR post-processing logic and supporting
functionality.

## Table of contents

- [Summary](#summary)
- [Requirements](#requirements)
- [Usage](#usage)

## Summary

The general processes `sash` runs include:

- gpgr for generating the summary Cancer Report
- PCGR to report processed small somatic variants (annotated, rescued, filtered, prioritised)
- CPSR to report processed small germline variants (filtered, annotated, prioritised)
- linxreport to collate SV annotations and plots from LINX
- MultiQC for reporting various WGS statistics / metrics for QC
- SAGE variant calling to supplement DRAGEN small somatic variants
- PURPLE for TMB, MSI, CNV calling, and purity / ploidy estimation
- HRDetect for HRD inference; CHORD predictions from oncoanalyser are ingested when present
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

The `dragen_somatic_vcf` and `dragen_germline_vcf` filetypes are optional overrides that take precedence over
the paths constructed from `dragen_somatic_dir`/`dragen_germline_dir`. Use these when DRAGEN output filenames
do not match the `tumor_id`/`normal_id` used by oncoanalyser (e.g. cohorts with differing sample ID prefixes):

```text
id,subject_name,sample_name,filetype,filepath
subject_a.example,subject_a,sample_germline,dragen_germline_dir,/path/to/dragen_germline/
subject_a.example,subject_a,sample_somatic,dragen_somatic_dir,/path/to/dragen_somatic/
subject_a.example,subject_a,sample_somatic,oncoanalyser_dir,/path/to/oncoanalyser/
subject_a.example,subject_a,sample_somatic,dragen_somatic_vcf,/path/to/dragen_somatic/SAMPLE.hard-filtered.vcf.gz
subject_a.example,subject_a,sample_germline,dragen_germline_vcf,/path/to/dragen_germline/SAMPLE.hard-filtered.vcf.gz
```

The `pcgr_variant_chunk_size` parameter can be set to run PCGR in chunked mode for hypermutated samples.
When set, the somatic annotation process splits the VCF into chunks of the specified size before running PCGR,
preventing memory and runtime issues on highly mutated tumours:

```bash
nextflow run umccr/sash \
  --pcgr_variant_chunk_size 50000 \
  ...
```

Execute analysis

```bash
nextflow run umccr/sash \
  -profile docker \
  --input samplesheet.csv \
  --ref_data_path /path/to/reference_data/ \
  --outdir output/
```

For detailed instructions, see [docs/usage.md](docs/usage.md).

## Documentation

The `sash` pipeline comes with documentation in the `docs/` directory:

- [Usage](docs/usage.md): Detailed instructions on how to run the pipeline.
- [Output](docs/output.md): Description of the output files and reports.
- [Details](docs/details.md): In-depth explanation of the pipeline steps and tools.
- [ADR](docs/adr.md): Architectural Decision Records.

## Citations

You can cite a specific version of `sash` from the Zenodo record [10.5281/zenodo.15833492](https://doi.org/10.5281/zenodo.15833492) such as:

> Watts, S. C., Savelyev, V., Diakumis, P., Clayssen, Q., Mitchell, H., & Hofmann, O. (2025). umccr/sash: 0.7.0 (0.7.0). Zenodo. <https://doi.org/10.5281/zenodo.15833493>
