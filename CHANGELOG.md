# umccr/sash: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of umccr/sash, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

## [0.7.0] - 2026-05-13

### Added

- Allow explicit DRAGEN VCF paths in the samplesheet via optional `dragen_somatic_vcf` and `dragen_germline_vcf` filetypes, falling back to the path constructed from `tumor_id` / `normal_id` when not provided. Enables samples where DRAGEN and oncoanalyser outputs use different sample ID prefixes (e.g. APGI cohort).
- VCF2MAF process for somatic PASS VCF conversion to MAF format
- SIGRAP_HRDETECT and SIGRAP_MUTPAT processes for HRDetect and mutational pattern analysis
- CUSTOM_EXTRACTTARBALL helper: PCGR and VEP reference data now shipped as tarballs and extracted at runtime
- `pcgr_variant_chunk_size` parameter for chunked PCGR runs on hypermutated samples
- `process_medium_memory` resource label
- Zenodo DOI

### Changed

- Updated to oncoanalyser v2.x data structure
- PCGR output filename: `*.pcgr_acmg.grch38.html` → `*.pcgr.grch38.html`
- BOLT_OTHER_CANCER_REPORT now consumes SIGRAP hrdetect and mutpat outputs
- `BOLT_OTHER_CANCER_REPORT` and `BOLT_SMLV_SOMATIC_ANNOTATE` relabelled to `process_medium_memory`
- Reference data: PCGR 20220203 → 20250314; VEP 113; PCGR/VEP data shipped as tarballs

### Dependencies

| Tool      | Old      | New      |
| --------- | -------- | -------- |
| bolt      | 0.2.17   | 0.3.0    |
| PCGR      | —        | v2.2.5   |
| PCGR data | 20220203 | 20250314 |
| VEP       | —        | 113      |

## [0.6.4] - 2026-03-31

### Fixed

- Revert hypermutated flag condition in cancer report (reverts >500k sage variant threshold, [gpgr#101](https://github.com/umccr/gpgr/pull/101), [bolt#28](https://github.com/umccr/bolt/pull/28)) ([#44](https://github.com/umccr/sash/issues/44))
- Include MANE Plus Clinical transcript ENST00000579755 for CDKN2A alongside MANE Select ENST00000304494, fixing missing gene entries in the Gene Somatic CNV Calls Table when only the NM_058195.4-specific exon has a copy number loss (`umccr_reference_data = '2--1'`) ([#45](https://github.com/umccr/sash/issues/45))

### Documentation

- Fix broken ToC anchor, absolute local path, duplicate word, and misplaced LINX section in `docs/details.md`
- Fix `scwatts/sash` → `umccr/sash` and stale CHORD reference in `README.md`
- Add `sigrap/` and `vcf2maf/` to output directory tree; fix inaccurate germline preparation description in `docs/output.md`

### Dependencies

| Tool | Old    | New    |
| ---- | ------ | ------ |
| bolt | 0.2.17 | 0.2.18 |

## [0.6.3] - 2025-10-07

### Dependencies

| Tool | Old    | New    |
| ---- | ------ | ------ |
| bolt | 0.2.15 | 0.2.17 |
| gpgr | 2.2.0  | 2.2.1  |

## [0.6.1] - 2025-09-16

### Added

- Support for oncoanalyser v2.2.0 data structure
- PAVE MNV filtering to remove variants with MNVTAG annotations (see [discussion](https://github.com/umccr/sash/issues/19))

### Changed

- **Reference data updates:**
  - Updated to hmf_pipeline_resources.38_v2.2.0--3
- **Process update:**
  - eSVee parameter updates
  - PAVE parameter updates for v1.8 compatibility
- Updated prepare_input paths for oncoanalyser v2.2.0 structure

### Removed

- Unused reference data files (gridss_region_blocklist, lilac, orange, sigs_signatures, disease_ontology)

### Dependencies

| Tool   | Old    | New    |
| ------ | ------ | ------ |
| eSVee  | 1.0.3  | 1.1.2  |
| LINX   | 2.0.2  | 2.1    |
| PAVE   | 1.7.1  | 1.8    |
| PURPLE | 4.1    | 4.2    |
| bolt   | 0.2.14 | 0.2.15 |
| gpgr   | 2.1.3  | 2.2.0  |

## [0.6.0] - 2025-06-04

### Added

- Software versions are now in `pipeline_info/software_versions.yml`

### Changed

- **SV caller:** GRIPSS → eSVee
  - SV counts (unmelted & melted)
  - CNV counts
  - TMB-SV counts
  - CopyNumberSegment counts
  - SV in Circos plots
  - Breakpoints & Breakends tables
  - Copy-number variants tables
  - Genome-wide somatic CNV segment tracks
  - SV Map visualisation
- Cancer-report Structural Variants summary plot:
  - SR (Split Read) → SF (Split Fragments)
  - PR (Paired-Read) → DF (Discordant Fragments)
- Linx v1.25 → v2.0 (affects all Linx reports/files)
- Purple v4.0.1 → v4.1.0
  - MSI calculation relay on SAGE-specific tags #7
  - (reverted) Circos have link sizes dependent on the size of SV #6
- Filter PoN SV in cancer report tables #8

### Removed

- Kataegis module
- CHORD HRD metrics
- MSI load and status from purple

### Deprecated

- Metric aliases `SR`, `PR`
- GRIDSS/GRIPSS modules

### Fixed

- Nextflow `Stub` run

### Security

- _none_

### Dependencies

| Tool   | Old   | New           |
| ------ | ----- | ------------- |
| Linx   | 1.25  | 2.0           |
| Purple | 4.0.1 | 4.1.0         |
| Bolt   | —     | umccr/bolt#6  |
| GPGR   | —     | umccr/gpgr#88 |
