# umccr/sash: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of umccr/sash, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

## [0.7.0] - Unreleased

### Added

- Reference-data preparation step that maps and auto-extracts bundled resources (PCGR/VEP tarballs) from `--ref_data_path`. (see PR #38)
- Sigrap HRDetect and MutationalPatterns now run as standalone stages; CHORD predictions from oncoanalyser plus HRDetect/MutPat outputs are integrated into the cancer report. (see PR #36)
- Somatic VCFs are exported as MAF files via `vcf2maf` for downstream compatibility.
- Optional `--pcgr_variant_chunk_size` parameter to tune PCGR chunking in BOLT somatic annotation.

### Changed

- BOLT updated to `0.3.0-dev-20` (and GPGR report image to `0.3.0-dev-14`), including PCGR HTML renamed to `<sample>.pcgr.grch38.html` and requiring an explicit VEP cache path.
- PCGR reference bundle refreshed to `20250314` and resolved from tarballs alongside the GRCh38 VEP cache.
- DRAGEN input directories are handled as file objects to better support staged/symlinked inputs. (see PR #12)

### Dependencies

| Tool | Old | New |
|------|-----|-----|
| bolt | 0.2.17 | 0.3.0-dev-20 |
| pcgr data | 20220203 | 20250314 |
| sigrap | — | 0.2.0-dev-7 |

## [0.6.3] - 2025-10-07

### Added

- Validate input files in `prepare_run`: warn when optional files are missing; fail when mandatory files are missing. (see PR #29)

## [0.6.2] - 2025-10-1

### Changed

- Dragen HRD file facultatif

### Dependencies

| Tool | Old | New |
|------|-----|-----|
| bolt | 0.2.15 | 0.2.17 |
| gpgr | 2.2.0 | 2.2.1 |

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

| Tool | Old | New |
|------|-----|-----|
| eSVee | 1.0.3 | 1.1.2 |
| LINX | 2.0.2 | 2.1 |
| PAVE | 1.7.1 | 1.8 |
| PURPLE | 4.1 | 4.2 |
| bolt | 0.2.14 | 0.2.15 |
| gpgr | 2.1.3 | 2.2.0 |

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

| Tool | Old | New |
|------|-----|-----|
| Linx | 1.25 | 2.0 |
| Purple | 4.0.1 | 4.1.0 |
| Bolt | — | umccr/bolt#6 |
| GPGR | — | umccr/gpgr#88 |
