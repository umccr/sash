# umccr/sash: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of umccr/sash, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

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
