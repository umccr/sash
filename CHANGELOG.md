# umccr/sash: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of umccr/sash, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

## [0.6.0] – 2025-05-19

### Added

- _none_

### Changed

- SV caller switched GRIPSS → eSVe
- Cancer-report Structural Variants plot: `SR`→`SF`, `PR`→`DF`
- Linx upgraded 1.25 → 2.0
- Purple upgraded 4.0.1 → 4.1.0
- Bolt & GPGR updated for the adpat to above change

### Removed

- Kataegis module
- CHORD HRD metrics
- MSI load and status from purple

### Deprecated

- Metric aliases `SR`, `PR`

### Fixed

- _none_

### Security

- _none_

### Dependencies

| Tool | Old | New |
|------|-----|-----|
| Linx | 1.25 | 2.0 |
| Purple | 4.0.1 | 4.1.0 |
| Bolt | — | umccr/bolt#6 |
| GPGR | — | umccr/gpgr#88 |
