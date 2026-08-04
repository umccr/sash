---
name: Release validation
about: Pre-merge validation checklist for a sash release or PR stack
title: "Release validation — sash <version>"
labels: validation
---

See [`docs/testing-validation.md`](../../docs/testing-validation.md) for the full procedure.

## Scope

PRs under test: #\_\_\_

## Test samples

| Sample                               | Mode      | S3 input                                                                                                                                         | Purpose                      |
| ------------------------------------ | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------- |
| SEQC-II SBJ00480 (L2300943/L2300950) | DRAGEN+OA | `s3://test-data-503977275616-ap-southeast-2/testdata/analysis/production/{dragen-wgts-dna/4.4.6,oncoanalyser-wgts-dna/2.2.0}/SEQC-II/`           | regression vs 0.6.4 baseline |
| L2600284/L2600285                    | OA-only   | `s3://pipeline-prod-cache-503977275616-ap-southeast-2/byob-icav2/production/analysis/oncoanalyser-wgts-dna/20260716eec549b2/L2600284__L2600285/` | OA-only e2e (real v2.2.0)    |
| synthetic 550k                       | PCGR-skip | see `docs/testing-validation.md`                                                                                                                 | >450k skip smoke             |

Samplesheets: [`assets/test/`](../../assets/test).

## Procedure

- Standard / OA-only run: `docs/testing-validation.md`
- Release comparison: Sash Regression Service (`umccr/service-sash-regression`, case SEQC-II-medium)
- NATA: `umccr/accreditation/wgs-variant-monitoring-workflow/docs/`

## Results

| PR  | Sample | Run ID | Result | Concordance / notes |
| --- | ------ | ------ | ------ | ------------------- |
|     |        |        |        |                     |

## Sign-off

- [ ] somatic PASS concordance ≥ 81.2% vs baseline
- [ ] PURPLE purity/ploidy identical to baseline
- [ ] all cancer reports render
