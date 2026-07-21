# Testing and validation

How to validate a sash change against real cloud data before merge. Covers the standard
DRAGEN+OA path, the OA-only path, the hypermutated edge case, and release regression.

Ready-to-run samplesheets live in [`assets/test/`](../assets/test). All input S3 paths were
verified readable on 2026-07-22 via `umccr-prod-operator` / `umccr-development-ro`.

## Test samples

| Sample | Mode | Samplesheet | Purpose |
|--------|------|-------------|---------|
| SEQC-II — SBJ00480 (L2300943/L2300950, HCC1395) | DRAGEN+OA | `assets/test/samplesheet_dragen_oa_seqcii.csv` | full-pipeline regression |
| SEQC-II — SBJ00480 | OA-only | `assets/test/samplesheet_oa_only_seqcii.csv` | OA-only clean reference |
| L2600284/L2600285 (real prod, OA v2.2.0) | OA-only | `assets/test/samplesheet_oa_only_l2600284.csv` | OA-only real-data e2e |
| Synthetic 550k-variant VCF | PCGR-skip | generated (see below) | >450k skip smoke test |

SEQC-II is the primary regression sample. Its full package is staged in the read-only
testdata bucket:

```
s3://test-data-503977275616-ap-southeast-2/testdata/analysis/production/
  dragen-wgts-dna/4.4.6/SEQC-II/20260330a606412c/      # DRAGEN somatic + germline
  oncoanalyser-wgts-dna/2.2.0/SEQC-II/20260401d6eff920/ # OA output
  sash/0.6.4/SEQC-II/20260428471bb6ed/                  # baseline for concordance
```

## Run a validation

1. Build a disposable bolt image off the bolt PR branch (`v<version>-<repo><issue>`, e.g.
   `v0.3.2-sash59`). Wait for all image variants before continuing.
2. Check out the sash branch and point its `container` refs in `modules/local/bolt/` at the
   eval tag.
3. Run with `nextflow-stack`'s `run.sh` and `nextflow_aws.template.config`, passing one of the
   samplesheets above as `--input`. Pin Nextflow to the version in the workflow config and use
   AWS CLI v2.

OA-only mode reads `tumor_id` and `normal_id` from the samplesheet directly and accepts only
`oncoanalyser_dir` rows — `run.sh` still expects DRAGEN args, so pass a hand-written samplesheet
to `nextflow run` rather than using its samplesheet generation.

## Choosing an OA-only sample

The OA-only path expects the oncoanalyser v2.2.0 layout: `esvee/`,
`sage_calling/somatic/{tumor_id}.sage.somatic.vcf.gz`, and `pave/{tumor_id}.pave.germline.vcf.gz`
(tumor-id prefix). Older OA outputs use `gridss`/`gripss`, `sage/`, and germline PAVE named
`*.sage.germline.pave.vcf.gz` — `prepare_input.nf` cannot detect these. Confirm `esvee/` and
`pave/{tumor}.pave.germline.vcf.gz` exist before picking a sample, or re-run OA at v2.2.0.

## Hypermutated edge case (>450k variants)

The PCGR-skip / VCF2MAF-skip logic (sash #59, tracked by #52) needs a genuinely hypermutated
sample.

- **Smoke test:** a synthetic 550k-variant VCF exercises the skip branch directly and confirms
  the sample exits 0 without a full run. Generate one with a script that emits >450k PASS
  records, then feed it to bolt's PCGR prep.
- **Real-data sign-off:** the known real hypermutators (L2100242, 595k PASS; SBJ02862 /
  L2201449, the OOM case) roll off the prod cache over time. Re-stage one from the portal
  (`api/v1/s3/attributes?portalRunId=…`); a DRAGEN+OA re-run may be needed before sash.

## Release regression

Compare a new release against the previous one with the Sash Regression Service
([`umccr/service-sash-regression`](https://github.com/umccr/service-sash-regression)) on case
SEQC-II-medium (L2301218/L2301217). The 0.6.4 baseline is already wired. This replaces manual
output diffing.

## Variant monitoring (NATA)

Germline variant-drift monitoring on GIAB control libraries has its own SOPs in
[`umccr/accreditation`](https://github.com/umccr/accreditation/tree/main/wgs-variant-monitoring-workflow/docs)
(`RUNBOOK.md`, `SOP_OUTLIER_REVIEW.md`, `SOP_OUT_OF_LIMITS_REVIEW.md`).

## Sign-off criteria

- somatic PASS concordance ≥ 81.2% vs the 0.6.4 baseline
- PURPLE purity and ploidy identical to baseline
- all three cancer reports render without error
