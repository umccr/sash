# Testing and validation

Four tiers, cheapest/fastest first. Use tier 1 for every change, tier 2 for workflow-logic
changes, tier 3 before merging anything touching clinically-relevant behavior, tier 4 only for
release candidates going into production.

| Tier | What | When |
|------|------|------|
| 1 — [Local test](#tier-1--local-test) | `nextflow run` on your machine | every change |
| 2 — [nf-test](#tier-2--nf-test) | unit/integration tests | workflow-logic changes |
| 3 — [Nextflow stack on AWS](#tier-3--nextflow-stack-on-aws) | real cloud data, AWS Batch | pre-merge, clinically-relevant changes |
| 4 — [ICA](#tier-4--ica-production-gate) | production deployment gate | release candidates |

## Tier 1 — Local test

```
nextflow run . -profile test,docker --outdir <dir>
```

**Known gap:** `conf/test.config` is still unmodified nf-core template boilerplate (points at
`nf-core/test-datasets` viralrecon/SARS-CoV-2 placeholder data, not real sash fixtures) — it
doesn't currently produce a meaningful run. Building a proper minimal `-profile test` fixture set
under `assets/test/` is tracked as follow-up work, not done as part of this doc.

For a fast structural/syntax check that doesn't need real data (e.g. while iterating on a Nextflow
or strict-syntax change), use `-stub-run`:

```
nextflow run . -profile test,docker -stub-run --outdir <dir>
```

This compiles the full include graph (`main.nf` → `workflows/sash.nf` → every subworkflow and
module it references) without executing real logic — it's usually enough to surface parser/syntax
errors and missing-file/param errors without waiting on a full pipeline run.

## Tier 2 — nf-test

```
nf-test test
```

`nf-test.config` (repo root) points nf-test at the whole tree, so it discovers every `*.nf.test`
file automatically. Existing coverage:

- `subworkflows/local/tests/prepare_input/`
- `subworkflows/local/tests/prepare_reference/`
- `subworkflows/nf-core/utils_nfcore_pipeline/tests/` (vendored nf-core test — 2 of 4 function
  tests currently fail on a pre-existing missing-fixture path, unrelated to sash code)
- `modules/local/pave/somatic/tests/`
- `modules/local/bolt/smlv_somatic/report/tests/`

Add new tests alongside the subworkflow/module they cover, following the `nextflow_process` /
`nextflow_workflow` block pattern already used in `prepare_input`'s test. nf-test covers
logic/wiring correctness (e.g. "does `prepare_input` pick the right code path for a given
samplesheet row") — it is **not** a substitute for tiers 3/4, which validate scientific output
correctness on real data.

**Known issue — strict syntax vs nf-test:** nf-test 0.9.2's generated test-harness wrapper scripts
use the Groovy spread operator (`PREPARE_INPUT(*input)`) to splat test inputs into the call. This
is itself banned under strict syntax, so running nf-test with the strict parser active (the
Nextflow 26.04+ default) fails every test with `Unexpected input: '*'`, regardless of whether the
subworkflow/module under test is strict-clean. This is an nf-test tooling limitation
([askimed/nf-test#326](https://github.com/askimed/nf-test/issues/326), open/unresolved as of
2026-07-22), not a bug in sash. **Workaround:** force the legacy parser for nf-test runs only —

```
NXF_SYNTAX_PARSER=v1 nf-test test
```

— sash's own `.nf` files stay strict-syntax-clean; only nf-test's own wrapper generation needs the
legacy parser. Revisit this pin once nf-test ships a fix for #326.

## Tier 3 — Nextflow stack on AWS

Validate against real cloud data before merge. Covers the standard DRAGEN+OA path, the OA-only
path, the hypermutated edge case, and release regression.

Ready-to-run samplesheets live in [`assets/test/`](../assets/test). All input S3 paths were
verified readable on 2026-07-22 via `umccr-prod-operator` / `umccr-development-ro`.

Pipeline runs in this tier go through the UMCCR
[`nextflow-stack`](https://github.com/umccr/nextflow-stack) CDK repo, which deploys sash as an AWS
Batch job:

- `application/pipeline-stacks/sash/stack.ts` — the Batch job definition/queue mapping and the
  `batch_job_submission` Lambda that dispatches jobs.
- `application/pipeline-stacks/sash/assets/run.sh` (and `run-v2.sh`) — runs inside the Batch
  container: stages ICA/GDS input data to S3 via `rclone`, builds a samplesheet, injects the Batch
  instance role into `nextflow_aws.template.config`, then runs
  `nextflow run software/sash/main.nf -profile docker`, uploading results to S3 on exit.
  Supports `--stub_run` (fast structural check, same idea as tier 1's `-stub-run`) and
  `--resume_nextflow_dir` (faster iteration by resuming a prior work dir) — use these for quick
  turnaround instead of always running the full pipeline end-to-end.

### Test samples

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

### Run a validation

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

### Choosing an OA-only sample

The OA-only path expects the oncoanalyser v2.2.0 layout: `esvee/`,
`sage_calling/somatic/{tumor_id}.sage.somatic.vcf.gz`, and `pave/{tumor_id}.pave.germline.vcf.gz`
(tumor-id prefix). Older OA outputs use `gridss`/`gripss`, `sage/`, and germline PAVE named
`*.sage.germline.pave.vcf.gz` — `prepare_input.nf` cannot detect these. Confirm `esvee/` and
`pave/{tumor}.pave.germline.vcf.gz` exist before picking a sample, or re-run OA at v2.2.0.

### Hypermutated edge case (>450k variants)

The PCGR-skip / VCF2MAF-skip logic (sash #59, tracked by #52) needs a genuinely hypermutated
sample. This is a bolt/PCGR problem, not sash-only: the failure it guards against is an OOM
inside PCGR (which runs in bolt), recorded for SBJ02862 in
[umccr/biodaily#199](https://github.com/umccr/biodaily/issues/199). bolt #33/#36 carry the
matching fixes.

- **Smoke test:** a synthetic 550k-variant VCF exercises the skip branch directly and confirms
  the sample exits 0 without a full run. Generate one with a script that emits >450k PASS
  records, then feed it to bolt's PCGR prep.
- **Real-data sign-off:** the known real hypermutators (L2100242, 595k PASS; SBJ02862 /
  L2201449, the OOM case) roll off the prod cache over time. Re-stage one from the portal
  (`api/v1/s3/attributes?portalRunId=…`); a DRAGEN+OA re-run may be needed before sash.

### Release regression

Compare a new release against the previous one with the Sash Regression Service
([`umccr/service-sash-regression`](https://github.com/umccr/service-sash-regression)) on case
SEQC-II-medium (L2301218/L2301217). The 0.6.4 baseline is already wired. This replaces manual
output diffing.

The prior 0.7.0-vs-0.6.x validation record — 11 samples, CPSR/PCGR concordance, tier analysis,
and the validated sash+bolt version table — is
[umccr/biodaily#199](https://github.com/umccr/biodaily/issues/199). Use it as the concordance
baseline rather than re-running those comparisons.

## Tier 4 — ICA (production gate)

sash and oncoanalyser are deployed and triggered via **ICA** (Illumina Connected Analytics) in
production, integrated with OrcaBus EventBridge events. A release candidate goes through a 4-stage
gate before it's trusted on real clinical data:

1. **Dev** — ICA dev workspace, arbitrary test runs, fastest iteration.
2. **Prod-orchestration / testdata** — promoted into the production orchestration layer but run
   against non-clinical test data (source: `github.com/umccr/research-projects`,
   `hofmann-lab/testdata`), proving the OrcaBus wiring end-to-end without clinical risk.
3. **Prod-orchestration / clinical** — same production orchestration path, now against real
   clinical samples.
4. **Side-by-side curator review** — a curator compares the new release's output against the
   previous version's for the same sample(s) and signs off before the version becomes the
   default.

Each promotion is a deliberate, reviewed step — a release doesn't reach stage 4 by default. This
process previously only existed as a personal note; it's captured here so it's discoverable
without depending on any one person's vault. Long-term goal (aspirational, not yet implemented):
automate the testdata stage as a CI gate so stage 2 runs automatically on every release candidate.

## Variant monitoring (NATA)

Germline variant-drift monitoring on GIAB control libraries has its own SOPs in
[`umccr/accreditation`](https://github.com/umccr/accreditation/tree/main/wgs-variant-monitoring-workflow/docs)
(`RUNBOOK.md`, `SOP_OUTLIER_REVIEW.md`, `SOP_OUT_OF_LIMITS_REVIEW.md`).

## Sign-off criteria

- somatic PASS concordance ≥ 81.2% vs the 0.6.4 baseline
- PURPLE purity and ploidy identical to baseline
- all three cancer reports render without error

**Refactor-only branches** (no scientific-logic change — e.g. the Nextflow version bump / strict
syntax migration) are held to a stricter bar than the above: outputs must be **content-identical**
to the same inputs run on `main`, not just above the concordance threshold. A ≥81.2% bar is
designed to tolerate expected drift from a behavior-changing release (new tool versions, new
logic) — it would silently mask a real regression introduced by a pure refactor. For this kind of
branch: run the same tier-3 samples on both `main` and the branch, diff outputs field-by-field
(decompressed content, not raw file hashes — gzip metadata/timestamps differ), and treat any
non-metadata diff as a blocker until it's root-caused.
