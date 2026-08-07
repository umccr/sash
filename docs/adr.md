# Architecture Decision Records

---

## ADR #1: PCGR Variant Limit Handling for Large Variant Sets

**Status**: Implemented (bolt ≥ 0.3.1, sash ≥ 0.7.0)  
**Date**: 2024-11-07 (original), 2026-07-22 (updated)  
**Deciders**: Oliver Hofmann, Stephen Watts, Quentin Clayssen  
**Technical Story**: PCGR cannot safely process ≥ 500,000 variants; hypermutated samples regularly exceed this limit.

### Context

[PCGR](https://sigven.github.io/pcgr/) has a 500,000 variant threshold that triggers two silent failure modes:

1. **Python side** (`pcgr/main.py`, `pcgr_vars.MAX_VARIANTS_FOR_REPORT = 500_000`):  
   When input > 500k, PCGR silently drops intergenic, intronic, upstream_gene, and downstream_gene variants. TMB becomes an underestimate; variant selection is indiscriminate.

2. **R side** (`pcgrr/R/main.R` ~line 954):  
   If variants are still ≥ 500k after step 1, the HTML report is **silently not generated** — no error, just no output file.

Neither triggers a hard failure. Both are silent. A sample could pass through bolt, enter PCGR, and produce either a misleading report (missing variants, wrong TMB) or no report at all — with no error in the logs.

### Decision

Set `MAX_SOMATIC_VARIANTS = 450_000` in bolt (`bolt/common/constants.py`). Bolt applies clinically-prioritised tiered filtering via `select_pcgr_variants()` to bring variant counts below this cap before sending to PCGR.

**Current implementation (bolt ≥ 0.3.1):**

1. Variants with clinical retention markers (SAGE_HOTSPOT, PANEL ±2kb, PCGR_MUTATION_HOTSPOT, HMF_HOTSPOT) bypass filtering entirely.
2. Remaining variants are classified by `(tier, impact, region)` into ordered categories.
3. Categories are dropped from lowest priority (NONCODING_INTERGENIC) to highest (TIER_1) until the total count ≤ 450k.
4. If retained variants alone exceed 450k (filtering cannot help), bolt raises `RuntimeError` which sash catches — PCGR is skipped entirely for that sample.

**Why 450k, not 500k?**

- Bolt's filter is coarse-grained (drops whole categories at once). Output is always ≤ `MAX_SOMATIC_VARIANTS`, but a sample at exactly 500k would trigger pcgrr's `< 500000` check and the HTML report would be silently skipped.
- The 50k margin ensures PCGR never triggers its own internal filtering, so bolt retains full control over clinical variant prioritisation.
- See also: `bolt/docs/adr/001-max-somatic-variants-450k.md` for the detailed boundary analysis.

**When PCGR is skipped:**

- Sash marks PCGR outputs as `optional: true` in the BOLT_SMLV_SOMATIC_REPORT module.
- VCF2MAF is skipped when its input VCF is absent.
- The cancer report (gpgr) is unaffected — it does not consume PCGR output.
- Real case: L2100242, 595k PASS variants (sash #52).

### Consequences

Positive:

- Bolt controls which variants reach PCGR, preserving clinical priority.
- PCGR's own filtering never fires — TMB is accurate, HTML always generated.
- For extreme hypermutators, graceful skip instead of silent failure.
- Full filtered VCF is published independently regardless of PCGR skip.

Negative:

- **Coverage loss in PCGR report**: for hypermutated samples, lower-priority variant categories are absent from the PCGR HTML. The published VCF is unaffected.
- **No PCGR output for extreme cases**: samples where retained variants alone > 450k get no PCGR HTML, no MAF. This is intentional — PCGR would produce a misleading report anyway.

### References

- bolt: `bolt/common/constants.py` — `MAX_SOMATIC_VARIANTS = 450_000`
- bolt: `bolt/workflows/smlv_somatic/report.py` — `select_pcgr_variants()`
- bolt: `bolt/docs/adr/001-max-somatic-variants-450k.md` — detailed boundary analysis
- PCGR Python filter: `sigven/pcgr`, `pcgr/main.py` ~line 559
- pcgrr HTML skip: `sigven/pcgr`, `pcgrr/R/main.R` ~line 954
- sash #52, bolt #26, bolt #35

### The `Hypermutated` flag in the Cancer Report is a coarser, separate signal

The `Hypermutated: TRUE/FALSE` row curators see in the gpgr Cancer Report QC table is **not** the same thing as the PCGR tiered-filtering cascade described above. A `TRUE` result does not by itself mean clinically-prioritised filtering (steps 2–3 of the cascade) actually ran against that sample's variants.

**How the flag is computed** (bolt ≥ 0.3.0, `bolt/workflows/smlv_somatic/report.py::count_variant_process()`):

```python
counts['is_hypermutated'] = counts['dragen'] > constants.MAX_SOMATIC_VARIANTS
```

`dragen` is the raw DRAGEN variant count — **PASS and non-PASS included** (only novel SAGE calls are excluded). gpgr (≥ commit `f51919d`) reads this flag directly rather than re-deriving it; this replaced an earlier, less reliable inference based on `sage != annotated` count divergence.

**Why a `TRUE` flag doesn't necessarily mean TMB was affected.** The variant-count-limiting cascade in `bolt smlv_somatic annotate` (`annotate.py::select_variants()`) runs independently of, and _before_, `select_pcgr_variants()`. It tries progressively more aggressive reductions and returns as soon as the count is under `MAX_SOMATIC_VARIANTS`:

1. Exclude non-PASS variants (retaining hotspots) — FILTER tag `max_variants_non_pass`. If this alone brings the count under the limit, **the cascade stops here**.
2. Exclude common gnomAD population variants — FILTER tag `max_variants_gnomad`.
3. Restrict to the cancer gene panel region — FILTER tag `max_variants_non_cancer_genes`.

Because `is_hypermutated` is derived purely from the raw (PASS + non-PASS) count, a sample with a lot of non-PASS noise (e.g. FFPE artifacts) can trip the flag while only needing step 1 above — i.e. ordinary PASS-only filtering, identical in effect to what any non-hypermutated sample gets via the separate `BOLT_SMLV_SOMATIC_FILTER` module. In that case PCGR annotates the full PASS(+hotspot) variant set and the resulting TMB is not an underestimate; no cancer-gene-panel-only restriction occurred.

**Real case, verified directly from the published VCF** (2026-08-03, `L2600475__L2600476`): `dragen: 531,556; sage: 531,556; annotated: 163,996; filter_pass: 62,883`. Pulling `smlv_somatic/annotate/L2600475.annotations.vcf.gz` and counting FILTER tags confirms this precisely — `max_variants_non_pass` appears on 367,560 records (531,556 − 367,560 = 163,996 = `annotated`), and `max_variants_gnomad`/`max_variants_non_cancer_genes` appear on **zero** records. The cascade stopped after step 1. The sample was flagged hypermutated, but PCGR's TMB was computed from PASS-only variants with no region-based restriction — i.e. the flag was a false positive with respect to TMB accuracy.

**Caveat for curators:** the flag alone cannot distinguish "raw non-PASS noise tripped the threshold, only ordinary PASS filtering applied" from "PCGR received a clinically-prioritised, region-restricted subset". Confirming which occurred currently requires checking the `annotated`/`filter_pass` counts in `{tumor_id}.somatic.variant_counts_process.json`, or grepping the FILTER column of `{tumor_id}.annotations.vcf.gz` for `max_variants_gnomad`/`max_variants_non_cancer_genes` as above — there is no single boolean for this yet. Changing `is_hypermutated` itself to be based on a PASS-only count (e.g. `filter_pass`, already computed in the same function) has been raised as a follow-up fix in `bolt` (Slack `#genomics-technical`, 2026-08-03 thread on `C025TLC7D`).

---

## ADR #2: SV Caller Replacement — GRIDSS/GRIPSS → eSVee

**Status**: Implemented
**Date**: 2025-06-04 (v0.6.0)
**Deciders**: Oliver Hofmann, Stephen Watts, Quentin Clayssen
**Technical Story**: HMF replaced their GRIDSS/GRIPSS SV calling stack with eSVee; sash follows the HMF tool ecosystem.

### Context

The original sash SV calling chain used GRIDSS (assembly-based SV detection) followed by GRIPSS (filtering and classification). HMF developed eSVee as a purpose-built short-read SV caller that consolidates prep, assembly, depth-counting, and filtering into a single tool with direct PURPLE integration. Staying in sync with HMF tool versions is important for reference data compatibility and downstream report consistency (LINX, PURPLE driver catalogue).

### Decision

Replace GRIDSS and GRIPSS with eSVee in [v0.6.0](https://github.com/umccr/sash/blob/main/CHANGELOG.md). sash imports the eSVee prep outputs from oncoanalyser (`${tumor_id}.esvee.ref_depth.vcf.gz` and the `esvee/` directory) and runs only the eSVee caller step ([modules/local/esvee/call/main.nf](modules/local/esvee/call/main.nf)), inheriting all upstream prep and depth annotation from oncoanalyser.

### Consequences

Positive:

- Single-tool SV calling stack actively maintained and tested by HMF.
- Direct compatibility with PURPLE 4.x, LINX 2.x, and current HMF reference data.
- `sv_prep_blacklist.38.bed` replaces the former `gridss_region_blocklist`, aligning blocklist management with HMF conventions.

Negative / Breaking:

- **Metric renames**: Split Read (SR) → Split Fragments (SF); Paired Read (PR) → Discordant Fragments (DF). Downstream tools or reports referencing SR/PR by name will break.
- **LINX v1.25 → v2.0 required**: eSVee output format is incompatible with LINX v1.x.
- **CHORD removed**: CHORD depended on the GRIDSS SV signature and was removed alongside GRIDSS.
- **Kataegis module removed**: Kataegis detection relied on GRIDSS output and was removed at the same time.

### Links

- [v0.6.0 changelog](https://github.com/umccr/sash/blob/main/CHANGELOG.md)
- [eSVee README](https://github.com/hartwigmedical/hmftools/tree/master/esvee)
- [modules/local/esvee/call/main.nf](modules/local/esvee/call/main.nf)

---

## ADR #3: PURPLE Germline Enrichment Disabled

**Status**: Implemented (limitation acknowledged)
**Date**: ongoing
**Deciders**: Stephen Watts
**Technical Story**: PURPLE supports germline variant enrichment to improve purity/ploidy fitting, but the required input format is incompatible with DRAGEN germline VCFs.

### Context

PURPLE's germline enrichment pipeline (`GermlineVariantEnrichment`, `GermlinePurityEnrichment`, `GermlineGenotypeEnrichment`) expects tumor allele depth (AD) and genotype (GT) fields in the germline VCF. These fields are added by SAGE when it jointly calls tumor and normal. DRAGEN germline VCFs do not carry tumor-sample FORMAT columns, so the required fields are absent.

### Decision

Pass an empty list in place of the germline VCF when calling `PURPLE_CALLING` ([workflows/sash.nf:258–259](workflows/sash.nf#L258)):

```groovy
// ch_smlv_germline_out,  // disabled — see below
ch_smlv_germline_out.map { meta, vcf -> return [meta, []] },
```

This disables germline enrichment entirely. PURPLE still runs and produces purity/ploidy estimates using somatic SNVs, copy number, and BAF data from AMBER/COBALT.

### Consequences

Positive:

- No crash or silent data corruption from passing an incompatible germline VCF to PURPLE.

Negative:

- PURPLE's purity/ploidy fitting operates without the germline variant constraint that SAGE-based pipelines benefit from. May reduce accuracy in edge cases (e.g. samples with low somatic variant counts or unusual ploidy).
- The limitation is silent at runtime; nothing in the pipeline output signals that germline enrichment is inactive.

### Remaining Challenges

- Evaluate whether a DRAGEN → SAGE germline re-call step could produce a compatible germline VCF to re-enable enrichment.
- Add an explicit runtime warning when germline enrichment is bypassed.

### Links

- [workflows/sash.nf:254–259](workflows/sash.nf#L254)
- [PURPLE GermlineVariantEnrichment.java](https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlineVariantEnrichment.java#L30)
- [PURPLE GermlinePurityEnrichment.java](https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlinePurityEnrichment.java#L50)
- [PURPLE GermlineGenotypeEnrichment.java](https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlineGenotypeEnrichment.java#L63)

---

## ADR #4: PAVE Runs Exclusively for PURPLE Driver Catalogue

**Status**: Implemented
**Date**: 2025-09-16 (v0.6.1)
**Deciders**: Stephen Watts
**Technical Story**: PURPLE requires PAVE-annotated somatic variants to produce its driver catalogue; sash clinical filtering uses bolt, not PAVE.

### Context

PURPLE's driver catalogue (identifying driver genes with somatic SNVs/indels) depends on PAVE functional annotations (consequence, impact, coding effect) being present in the somatic VCF. Without PAVE, PURPLE produces an incomplete driver catalogue. At the same time, sash uses bolt for all clinical filtering and reporting; PAVE annotations are neither required nor consumed by the bolt pipeline or the cancer report.

### Decision

Run PAVE on the bolt-filtered VCF ([workflows/sash.nf:171–185](workflows/sash.nf#L171)) and feed its output exclusively to `PURPLE_CALLING`. The PAVE output VCF is not published as a standalone artifact and does not flow into the clinical report chain.

Before running PAVE, MNV-component variants are excluded ([modules/local/pave/somatic/main.nf:32](modules/local/pave/somatic/main.nf#L32)):

```sh
bcftools view --exclude 'INFO/MNVTAG!="."' ... ${vcf}
```

This prevents PAVE from annotating individual substitution components already represented as a combined MNV record, which would cause double-counting in PURPLE's driver catalogue. Added in [v0.6.1](https://github.com/umccr/sash/blob/main/CHANGELOG.md) ([sash#19](https://github.com/umccr/sash/issues/19)).

### Consequences

Positive:

- PURPLE produces a complete driver catalogue including SNV/indel drivers.
- PAVE runs on the already-filtered VCF, keeping its input size small.

Negative:

- **PAVE annotations absent from the published VCF**: the clinical output (`smlv_somatic/filter/{tumor_id}.pass.vcf.gz`) does not carry PAVE consequence fields. Any downstream tool that needs PAVE annotations must re-run PAVE independently.
- The pipeline runs PAVE on a temporary intermediate that is discarded; compute cost is paid but the output is invisible to the analyst.

### Links

- [workflows/sash.nf:171](workflows/sash.nf#L171)
- [modules/local/pave/somatic/main.nf](modules/local/pave/somatic/main.nf)
- [sash#19: PAVE MNV filtering discussion](https://github.com/umccr/sash/issues/19)
- [v0.6.1 changelog](https://github.com/umccr/sash/blob/main/CHANGELOG.md)
