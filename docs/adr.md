# Architecture Decision Records

---

## ADR #1: PCGR Variant Limit Handling for Large Variant Sets

**Status**: Partially Implemented
**Date**: 2024-11-07
**Deciders**: Oliver Hofmann, Stephen Watts, Quentin Clayssen
**Technical Story**: PCGR cannot process more than 500,000 variants per run; hypermutated samples regularly exceed this limit.

### Context

[PCGR](https://sigven.github.io/pcgr/) has a hard variant processing limit of 500,000 variants per run. Hypermutated samples in the sash workflow often produce annotated VCFs that exceed this limit, causing PCGR to truncate or fail. Two approaches were considered: chunking the VCF into parallel PCGR runs, or progressively reducing the variant set upstream of PCGR.

### Decision

Original plan (not implemented): split VCF files into chunks of ≤ 500,000 variants, process each chunk in parallel through PCGR, then merge annotated outputs.

Current implementation — bolt applies a `select_variants` function that progressively excludes lower-priority variants until the count falls below 500,000 ([bolt/annotate.py:182–200](https://github.com/umccr/bolt/blob/v0.2.18/bolt/workflows/smlv_somatic/annotate.py#L182)):

1. If total ≤ 500,000: pass all variants unchanged
2. Remove non-PASS variants; if now ≤ 500,000: stop
3. Remove population-common variants (gnomAD AF ≥ 1%); if now ≤ 500,000: stop
4. Remove non-cancer-gene variants; pass whatever remains to PCGR

This reduction affects only the PCGR report input. The full filtered VCF is published independently. TMB and MSI are derived from PURPLE and are unaffected by PCGR input size.

Open TODOs remain: logging of the reduction steps and tempdir cleanup ([bolt/annotate.py:177–179](https://github.com/umccr/bolt/blob/v0.2.18/bolt/workflows/smlv_somatic/annotate.py#L177)).

A "hypermutated" flag was added to the gpgr cancer report to surface this condition to analysts. It was reverted in [v0.6.4](https://github.com/umccr/sash/blob/main/CHANGELOG.md) due to incorrect triggering conditions ([gpgr#101](https://github.com/umccr/gpgr/pull/101), [bolt#28](https://github.com/umccr/bolt/pull/28)). A redesigned flag may be reintroduced in a future release.

### Consequences

Positive:

- No merge step required; a single PCGR run is preserved per sample.
- The exclusion order ensures the most clinically relevant variants are the last to be removed.

Negative:

- **Coverage loss in PCGR report**: for extremely high variant counts, non-cancer-gene variants may be absent from the PCGR HTML even when they pass all sash filters. The published VCF is unaffected.
- **No parallelization**: the original goal of parallel PCGR runs is not met; PCGR still runs serially on a single reduced input.

### Remaining Challenges

- True VCF chunking and parallel PCGR execution remain viable as a future improvement.
- The `select_variants` logging and tempdir TODOs need resolution.
- Reintroduction of a correctly scoped hypermutated indicator in the cancer report.

### Links

- [bolt/annotate.py — select_variants](https://github.com/umccr/bolt/blob/v0.2.18/bolt/workflows/smlv_somatic/annotate.py#L182)
- [PCGR — large input sets](https://sigven.github.io/pcgr/articles/running.html#large-input-sets-vcf)
- [v0.6.4 — hypermutated flag revert](https://github.com/umccr/sash/blob/main/CHANGELOG.md)
- [gpgr#101](https://github.com/umccr/gpgr/pull/101), [bolt#28](https://github.com/umccr/bolt/pull/28)

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
- The limitation is silent at runtime — nothing in the pipeline output signals that germline enrichment is inactive.

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

PURPLE's driver catalogue (identifying driver genes with somatic SNVs/indels) depends on PAVE functional annotations (consequence, impact, coding effect) being present in the somatic VCF. Without PAVE, PURPLE produces an incomplete driver catalogue. At the same time, sash uses bolt for all clinical filtering and reporting — PAVE annotations are neither required nor consumed by the bolt pipeline or the cancer report.

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
- The pipeline runs PAVE on a temporary intermediate that is discarded — compute cost is paid but the output is invisible to the analyst.

### Links

- [workflows/sash.nf:171](workflows/sash.nf#L171)
- [modules/local/pave/somatic/main.nf](modules/local/pave/somatic/main.nf)
- [sash#19 — PAVE MNV filtering discussion](https://github.com/umccr/sash/issues/19)
- [v0.6.1 changelog](https://github.com/umccr/sash/blob/main/CHANGELOG.md)
