# Sash Output

## Introduction

This document outlines the key results and files produced by the UMCCR SASH (post-processing WGS tumor/normal) pipeline. After a run, the pipeline organizes output files by analysis module under a directory for each tumor/normal pair (identified by run ID and sample names). The main outputs include annotated variant reports for somatic and germline variants, copy number and structural variant analyses, and a comprehensive MultiQC report for quality control. All paths below are relative to the top-level results directory of a given run.
## Pipeline overview

- [Sash Output](#sash-output)
  - [Introduction](#introduction)
  - [Pipeline overview](#pipeline-overview)
  - [Directory Structure](#directory-structure)
  - [Summary](#summary)
  - [Workflows](#workflows)
    - [Somatic Small Variants](#somatic-small-variants)
      - [General](#general)
      - [Summary](#summary-1)
    - [Details](#details)
  - [bolt smlv somatic rescue](#bolt-smlv-somatic-rescue)
      - [BOLT\_SMLV\_SOMATIC\_ANNOTATE](#bolt_smlv_somatic_annotate)
      - [BOLT\_SMLV\_SOMATIC\_FILTER](#bolt_smlv_somatic_filter)
      - [SOMATIC\_SNV\_REPORTS](#somatic_snv_reports)
    - [Somatic Structural Variants](#somatic-structural-variants)
      - [General](#general-1)
      - [Summary](#summary-2)
      - [SV Annotation](#sv-annotation)
      - [SV Prioritization](#sv-prioritization)
    - [Germline Variants](#germline-variants)
      - [General](#general-2)
      - [Summary](#summary-3)
      - [Germline Preparation](#germline-preparation)
      - [Germline Reports](#germline-reports)
    - [Reports](#reports)
      - [Cancer Report](#cancer-report)
      - [LINX Reports](#linx-reports)
      - [PURPLE Reports](#purple-reports)
      - [PCGR Reports](#pcgr-reports)
      - [SIGRAP Reports](#sigrap-reports)
      - [CPSR Reports](#cpsr-reports)
      - [MultiQC Reports](#multiqc-reports)

## Directory Structure

```bash
[RUN_ID]/[sample]/
├── <sample>.cancer_report.html
├── <sample>.cpsr.html
├── <sample>.pcgr.html
├── <sample>_linx.html
├── <sample>.multiqc.html
├── cancer_report/
│ ├── img/
│ └── cancer_report_tables/
│   ├── hrd/
│   ├── json/
│   ├── purple/
│   └── sigs/ 
├── linx/
│ ├── germline_annotations/
│ ├── somatic_annotations/
│ └── somatic_plots/
├── multiqc_data/
├── purple/
├── smlv_germline/
│ └── prepare/
│ └── report/
├── smlv_somatic/
│ ├── report/
│ ├── annotate/
│ ├── filter/
│ └── rescue/
└── sv_somatic/
   ├── annotate/
   └── prioritise/
```
## Summary

The **Sash Workflow** comprises three primary pipelines: **Somatic Small Variants**, **Somatic Structural Variants**, and **Germline Variants**. These pipelines utilize **Bolt**, a Python package designed for modular processing, and leverage outputs from the **DRAGEN Variant Caller** alongside **HMFtools in Oncoanalyser**. Each pipeline is tailored to a specific type of genomic variant, incorporating filtering, annotation, and HTML reports for research and curation.

## Workflows

### Somatic Small Variants

#### General

In the **Somatic Small Variants** workflow, variant detection is performed using the **DRAGEN Variant Caller** and **Oncoanalyser (SAGE, Purple)** outputs. It's structured into four steps: **Rescue**, **Annotation**, **Filter**, and **Report**. The final outputs include an **HTML report** summarizing the results.

#### Summary

1. **Rescue** variants using SAGE to recover low-frequency alterations in clinically important hotspots.
2. **Annotate** variants with clinical and functional information using PCGR.
3. **Filter** variants based on allele frequency (AF), supporting reads (AD), and population frequency (gnomAD AF), removing low-confidence and common variants.
4. **Report** final annotated variants in a comprehensive HTML format.

### Details

## bolt smlv somatic rescue

<details markdown="1">
<summary>Output files</summary>

- `smlv_somatic/rescue/`
  - `<tumor_id>.rescued.vcf.gz`: Rescued somatic VCF file containing previously filtered variants at known hotspots.

</details>

The `BOLT_SMLV_SOMATIC_RESCUE` process rescues somatic variants using the BOLT tool. The output includes the rescued VCF file that recovers potentially important variants that may have been filtered in earlier steps due to borderline quality metrics.

#### BOLT_SMLV_SOMATIC_ANNOTATE

<details markdown="1">
<summary>Output files</summary>

- `smlv_somatic/annotate/`
  - `<tumor_id>.annotations.vcf.gz`: Annotated somatic VCF file with functional and clinical annotations.

</details>

The `BOLT_SMLV_SOMATIC_ANNOTATE` process annotates somatic variants using the BOLT tool. The output includes the annotated VCF file enriched with gene information, variant effect predictions, and other annotations to aid in variant interpretation.

#### BOLT_SMLV_SOMATIC_FILTER

<details markdown="1">
<summary>Output files</summary>

- `smlv_somatic/filter/`
  - `<tumor_id>.filters_set.vcf.gz`: VCF file with filters set but all variants retained.
  - `<tumor_id>.pass.vcf.gz`: Filtered somatic VCF file containing only PASS variants.
  - `<tumor_id>.pass.vcf.gz.tbi`: Index file for the filtered VCF.

</details>

The `BOLT_SMLV_SOMATIC_FILTER` process filters somatic variants using the BOLT tool. The output includes both a VCF with all variants but filter tags applied, and a filtered VCF containing only variants that pass all quality filters.

#### SOMATIC_SNV_REPORTS

<details markdown="1">
<summary>Output files</summary>

- `smlv_somatic/report/`
  - `<tumor_id>.somatic.bcftools_stats.txt`: Statistical summary of somatic variants.
  - `<tumor_id>.somatic.variant_counts_process.json`: Variant count metrics at each processing step.
  - `<tumor_id>.somatic.variant_counts_type.yaml`: Variant counts by variant type.
  - `af_tumor.txt`: Information about variant allele frequencies in tumor.
  - `af_tumor_keygenes.txt`: Variant allele frequencies in key cancer-related genes.
  - `pcgr/`: Directory containing PCGR report outputs.

</details>

The reporting process generates statistical summaries and specialized reports for somatic SNVs, including PCGR HTML reports for clinical interpretation.

### Somatic Structural Variants

#### General

The **Somatic Structural Variants** workflow identifies and analyzes large genomic rearrangements such as deletions, duplications, inversions, and translocations. It processes outputs from GRIDSS, PURPLE, and LINX to provide comprehensive SV analysis.

#### Summary

1. **Annotate** SVs with gene context and potential functional impacts.
2. **Prioritize** SVs based on cancer relevance and gene disruption potential.
3. **Report** clinically relevant SVs with gene fusion predictions and visualization.

#### SV Annotation

<details markdown="1">
<summary>Output files</summary>

- `sv_somatic/annotate/`
  - `<tumor_id>.annotated.vcf.gz`: Annotated structural variant VCF file.

</details>

This process adds gene annotations and functional impact predictions to structural variants. The annotated VCF contains information about genes affected by breakpoints, potential fusion events, and other biologically relevant details.

#### SV Prioritization

<details markdown="1">
<summary>Output files</summary>

- `sv_somatic/prioritise/`
  - `<tumor_id>.cnv.prioritised.tsv`: Prioritized copy number variations in tabular format.
  - `<tumor_id>.sv.prioritised.tsv`: Prioritized structural variants in tabular format.
  - `<tumor_id>.sv.prioritised.vcf.gz`: Prioritized structural variants in VCF format.

</details>

This process ranks structural variants based on their potential clinical relevance, creating filterable lists for review. It separately handles copy number variations and other structural variants for easier interpretation.

### Germline Variants

#### General

The **Germline Variants** workflow analyzes inherited variants from the normal sample to identify potential cancer predisposition genes and variants that may influence treatment decisions.

#### Summary

1. **Prepare** germline variants from DRAGEN normal sample outputs.
2. **Report** potentially actionable germline variants through CPSR.

#### Germline Preparation

<details markdown="1">
<summary>Output files</summary>

- `smlv_germline/prepare/`
  - `<normal_id>.prepared.vcf.gz`: Prepared germline VCF file for annotation.

</details>

This process prepares germline variants for downstream annotation and reporting. It applies normalization, left-alignment, and other preprocessing steps to ensure consistent variant representation.

#### Germline Reports

<details markdown="1">
<summary>Output files</summary>

- `smlv_germline/report/`
  - `<normal_id>.annotations.vcf.gz`: Annotated germline VCF.
  - `<normal_id>.germline.bcftools_stats.txt`: Statistical summary of germline variants.
  - `<normal_id>.germline.variant_counts_type.yaml`: Variant counts by type.
  - `cpsr/`: Directory containing CPSR outputs.
    - `<normal_id>.cpsr.grch38.json.gz`: Structured CPSR data.
    - `<normal_id>.cpsr.grch38.pass.tsv.gz`: Filtered CPSR variants in tabular format.
    - `<normal_id>.cpsr.grch38.snvs_indels.tiers.tsv`: Tiered variants by clinical significance.
    - Other CPSR output files.

</details>

The germline reporting process focuses on identifying variants in cancer predisposition genes and producing a comprehensive CPSR (Cancer Predisposition Sequencing Reporter) report.

### Reports

#### Cancer Report

<details markdown="1">
<summary>Output files</summary>

- `<tumor_id>.cancer_report.html`: Main cancer report HTML file.
- `cancer_report/`
  - `<tumor_id>.snvs.normalised.vcf.gz`: Normalized SNVs used in the report.
  - `img/`: Images used in the cancer report.
  - `cancer_report_tables/`: Tabular data supporting the report.
    - `<normal_id>_<tumor_id>-qc_summary.tsv.gz`: Quality control summary.
    - `<normal_id>_<tumor_id>-report_inputs.tsv.gz`: Report configuration inputs.
    - `hrd/`: Homologous Recombination Deficiency analysis from multiple methods.
    - `json/`: JSON-formatted report data.
    - `purple/`: Copy number information.
    - `sigs/`: Mutational signature analysis (SBS, DBS, indels).

</details>

The cancer report integrates findings from all analysis modules into a comprehensive HTML report for clinical interpretation. It includes tumor characteristics, key somatic alterations, mutational signatures, and therapy recommendations.

#### LINX Reports

<details markdown="1">
<summary>Output files</summary>

- `<tumor_id>_linx.html`: LINX visualization report.
- `linx/`
  - `germline_annotations/`: Germline structural variant analysis.
    - `<tumor_id>.linx.germline.breakend.tsv`: Germline breakend annotations.
    - `<tumor_id>.linx.germline.clusters.tsv`: Germline SV clusters.
    - `<tumor_id>.linx.germline.disruption.tsv`: Gene disruptions by germline SVs.
    - `<tumor_id>.linx.germline.driver.catalog.tsv`: Potential driver germline SVs.
    - `<tumor_id>.linx.germline.links.tsv`: Links between germline SVs.
    - `<tumor_id>.linx.germline.svs.tsv`: Germline structural variants.
    - `linx.version`: LINX version information.
  - `somatic_annotations/`: Somatic structural variant analysis.
    - `<tumor_id>.linx.breakend.tsv`: Somatic breakend annotations.
    - `<tumor_id>.linx.clusters.tsv`: Somatic SV clusters.
    - `<tumor_id>.linx.driver.catalog.tsv`: Potential driver somatic SVs.
    - `<tumor_id>.linx.drivers.tsv`: Driver SV details.
    - `<tumor_id>.linx.fusion.tsv`: Gene fusion predictions.
    - `<tumor_id>.linx.links.tsv`: Links between somatic SVs.
    - `<tumor_id>.linx.svs.tsv`: Somatic structural variants.
    - Visualization data files (vis_*).
    - `linx.version`: LINX version information.
  - `somatic_plots/`: Visualizations of somatic structural variants.

</details>

LINX reports provide detailed analysis of structural variants, including gene fusions, disruptions, and visualization of complex rearrangements. The HTML visualization report offers interactive exploration of structural variants and their potential functional impacts.

#### PURPLE Reports

<details markdown="1">
<summary>Output files</summary>

- `purple/`
  - `<tumor_id>.purple.cnv.gene.tsv`: Gene-level copy number variations.
  - `<tumor_id>.purple.cnv.somatic.tsv`: Segment-level somatic copy number variations.
  - `<tumor_id>.purple.driver.catalog.germline.tsv`: Potential germline driver variants.
  - `<tumor_id>.purple.driver.catalog.somatic.tsv`: Potential somatic driver variants.
  - `<tumor_id>.purple.germline.deletion.tsv`: Germline deletion information.
  - `<tumor_id>.purple.purity.range.tsv`: Range of possible purity values.
  - `<tumor_id>.purple.purity.tsv`: Tumor purity, ploidy, and microsatellite status.
  - `<tumor_id>.purple.qc`: Quality control metrics.
  - `<tumor_id>.purple.segment.tsv`: Genomic segmentation data.
  - `<tumor_id>.purple.somatic.clonality.tsv`: Clonality analysis of somatic variants.
  - `<tumor_id>.purple.somatic.hist.tsv`: Somatic variant histograms.
  - `<tumor_id>.purple.somatic.vcf.gz`: Somatic variants VCF with copy number annotations.
  - `<tumor_id>.purple.sv.germline.vcf.gz`: Germline structural variants.
  - `<tumor_id>.purple.sv.vcf.gz`: Somatic structural variants.
  - `circos/`: Circos visualization data files.
    - `<normal_id>.ratio.circos`: Normal sample coverage ratio data.
    - `<tumor_id>.baf.circos`: B-allele frequency data.
    - `<tumor_id>.cnv.circos`: Copy number data.
    - `<tumor_id>.indel.circos`: Indel visualization data.
    - `<tumor_id>.link.circos`: SV links visualization data.
    - `<tumor_id>.snp.circos`: SNP visualization data.
    - Configuration and input files for Circos.
  - `plot/`: Additional visualization data.
  - `purple.version`: PURPLE version information.

</details>

PURPLE reports provide copy number analysis, tumor purity estimation, and whole genome doubling assessment. The circos directory contains data for generating circular genome plots that visualize genomic alterations across the entire genome.

#### PCGR Reports

<details markdown="1">
<summary>Output files</summary>

- `<tumor_id>.pcgr.html`: PCGR HTML report.
- `smlv_somatic/report/pcgr/`
  - `<tumor_id>.pcgr_acmg.grch38.flexdb.html`: Flexible database PCGR report.
  - `<tumor_id>.pcgr_acmg.grch38.json.gz`: Structured PCGR data in JSON format.
  - `<tumor_id>.pcgr_acmg.grch38.mp_input.vcf.gz`: Input VCF for mutational pattern analysis.
  - `<tumor_id>.pcgr_acmg.grch38.mutational_signatures.tsv`: Mutational signature analysis.
  - `<tumor_id>.pcgr_acmg.grch38.pass.tsv.gz`: Filtered variants in tabular format.
  - `<tumor_id>.pcgr_acmg.grch38.pass.vcf.gz`: Filtered variants in VCF format.
  - `<tumor_id>.pcgr_acmg.grch38.snvs_indels.tiers.tsv`: Tiered variants by clinical significance.
  - `<tumor_id>.pcgr_acmg.grch38.vcf.gz`: All variants in VCF format.
  - `<tumor_id>.pcgr_config.rds`: PCGR configuration.

</details>

PCGR (Personal Cancer Genome Reporter) reports provide clinical interpretation of somatic variants, including therapy matches, clinical trial eligibility, and tumor mutational burden assessment.

#### SIGRAP Reports

<details markdown="1">
<summary>Output files</summary>

- `cancer_report/cancer_report_tables/sigs/`
  - `<normal_id>_<tumor_id>-dbs.tsv.gz`: Double base substitution signature analysis.
  - `<normal_id>_<tumor_id>-indel.tsv.gz`: Indel signature analysis.
  - `<normal_id>_<tumor_id>-snv_2015.tsv.gz`: SNV signature analysis using 2015 signatures.
  - `<normal_id>_<tumor_id>-snv_2020.tsv.gz`: SNV signature analysis using 2020 signatures.
- `cancer_report/cancer_report_tables/json/sigs/`: JSON-formatted signature data.

</details>

SIGRAP reports provide mutational signature analysis, identifying patterns associated with specific mutational processes or exposures. The pipeline analyzes single base substitutions (SBS), double base substitutions (DBS), and indel signatures using both the 2015 and 2020 reference signature sets.

#### CPSR Reports

<details markdown="1">
<summary>Output files</summary>

- `<normal_id>.cpsr.html`: CPSR HTML report.
- `smlv_germline/report/cpsr/`
  - `<normal_id>.cpsr.grch38.custom_list.bed`: Custom gene list in BED format.
  - `<normal_id>.cpsr.grch38.json.gz`: Structured CPSR data in JSON format.
  - `<normal_id>.cpsr.grch38.pass.tsv.gz`: Filtered variants in tabular format.
  - `<normal_id>.cpsr.grch38.pass.vcf.gz`: Filtered variants in VCF format.
  - `<normal_id>.cpsr.grch38.snvs_indels.tiers.tsv`: Tiered variants by clinical significance.
  - `<normal_id>.cpsr.grch38.vcf.gz`: All variants in VCF format.
  - `<normal_id>.cpsr_config.rds`: CPSR configuration.

</details>

CPSR (Cancer Predisposition Sequencing Reporter) focuses on germline variants in known cancer predisposition genes, providing a comprehensive report of inherited cancer risk variants.

#### MultiQC Reports

<details markdown="1">
<summary>Output files</summary>

- `<tumor_id>.multiqc.html`: Main MultiQC report.
- `multiqc_data/`: Supporting data for the MultiQC report.
  - `dragen_frag_len.txt`: Fragment length metrics.
  - `dragen_map_metrics.txt`: Mapping metrics.
  - `dragen_ploidy.txt`: Ploidy estimation metrics.
  - `dragen_time_metrics.txt`: Processing time metrics.
  - `dragen_trimmer_metrics.txt`: Read trimming metrics.
  - `dragen_vc_metrics.txt`: Variant calling metrics.
  - `dragen_wgs_cov_metrics.txt`: WGS coverage metrics.
  - `multiqc.log`: MultiQC log file.
  - `multiqc_bcftools_stats.txt`: BCFtools statistics.
  - `multiqc_data.json`: MultiQC data in JSON format.
  - `multiqc_general_stats.txt`: General statistics.
  - `purple.txt`: PURPLE metrics.

</details>

MultiQC aggregates quality metrics from all pipeline components into a single HTML report, providing an overview of sample quality and analysis performance.
