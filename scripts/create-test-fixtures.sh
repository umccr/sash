#!/usr/bin/env bash
# Creates minimal stub fixtures so `nextflow run . -profile test,docker -stub-run` works.
#
# Usage:
#   bash scripts/create-test-fixtures.sh
#
# This creates:
#   assets/test/stub/samples/       — mock DRAGEN+OA sample output directories
#   assets/test/stub/refdata/       — empty placeholder reference files
#   assets/test/stub/test_samplesheet.csv  — samplesheet pointing at the stubs
#
# It also overwrites conf/test.config with a test profile that points at these stubs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
STUB_DIR="$REPO_ROOT/assets/test/stub"
SAMPLE_DIR="$STUB_DIR/samples"
REFDATA_DIR="$STUB_DIR/refdata"
SAMPLE_ID="SBJ0000"
TUMOR_ID="TUMOR01"
NORMAL_ID="NORMAL01"

echo "Creating stub fixtures in $STUB_DIR"

# ---- Clean slate ----
rm -rf "$STUB_DIR"
mkdir -p "$SAMPLE_DIR" "$REFDATA_DIR"

# ============================================================================
# 1. Sample data stubs (for prepare_input.nf)
# ============================================================================
# Layout matches what dragen_oa_seqcii samplesheet expects

# --- oncoanalyser output stub ---
OA_DIR="$SAMPLE_DIR/oncoanalyser"
mkdir -p "$OA_DIR"

# AMBER
mkdir -p "$OA_DIR/amber"
touch "$OA_DIR/amber/.stub"

# COBALT
mkdir -p "$OA_DIR/cobalt"
touch "$OA_DIR/cobalt/.stub"

# eSVee
mkdir -p "$OA_DIR/esvee"
touch "$OA_DIR/esvee/${TUMOR_ID}.esvee.ref_depth.vcf.gz"
touch "$OA_DIR/esvee/.stub"

# SAGE somatic
mkdir -p "$OA_DIR/sage_calling/somatic"
touch "$OA_DIR/sage_calling/somatic/${TUMOR_ID}.sage.somatic.vcf.gz"
touch "$OA_DIR/sage_calling/somatic/${TUMOR_ID}.sage.somatic.vcf.gz.tbi"

# VirusBreakend
mkdir -p "$OA_DIR/virusbreakend"
touch "$OA_DIR/virusbreakend/.stub"

# CHORD
mkdir -p "$OA_DIR/chord"
touch "$OA_DIR/chord/${TUMOR_ID}.chord.prediction.tsv"

# --- DRAGEN somatic output stub ---
DRAGEN_SOM_DIR="$SAMPLE_DIR/dragen_somatic"
mkdir -p "$DRAGEN_SOM_DIR"
touch "$DRAGEN_SOM_DIR/${TUMOR_ID}.hard-filtered.vcf.gz"
touch "$DRAGEN_SOM_DIR/${TUMOR_ID}.hard-filtered.vcf.gz.tbi"
# Optional HRD score — include it for realism
touch "$DRAGEN_SOM_DIR/${TUMOR_ID}.hrdscore.csv"

# --- DRAGEN germline output stub ---
DRAGEN_GER_DIR="$SAMPLE_DIR/dragen_germline"
mkdir -p "$DRAGEN_GER_DIR"
touch "$DRAGEN_GER_DIR/${NORMAL_ID}.hard-filtered.vcf.gz"

echo "  Sample stubs created under $SAMPLE_DIR"

# ============================================================================
# 2. Reference data stubs (for prepare_reference.nf)
# ============================================================================
# Each path from conf/refdata.config needs a placeholder file.

# --- Genome files ---
mkdir -p "$REFDATA_DIR/genomes/GRCh38_umccr"
mkdir -p "$REFDATA_DIR/genomes/GRCh38_umccr/samtools_index/1.16"
touch "$REFDATA_DIR/genomes/GRCh38_umccr/GRCh38_full_analysis_set_plus_decoy_hla.fa"
touch "$REFDATA_DIR/genomes/GRCh38_umccr/samtools_index/1.16/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
touch "$REFDATA_DIR/genomes/GRCh38_umccr/samtools_index/1.16/GRCh38_full_analysis_set_plus_decoy_hla.fa.dict"

# --- HMF reference data ---
HMF_DIR="$REFDATA_DIR/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.2.0--3"

mkdir -p "$HMF_DIR/dna/copy_number"
touch "$HMF_DIR/dna/copy_number/AmberGermlineSites.38.tsv.gz"
touch "$HMF_DIR/dna/copy_number/GC_profile.1000bp.38.cnp"
touch "$HMF_DIR/dna/copy_number/DiploidRegions.38.bed.gz"
touch "$HMF_DIR/dna/copy_number/cohort_germline_del_freq.38.csv"

mkdir -p "$HMF_DIR/misc/cuppa"
touch "$HMF_DIR/misc/cuppa/.stub"

mkdir -p "$HMF_DIR/dna/sv"
touch "$HMF_DIR/dna/sv/sv_prep_blacklist.38.bed"
touch "$HMF_DIR/dna/sv/sgl_pon.38.bed.gz"
touch "$HMF_DIR/dna/sv/sv_pon.38.bedpe.gz"
touch "$HMF_DIR/dna/sv/repeat_mask_data.38.fa.gz"
touch "$HMF_DIR/dna/sv/gridss.properties"

mkdir -p "$HMF_DIR/rna"
touch "$HMF_DIR/rna/read_151_exp_counts.38.csv"
touch "$HMF_DIR/rna/read_100_exp_gc_ratios.38.csv"

mkdir -p "$HMF_DIR/dna/variants"
touch "$HMF_DIR/dna/variants/clinvar.38.vcf.gz"
touch "$HMF_DIR/dna/variants/KnownBlacklist.germline.38.bed"
touch "$HMF_DIR/dna/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz"
touch "$HMF_DIR/dna/variants/KnownHotspots.somatic.38.vcf.gz"
touch "$HMF_DIR/dna/variants/hmf_wgs_sage_pon_1000.38.tsv.gz"
mkdir -p "$HMF_DIR/dna/variants/gnomad"
touch "$HMF_DIR/dna/variants/gnomad/.stub"

mkdir -p "$HMF_DIR/misc/virusinterpreter"
touch "$HMF_DIR/misc/virusinterpreter/virus_reporting_db.tsv"
touch "$HMF_DIR/misc/virusinterpreter/taxonomy_db.tsv"

# --- UMCCR reference data ---
UMCCR_DIR="$REFDATA_DIR/umccr_reference_data/v2--1"

# Gene panels
mkdir -p "$UMCCR_DIR/gene_panels/somatic"
touch "$UMCCR_DIR/gene_panels/somatic/umccr_cancer_genes.tsv"
touch "$UMCCR_DIR/gene_panels/somatic/umccr_cancer_genes.gene_symbols.all.txt"
touch "$UMCCR_DIR/gene_panels/somatic/umccr_cancer_genes.gene_symbols.tsgenes.txt"
touch "$UMCCR_DIR/gene_panels/somatic/umccr_cancer_genes.gene_regions.bed"
touch "$UMCCR_DIR/gene_panels/somatic/umccr_cancer_genes.cds_regions.bed"

mkdir -p "$UMCCR_DIR/gene_panels/fusion"
touch "$UMCCR_DIR/gene_panels/fusion/hartwig_known_pairs.csv"
touch "$UMCCR_DIR/gene_panels/fusion/hartwig_promiscuous_five.csv"
touch "$UMCCR_DIR/gene_panels/fusion/hartwig_promiscuous_three.csv"

mkdir -p "$UMCCR_DIR/gene_panels/germline"
touch "$UMCCR_DIR/gene_panels/germline/umccr_predisposition_genes.gene_symbols.txt"
touch "$UMCCR_DIR/gene_panels/germline/umccr_predisposition_genes.transcript_regions.bed"

# Annotations
mkdir -p "$UMCCR_DIR/annotations/hotspots"
touch "$UMCCR_DIR/annotations/hotspots/hotspots.combined.vcf.gz"
touch "$UMCCR_DIR/annotations/.stub"

# HMFtools
mkdir -p "$UMCCR_DIR/hmftools/small_variants"
touch "$UMCCR_DIR/hmftools/small_variants/DriverGenePanel.38.tsv"
touch "$UMCCR_DIR/hmftools/small_variants/KnownHotspots.germline.38.vcf.gz"

mkdir -p "$UMCCR_DIR/hmftools/ensembl_data_cache"
touch "$UMCCR_DIR/hmftools/ensembl_data_cache/.stub"

mkdir -p "$UMCCR_DIR/hmftools/fusions"
touch "$UMCCR_DIR/hmftools/fusions/known_fusions.38.bedpe"
touch "$UMCCR_DIR/hmftools/fusions/known_fusion_data.38.csv"

mkdir -p "$UMCCR_DIR/other/purple_baf_circos"
touch "$UMCCR_DIR/other/purple_baf_circos/circos_baf.conf"
touch "$UMCCR_DIR/other/purple_baf_circos/gaps.38.txt"

# --- Misc data (tarballs and directory stubs) ---
PCGR_DIR="$REFDATA_DIR/databases/pcgr"
mkdir -p "$PCGR_DIR"
# Create minimal valid tar files
tar -cf "$PCGR_DIR/pcgr_ref_data.20250314.grch38.tgz" -T /dev/null 2>/dev/null || touch "$PCGR_DIR/pcgr_ref_data.20250314.grch38.tgz"
tar -cf "$PCGR_DIR/homo_sapiens_vep_113_GRCh38.tar.gz" -T /dev/null 2>/dev/null || touch "$PCGR_DIR/homo_sapiens_vep_113_GRCh38.tar.gz"

SNPEFF_DIR="$REFDATA_DIR/databases/snpeff/v5_1"
mkdir -p "$SNPEFF_DIR"
touch "$SNPEFF_DIR/.stub"

ONCOKB_DIR="$REFDATA_DIR/databases/oncokb/v4_12"
mkdir -p "$ONCOKB_DIR"
touch "$ONCOKB_DIR/cancerGeneList.txt"

APPRIS_DIR="$REFDATA_DIR/databases/appris/v46"
mkdir -p "$APPRIS_DIR"
touch "$APPRIS_DIR/appris_data.tsv"

PON_DIR="$REFDATA_DIR/panel_of_normals/umccr/tertiary_workflow/v2019_Feb"
mkdir -p "$PON_DIR"
touch "$PON_DIR/.stub"

echo "  Reference stubs created under $REFDATA_DIR"

# ============================================================================
# 3. Samplesheet CSV
# ============================================================================
SAMPLESHEET="$STUB_DIR/test_samplesheet.csv"
cat > "$SAMPLESHEET" << EOF
filepath,filetype,id,sample_name,subject_name
${SAMPLE_DIR}/dragen_somatic,dragen_somatic_dir,${SAMPLE_ID},${TUMOR_ID},SBJ00480
${SAMPLE_DIR}/dragen_germline,dragen_germline_dir,${SAMPLE_ID},${NORMAL_ID},SBJ00480
${SAMPLE_DIR}/oncoanalyser,oncoanalyser_dir,${SAMPLE_ID},${TUMOR_ID},SBJ00480
EOF
echo "  Samplesheet created at $SAMPLESHEET"

# # ============================================================================
# # 4. Update conf/test.config (commented out — update manually)
# # ============================================================================
# echo ""
# echo "Done. To use these fixtures, update conf/test.config with:"
# echo ""
# echo "  params.ref_data_path = '${REFDATA_DIR}'"
# echo "  params.input         = '${SAMPLESHEET}'"
# echo ""

echo ""
echo "==========================================="
echo "Fixture setup complete!"
echo "==========================================="
echo ""
echo "Stub directory: $STUB_DIR"
echo ""
echo "To use, update conf/test.config:"
echo ""
echo "  params {"
echo "    ref_data_path = '${REFDATA_DIR}'"
echo "    input         = '${SAMPLESHEET}'"
echo "  }"
echo ""
echo "Then run:"
echo "  nextflow run . -profile test,docker -stub-run --outdir /tmp/sash-stub-test"
