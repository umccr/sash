import Utils

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSash.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.ref_data_path,
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// TODO(SW): place this into parameter validation
// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BOLT_OTHER_CANCER_REPORT   } from '../modules/local/bolt/other/cancer_report/main'
include { BOLT_OTHER_MULTIQC_REPORT  } from '../modules/local/bolt/other/multiqc_report/main'
include { BOLT_OTHER_PURPLE_BAF_PLOT } from '../modules/local/bolt/other/purple_baf_plot/main'
include { BOLT_SMLV_GERMLINE_PREPARE } from '../modules/local/bolt/smlv_germline/prepare/main'
include { BOLT_SMLV_GERMLINE_REPORT  } from '../modules/local/bolt/smlv_germline/report/main'
include { BOLT_SMLV_SOMATIC_ANNOTATE } from '../modules/local/bolt/smlv_somatic/annotate/main'
include { BOLT_SMLV_SOMATIC_FILTER   } from '../modules/local/bolt/smlv_somatic/filter/main'
include { BOLT_SMLV_SOMATIC_REPORT   } from '../modules/local/bolt/smlv_somatic/report/main'
include { BOLT_SMLV_SOMATIC_RESCUE   } from '../modules/local/bolt/smlv_somatic/rescue/main'
include { BOLT_SV_SOMATIC_ANNOTATE   } from '../modules/local/bolt/sv_somatic/annotate/main'
include { BOLT_SV_SOMATIC_PRIORITISE } from '../modules/local/bolt/sv_somatic/prioritise/main'

include { LINX_ANNOTATION            } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING              } from '../subworkflows/local/linx_plotting'
include { PREPARE_INPUT              } from '../subworkflows/local/prepare_input'
include { PREPARE_REFERENCE          } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING             } from '../subworkflows/local/purple_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SASH {
    // Create channel for versions
    // channel: [versions.yml]
    ch_versions = Channel.empty()




    //
    // SUBWORKFLOW: Prepare inputs from samplesheet
    //

    // channel: [val(meta)]
    PREPARE_INPUT(
        file(params.input),
    )
    ch_inputs = PREPARE_INPUT.out.metas

    // TODO(SW): place into subworkflow; simplify with map of variable name -> subpath
    ch_amber = ch_inputs.map { meta -> [meta, file(meta.oncoanalyser_dir).toUriString() + '/amber/'] }
    ch_cobalt = ch_inputs.map { meta -> [meta, file(meta.oncoanalyser_dir).toUriString() + '/cobalt/'] }

    ch_gripss_somatic = ch_inputs
        .map { meta ->
            def subpath = "/gripss/somatic/${meta.tumor_id}.gripss.filtered.somatic.vcf.gz"
            def vcf = file(meta.oncoanalyser_dir).toUriString() + subpath
            return [meta, vcf, "${vcf}.tbi"]
        }

    ch_gripss_somatic_unfiltered = ch_inputs
        .map { meta ->
            def subpath = "/gripss/somatic/${meta.tumor_id}.gripss.somatic.vcf.gz"
            def vcf = file(meta.oncoanalyser_dir).toUriString() + subpath
            return [meta, vcf, "${vcf}.tbi"]
        }

    ch_gripss_germline = ch_inputs
        .map { meta ->
            def subpath = "/gripss/germline/${meta.tumor_id}.gripss.filtered.germline.vcf.gz"
            def vcf = file(meta.oncoanalyser_dir).toUriString() + subpath
            return [meta, vcf, "${vcf}.tbi"]
        }

    ch_sage_somatic = ch_inputs
        .map { meta ->
            def subpath = "/sage/somatic/${meta.tumor_id}.sage.somatic.filtered.vcf.gz"
            def vcf = file(meta.oncoanalyser_dir).toUriString() + subpath
            return [meta, vcf, "${vcf}.tbi"]
        }

    ch_virusbreakend = ch_inputs
        .map { meta ->
            def subpath = '/virusbreakend/'
            def virusbreakend_dir = file(meta.oncoanalyser_dir).toUriString() + subpath
            return [meta, virusbreakend_dir]
        }




    //
    // SUBWORKFLOW: Prepare reference data
    //

    // channel: [val(meta)]
    PREPARE_REFERENCE()
    genome = PREPARE_REFERENCE.out.genome
    umccr_data = PREPARE_REFERENCE.out.umccr_data
    hmf_data = PREPARE_REFERENCE.out.hmf_data




    //
    // MODULE: Somatic small variants
    //

    // channel: [meta_bolt, dragen_vcf, dragen_tbi]
    ch_input_vcf_somatic = ch_inputs
        .map { meta ->
            def vcf = file(meta.dragen_somatic_dir).toUriString() + "/${meta.tumor_id}.hard-filtered.vcf.gz"
            return [meta, vcf, "${vcf}.tbi"]
        }


    // channel: [meta_bolt, dragen_vcf, dragen_tbi, sage_vcf, sage_tbi]
    ch_smlv_somatic_rescue_inputs = WorkflowSash.groupByMeta(
        ch_input_vcf_somatic,
        ch_sage_somatic,
    )
        .map {
            def meta = it[0]
            def meta_bolt = [
                key: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
            ]
            return [meta_bolt, *it[1..-1]]
        }

    BOLT_SMLV_SOMATIC_RESCUE(
        ch_smlv_somatic_rescue_inputs,
        umccr_data.hotspots,
    )

    BOLT_SMLV_SOMATIC_ANNOTATE(
        BOLT_SMLV_SOMATIC_RESCUE.out.vcf,
        umccr_data.somatic_driver_panel_regions_gene,
        umccr_data.annotations_dir,
        umccr_data.pon_dir,
        umccr_data.pcgr_dir,
    )

    BOLT_SMLV_SOMATIC_FILTER(
        BOLT_SMLV_SOMATIC_ANNOTATE.out.vcf,
    )

    // Set outputs, restoring original meta
    ch_versions = ch_versions.mix(
        BOLT_SMLV_SOMATIC_RESCUE.out.versions,
        BOLT_SMLV_SOMATIC_ANNOTATE.out.versions,
        BOLT_SMLV_SOMATIC_FILTER.out.versions,
    )

    // channel: [meta, smlv_somatic_vcf]
    ch_smlv_somatic_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_FILTER.out.vcf, ch_inputs)
    ch_smlv_somatic_unfiltered_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_FILTER.out.vcf_unfiltered, ch_inputs)




    //
    // Germline small variants
    //

    // channel: [meta_bolt, dragen_vcf]
    ch_input_vcf_germline = ch_inputs
        .map { meta ->
            def vcf = file(meta.dragen_germline_dir).toUriString() + "/${meta.normal_id}.hard-filtered.vcf.gz"
            return [meta, vcf]
        }

    // channel: [meta_bolt, dragen_vcf]
    ch_smlv_germline_prepare_inputs = ch_input_vcf_germline
        .map { meta, dragen_vcf ->

            def meta_bolt = [
                key: meta.id,
                normal_id: meta.normal_id,
            ]

            return [meta_bolt, dragen_vcf]
        }

    BOLT_SMLV_GERMLINE_PREPARE(
        ch_smlv_germline_prepare_inputs,
        umccr_data.germline_predisposition_panel_regions_transcript,
    )

    // Set outputs, restoring original meta
    ch_versions = ch_versions.mix(
        BOLT_SMLV_GERMLINE_PREPARE.out.versions,
    )

    // channel: [meta, smlv_germline_vcf]
    ch_smlv_germline_out = WorkflowSash.restoreMeta(BOLT_SMLV_GERMLINE_PREPARE.out.vcf, ch_inputs)




    //
    // CNV calling using UMCCR postprocessed variants
    //

    PURPLE_CALLING(
        ch_inputs,
        ch_amber,
        ch_cobalt,
        ch_smlv_somatic_out,


        // NOTE(SW): PURPLE germline enrichment requires tumor AD and GT information in the germline calls but DRAGEN does not generate such a file, see:
        //   * https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlineVariantEnrichment.java#L30
        //   * https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlinePurityEnrichment.java#L50
        //   * https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlineGenotypeEnrichment.java#L63
        //ch_smlv_germline_out,
        ch_smlv_germline_out.map { meta, vcf -> return [meta, []] },


        ch_gripss_somatic,
        ch_gripss_germline,
        ch_gripss_somatic_unfiltered,
        genome.fasta,
        genome.version,
        genome.fai,
        genome.dict,
        hmf_data.gc_profile,
        hmf_data.sage_known_hotspots_somatic,
        hmf_data.sage_known_hotspots_germline,
        hmf_data.driver_gene_panel,
        hmf_data.ensembl_data_resources,
        hmf_data.purple_germline_del,
    )




    //
    // Small variant reporting (PCGR, CPSR, stats)
    //

    ch_smlv_somatic_report_inputs = WorkflowSash.groupByMeta(
        ch_smlv_somatic_out,
        ch_input_vcf_somatic,
        PURPLE_CALLING.out.purple_dir,
    )
        .map { meta, vcf, vcf_unfiltered, tbi_unfiltered, purple_dir ->
            def purity_tsv = file(purple_dir).resolve("${meta.tumor_id}.purple.purity.tsv")

            // Require small variant VCF files and purity file from the PURPLE directory
            if (!purity_tsv.exists()) {


                // TODO(SW): consider how to replace this, not actually handled in groupByMeta right now


                //return Constants.META_PLACEHOLDER
                return 'Constants.META_PLACEHOLDER'
            }

            def meta_bolt = [
                key: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
            ]

            return [meta_bolt, vcf, vcf_unfiltered, purity_tsv]
        }

    BOLT_SMLV_SOMATIC_REPORT(
        ch_smlv_somatic_report_inputs,
        umccr_data.pcgr_dir,
        umccr_data.somatic_driver_panel_regions_coding,
        hmf_data.sage_highconf_regions,
        genome.fasta,
    )

    ch_smlv_germline_report_inputs = WorkflowSash.groupByMeta(
        ch_smlv_germline_out,
        ch_input_vcf_germline,
    )
        .map { meta, vcf, vcf_unfiltered ->
            def meta_bolt = [
                key: meta.id,
                id: meta.normal_id,
            ]
            return [meta_bolt, vcf, vcf_unfiltered]
        }

    BOLT_SMLV_GERMLINE_REPORT(
        ch_smlv_germline_report_inputs,
        umccr_data.germline_predisposition_panel_genes,
        umccr_data.pcgr_dir,
    )

    ch_versions = ch_versions.mix(
        BOLT_SMLV_SOMATIC_REPORT.out.versions,
        BOLT_SMLV_GERMLINE_REPORT.out.versions,
    )

    ch_smlv_somatic_report_af_global_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.af_global, ch_inputs)
    ch_smlv_somatic_report_af_keygenes_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.af_keygenes, ch_inputs)
    ch_smlv_somatic_report_stats_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.bcftools_stats, ch_inputs)
    ch_smlv_somatic_report_counts_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.variant_counts, ch_inputs)

    ch_smlv_germline_report_stats_out = WorkflowSash.restoreMeta(BOLT_SMLV_GERMLINE_REPORT.out.bcftools_stats, ch_inputs)
    ch_smlv_germline_report_counts_out = WorkflowSash.restoreMeta(BOLT_SMLV_GERMLINE_REPORT.out.variant_counts, ch_inputs)




    //
    // Somatic structural variants
    //

    // channel: [meta_bolt, sv_vcf, cnv_tsv]
    ch_sv_somatic_inputs = PURPLE_CALLING.out.purple_dir
        .map { meta, purple_dir ->

            def sv_vcf = file(purple_dir).resolve("${meta.tumor_id}.purple.sv.vcf.gz")
            def cnv_tsv = file(purple_dir).resolve("${meta.tumor_id}.purple.cnv.somatic.tsv")

            // Require SV VCF and CNV TSV from the PURPLE directory
            if (!sv_vcf.exists() || !cnv_tsv.exists()) {


                // TODO(SW): consider how to replace this, not actually handled in groupByMeta right now


                //return Constants.META_PLACEHOLDER
                return 'Constants.META_PLACEHOLDER'
            }

            def meta_bolt = [
                key: meta.id,
                id: meta.tumor_id,
            ]
            return [meta_bolt, sv_vcf, cnv_tsv]
        }

    BOLT_SV_SOMATIC_ANNOTATE(
        ch_sv_somatic_inputs,
        genome.fasta,
        genome.fai,
        umccr_data.snpeff_dir,
    )

    BOLT_SV_SOMATIC_PRIORITISE(
        BOLT_SV_SOMATIC_ANNOTATE.out.vcf,
        umccr_data.known_fusion_pairs,
        umccr_data.known_fusion_heads,
        umccr_data.known_fusion_tails,
        umccr_data.known_fusioncatcher_pairs,
        umccr_data.somatic_driver_panel_genes,
        umccr_data.somatic_driver_panel_genes_ts,
        umccr_data.appris,
    )

    // Set outputs, restoring original meta
    ch_versions = ch_versions.mix(
        BOLT_SV_SOMATIC_ANNOTATE.out.versions,
        BOLT_SV_SOMATIC_PRIORITISE.out.versions,
    )

    ch_sv_somatic_sv_tsv_out = WorkflowSash.restoreMeta(BOLT_SV_SOMATIC_PRIORITISE.out.sv_tsv, ch_inputs)
    ch_sv_somatic_sv_vcf_out = WorkflowSash.restoreMeta(BOLT_SV_SOMATIC_PRIORITISE.out.sv_vcf, ch_inputs)
    ch_sv_somatic_cnv_tsv_out = WorkflowSash.restoreMeta(BOLT_SV_SOMATIC_PRIORITISE.out.cnv_tsv, ch_inputs)




    //
    // Generate custom PURPLE β-allele frequency circos plot
    //

    // channel: [meta_bolt, purple_dir]
    ch_other_purple_baf_plot_inputs = PURPLE_CALLING.out.purple_dir
        .map { meta, purple_dir ->
            def meta_bolt = [
                key: meta.id,
                tumor_id: meta.tumor_id,
            ]
            return [meta_bolt, purple_dir]
        }

    BOLT_OTHER_PURPLE_BAF_PLOT(
        ch_other_purple_baf_plot_inputs,
        umccr_data.purple_baf_circos_config,
        umccr_data.purple_baf_circos_gaps,
    )

    ch_versions = ch_versions.mix(BOLT_OTHER_PURPLE_BAF_PLOT.out.versions)

    ch_purple_baf_plot_out = WorkflowSash.restoreMeta(BOLT_OTHER_PURPLE_BAF_PLOT.out.plot, ch_inputs)




    //
    // Generate the cancer report
    //

    // channel: [ meta, smlv_somatic_vcf, sv_somatic_tsv, sv_somatic_vcf, cnv_somatic_tsv, af_global, af_keygenes, purple_baf_plot, purple_dir, virusbreakend_dir ]
    ch_cancer_report_inputs = WorkflowSash.groupByMeta(
        ch_smlv_somatic_out,
        ch_sv_somatic_sv_tsv_out,
        ch_sv_somatic_sv_vcf_out,
        ch_sv_somatic_cnv_tsv_out,
        ch_smlv_somatic_report_af_global_out,
        ch_smlv_somatic_report_af_keygenes_out,
        ch_purple_baf_plot_out,
        PURPLE_CALLING.out.purple_dir,
        ch_virusbreakend,
    )
        .map {
            def meta = it[0]
            def meta_bolt = [
                key: meta.id,
                subject_id: meta.subject_id,
                tumor_id: meta.tumor_id,
            ]
            return [meta_bolt, *it[1..-1]]
        }

    BOLT_OTHER_CANCER_REPORT(
        ch_cancer_report_inputs,
        umccr_data.somatic_driver_panel,
    )




    //
    // Generate MultiQC report
    //

    ch_input_dragen_somatic_dir = ch_inputs
        .map { meta -> [meta, meta.dragen_somatic_dir] }

    ch_input_dragen_germline_dir = ch_inputs
        .map { meta -> [meta, meta.dragen_germline_dir] }

    // channel: [ meta, [somatic_dragen_dir, germline_dragen_dir, somatic_bcftools_stats, germline_bcftools_stats, somatic_variant_counts, germline_variant_counts, purple_dir] ]
    ch_multiqc_report_inputs = WorkflowSash.groupByMeta(
        ch_input_dragen_somatic_dir,
        ch_input_dragen_germline_dir,
        ch_smlv_somatic_report_stats_out,
        ch_smlv_germline_report_stats_out,
        ch_smlv_somatic_report_counts_out,
        ch_smlv_germline_report_counts_out,
        PURPLE_CALLING.out.purple_dir,
    )
        .map {
            def meta = it[0]
            def other = it[1..-1]

            def meta_multiqc = [
                key: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
            ]

            return [meta_multiqc, other]
        }

    BOLT_OTHER_MULTIQC_REPORT(
        ch_multiqc_report_inputs,
    )




    //
    // Annotate post processed strucutral variant events
    //

    LINX_ANNOTATION(
        ch_inputs,
        PURPLE_CALLING.out.purple_dir,
        genome.version,
        hmf_data.ensembl_data_resources,
        hmf_data.known_fusion_data,
        hmf_data.driver_gene_panel,
    )
    ch_versions = ch_versions.mix(LINX_ANNOTATION.out.versions)

    LINX_PLOTTING(
        ch_inputs,
        LINX_ANNOTATION.out.somatic,
        genome.version,
        hmf_data.ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(LINX_PLOTTING.out.versions)
    ch_linx_somatic_plots_out = WorkflowSash.restoreMeta(LINX_PLOTTING.out.visualiser_dir, ch_inputs)




    //
    // MODULE: Collect software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
