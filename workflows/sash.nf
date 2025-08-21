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
include { ESVEE_CALL                 } from '../modules/local/esvee/call/main'
include { PAVE_SOMATIC               } from '../modules/local/pave/somatic/main'
include { VCF2MAF                    } from '../modules/local/vcf2maf/main'

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

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SASH {
    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()




    //
    // Prepare inputs from samplesheet
    //

    PREPARE_INPUT(
        file(params.input),
    )

    ch_inputs              = PREPARE_INPUT.out.metas             // channel: [ meta ]

    // OncoAnalyser inputs
    ch_amber               = PREPARE_INPUT.out.amber             // channel: [ meta, amber_dir ]
    ch_cobalt              = PREPARE_INPUT.out.cobalt            // channel: [ meta, cobalt_dir ]
    ch_sage_somatic        = PREPARE_INPUT.out.sage_somatic      // channel: [ meta, sage_somatic_vcf, sage_somatic_tbi ]
    ch_virusbreakend       = PREPARE_INPUT.out.virusbreakend     // channel: [ meta, virusbreakend_dir ]
    ch_call_inputs         = PREPARE_INPUT.out.call_inputs       // channel: [ meta_esvee, esvee_ref_depth_vcf, esvee_prep_dir ]

    // DRAGEN inputs
    ch_input_hrd           = PREPARE_INPUT.out.hrd               // channel: [ meta, hrdscore_csv ]
    ch_input_vcf_germline  = PREPARE_INPUT.out.vcf_germline      // channel: [ meta, dragen_germline_vcf ]
    ch_input_vcf_somatic   = PREPARE_INPUT.out.vcf_somatic       // channel: [ meta, dragen_somatic_vcf, dragen_somatic_tbi ]




    //
    // Prepare reference data
    //

    // channel: [ meta ]
    PREPARE_REFERENCE()
    genome               = PREPARE_REFERENCE.out.genome
    umccr_data           = PREPARE_REFERENCE.out.umccr_data
    hmf_data             = PREPARE_REFERENCE.out.hmf_data
    misc_data            = PREPARE_REFERENCE.out.misc_data




    //
    // Somatic small variants
    //



    // channel: [ meta_bolt, dragen_somatic_vcf, dragen_somatic_tbi, sage_somatic_vcf, sage_somatic_tbi ]
    ch_smlv_somatic_rescue_inputs = WorkflowSash.groupByMeta(
        ch_input_vcf_somatic,
        ch_sage_somatic,
    )
        .map {
            def meta = it[0]
            def meta_bolt = [
                key: meta.id,
                id: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
                sample_id: meta.tumor_id
            ]
            return [meta_bolt, *it[1..-1]]
        }

    BOLT_SMLV_SOMATIC_RESCUE(
        ch_smlv_somatic_rescue_inputs,
        umccr_data.hotspots,
    )

    ch_versions = ch_versions.mix(BOLT_SMLV_SOMATIC_RESCUE.out.versions)

    BOLT_SMLV_SOMATIC_ANNOTATE(
        BOLT_SMLV_SOMATIC_RESCUE.out.vcf,
        umccr_data.somatic_panel_regions_gene,
        umccr_data.annotations_dir,
        misc_data.pon_dir,
        misc_data.pcgr_dir,
    )

    ch_versions = ch_versions.mix(BOLT_SMLV_SOMATIC_ANNOTATE.out.versions)

    BOLT_SMLV_SOMATIC_FILTER(
        BOLT_SMLV_SOMATIC_ANNOTATE.out.vcf,
    )

    ch_versions = ch_versions.mix(BOLT_SMLV_SOMATIC_FILTER.out.versions)

    // channel: [ meta, smlv_somatic_vcf ]
    ch_smlv_somatic_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_FILTER.out.vcf, ch_inputs)
        .map { meta, vcf, tbi -> [meta, vcf] }
    // channel: [ meta, smlv_somatic_filters_vcf ]
    ch_smlv_somatic_filters_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_FILTER.out.vcf_filters, ch_inputs)

    // NOTE(SW): PAVE is run so that we obtain a complete PURPLE driver catalogue
    PAVE_SOMATIC(
        BOLT_SMLV_SOMATIC_FILTER.out.vcf,
        genome.fasta,
        genome.version,
        genome.fai,
        hmf_data.clinvar_annotations,
        hmf_data.segment_mappability,
        umccr_data.driver_gene_panel,
        umccr_data.ensembl_data_resources,
    )

    // channel: [ meta, pave_somatic_vcf ]
    ch_pave_somatic_out = WorkflowSash.restoreMeta(PAVE_SOMATIC.out.vcf, ch_inputs)




    //
    // Convert somatic VCF to MAF format
    //

    VCF2MAF(
        ch_smlv_somatic_out,
        genome.fasta
    )

    ch_versions = ch_versions.mix(VCF2MAF.out.versions)

    // channel: [ meta, somatic_maf ]
    ch_vcf2maf_out = VCF2MAF.out.maf




    //
    // Germline small variants
    //

    // channel: [ meta_bolt, dragen_germline_vcf ]
    ch_smlv_germline_prepare_inputs = ch_input_vcf_germline
        .map { meta, dragen_vcf ->

            def meta_bolt = [
                key: meta.id,
                id: meta.id,
                normal_id: meta.normal_id,
            ]

            return [meta_bolt, dragen_vcf]
        }

    BOLT_SMLV_GERMLINE_PREPARE(
        ch_smlv_germline_prepare_inputs,
        umccr_data.germline_panel_regions_transcript,
    )

    ch_versions = ch_versions.mix(BOLT_SMLV_GERMLINE_PREPARE.out.versions)

    // channel: [ meta, smlv_germline_vcf ]
    ch_smlv_germline_out = WorkflowSash.restoreMeta(BOLT_SMLV_GERMLINE_PREPARE.out.vcf, ch_inputs)




    //
    // Somatic structural variants

    //

    ESVEE_CALL(
        ch_call_inputs,
        genome.fasta,
        genome.version,
        hmf_data.gridss_pon_breakends,
        hmf_data.gridss_pon_breakpoints,
        umccr_data.known_fusions,
        hmf_data.repeatmasker_annotations,
    )

    ch_versions = ch_versions.mix(ESVEE_CALL.out.versions)

    ch_esvee_somatic_out = WorkflowSash.restoreMeta(ESVEE_CALL.out.somatic_vcf, ch_inputs)
    ch_esvee_germline_out =  WorkflowSash.restoreMeta(ESVEE_CALL.out.germline_vcf, ch_inputs)




    //
    // CNV calling using UMCCR postprocessed variants
    //

    PURPLE_CALLING(
        ch_inputs,
        ch_amber,
        ch_cobalt,
        ch_pave_somatic_out,


        // NOTE(SW): PURPLE germline enrichment requires tumor AD and GT information in the germline calls but DRAGEN does not generate such a file, see:
        //   * https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlineVariantEnrichment.java#L30
        //   * https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlinePurityEnrichment.java#L50
        //   * https://github.com/hartwigmedical/hmftools/blob/a2f82e5/purple/src/main/java/com/hartwig/hmftools/purple/germline/GermlineGenotypeEnrichment.java#L63
        //ch_smlv_germline_out,
        ch_smlv_germline_out.map { meta, vcf -> return [meta, []] },

        ch_esvee_somatic_out,
        ch_esvee_germline_out,
        genome.fasta,
        genome.version,
        genome.fai,
        genome.dict,
        hmf_data.gc_profile,
        hmf_data.sage_known_hotspots_somatic,
        umccr_data.sage_known_hotspots_germline,
        umccr_data.driver_gene_panel,
        umccr_data.ensembl_data_resources,
        hmf_data.purple_germline_del,
    )

    ch_versions = ch_versions.mix(PURPLE_CALLING.out.versions)




    //
    // Small variant reporting (PCGR, CPSR, stats)
    //

    // channel: [ meta_bolt, smlv_somatic_vcf, smlv_somatic_filters_vcf, dragen_somatic_vcf, dragen_somatic_tbi, purple_dir ]
    ch_smlv_somatic_report_inputs = WorkflowSash.groupByMeta(
        ch_smlv_somatic_out,
        ch_smlv_somatic_filters_out,
        ch_input_vcf_somatic,
        PURPLE_CALLING.out.purple_dir,
    )
        .map { meta, vcf, vcf_filters, vcf_dragen, tbi_dragen, purple_dir ->
            def purity_tsv = file(purple_dir).resolve("${meta.tumor_id}.purple.purity.tsv")

            def meta_bolt = [
                key: meta.id,
                id: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
            ]

            return [meta_bolt, vcf, vcf_filters, vcf_dragen, purity_tsv]
        }

    BOLT_SMLV_SOMATIC_REPORT(
        ch_smlv_somatic_report_inputs,
        misc_data.pcgr_dir,
        umccr_data.somatic_panel_regions_cds,
        hmf_data.sage_highconf_regions,
        genome.fasta,
        genome.fai,
    )

    ch_versions = ch_versions.mix(BOLT_SMLV_SOMATIC_REPORT.out.versions)

    // channel: [ meta_bolt, smlv_germline_vcf, dragen_germline_vcf ]
    ch_smlv_germline_report_inputs = WorkflowSash.groupByMeta(
        ch_smlv_germline_out,
        ch_input_vcf_germline,
    )
        .map { meta, vcf, vcf_dragen ->
            def meta_bolt = [
                key: meta.id,
                id: meta.id,
                normal_id: meta.normal_id,
            ]
            return [meta_bolt, vcf, vcf_dragen]
        }

    BOLT_SMLV_GERMLINE_REPORT(
        ch_smlv_germline_report_inputs,
        umccr_data.germline_panel_genes,
        misc_data.pcgr_dir,
    )

    ch_versions = ch_versions.mix(BOLT_SMLV_GERMLINE_REPORT.out.versions)

    // channel: [ meta, af_global ]
    ch_smlv_somatic_report_af_global_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.af_global, ch_inputs)
    // channel: [ meta, af_keygenes ]
    ch_smlv_somatic_report_af_keygenes_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.af_keygenes, ch_inputs)
    // channel: [ meta, bcftools_stats_somatic ]
    ch_smlv_somatic_report_stats_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.bcftools_stats, ch_inputs)
    // channel: [ meta, variant_counts_type_somatic ]
    ch_smlv_somatic_report_counts_type_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.counts_type, ch_inputs)
    // channel: [ meta, variant_counts_process_somatic ]
    ch_smlv_somatic_report_counts_process_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_REPORT.out.counts_process, ch_inputs)

    // channel: [ meta, bcftools_stats_germline ]
    ch_smlv_germline_report_stats_out = WorkflowSash.restoreMeta(BOLT_SMLV_GERMLINE_REPORT.out.bcftools_stats, ch_inputs)
    // channel: [ meta, variant_counts_type_germline ]
    ch_smlv_germline_report_counts_type_out = WorkflowSash.restoreMeta(BOLT_SMLV_GERMLINE_REPORT.out.counts_type, ch_inputs)




    //
    // Somatic structural variants
    //

    // channel: [ meta_bolt, sv_vcf, cnv_tsv ]
    ch_sv_somatic_inputs = PURPLE_CALLING.out.purple_dir
        .map { meta, purple_dir ->
            def sv_vcf = file(purple_dir).resolve("${meta.tumor_id}.purple.sv.vcf.gz")
            def cnv_tsv = file(purple_dir).resolve("${meta.tumor_id}.purple.cnv.somatic.tsv")

            def meta_bolt = [
                key: meta.id,
                id: meta.id,
                tumor_id: meta.tumor_id,
            ]
            return [meta_bolt, sv_vcf, cnv_tsv]
        }

    BOLT_SV_SOMATIC_ANNOTATE(
        ch_sv_somatic_inputs,
        genome.fasta,
        genome.fai,
        misc_data.snpeff_dir,
    )

    ch_versions = ch_versions.mix(BOLT_SV_SOMATIC_ANNOTATE.out.versions)

    BOLT_SV_SOMATIC_PRIORITISE(
        BOLT_SV_SOMATIC_ANNOTATE.out.vcf,
        umccr_data.known_fusion_pairs,
        umccr_data.known_fusion_five,
        umccr_data.known_fusion_three,
        umccr_data.somatic_panel_genes,
        umccr_data.somatic_panel_genes_ts,
        misc_data.appris,
    )

    ch_versions = ch_versions.mix(BOLT_SV_SOMATIC_PRIORITISE.out.versions)

    // channel: [ meta, sv_tsv ]
    ch_sv_somatic_sv_tsv_out = WorkflowSash.restoreMeta(BOLT_SV_SOMATIC_PRIORITISE.out.sv_tsv, ch_inputs)
    // channel: [ meta, sv_vcf ]
    ch_sv_somatic_sv_vcf_out = WorkflowSash.restoreMeta(BOLT_SV_SOMATIC_PRIORITISE.out.sv_vcf, ch_inputs)
    // channel: [ meta, cnv_tsv ]
    ch_sv_somatic_cnv_tsv_out = WorkflowSash.restoreMeta(BOLT_SV_SOMATIC_PRIORITISE.out.cnv_tsv, ch_inputs)



    //
    // Generate custom PURPLE β-allele frequency circos plot
    //

    // channel: [ meta_bolt, purple_dir ]
    ch_other_purple_baf_plot_inputs = PURPLE_CALLING.out.purple_dir
        .map { meta, purple_dir ->
            def meta_bolt = [
                key: meta.id,
                id: meta.id,
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

    // channel: [ meta, purple_baf_circos_plot ]
    ch_purple_baf_plot_out = WorkflowSash.restoreMeta(BOLT_OTHER_PURPLE_BAF_PLOT.out.plot, ch_inputs)


    // VCF2MAF
    VCF2MAF(
        ch_smlv_somatic_out,
        umccr_data.vep_cache,
        ref

    )


    //
    // Generate the cancer report
    //

    // channel: [ meta_bolt, smlv_somatic_vcf, smlv_somatic_bcftools_stats, smlv_somatic_counts_process, sv_tsv, sv_vcf, cnv_tsv, af_global, af_keygenes, purple_baf_circos_plot, purple_dir, virusbreakend_dir, dragen_hrd ]
    ch_cancer_report_inputs = WorkflowSash.groupByMeta(
        ch_smlv_somatic_out,
        ch_smlv_somatic_report_stats_out,
        ch_smlv_somatic_report_counts_process_out,
        ch_sv_somatic_sv_tsv_out,
        ch_sv_somatic_sv_vcf_out,
        ch_sv_somatic_cnv_tsv_out,
        ch_smlv_somatic_report_af_global_out,
        ch_smlv_somatic_report_af_keygenes_out,
        ch_purple_baf_plot_out,
        PURPLE_CALLING.out.purple_dir,
        ch_virusbreakend,
        ch_input_hrd,
    )
        .map {
            def meta = it[0]
            def meta_bolt = [
                key: meta.id,
                id: meta.id,
                subject_id: meta.subject_id,
                tumor_id: meta.tumor_id,
            ]
            return [meta_bolt, *it[1..-1]]
        }

    BOLT_OTHER_CANCER_REPORT(
        ch_cancer_report_inputs,
        umccr_data.somatic_panel,
        misc_data.oncokb_genes,
    )

    ch_versions = ch_versions.mix(BOLT_OTHER_CANCER_REPORT.out.versions)




    //
    // Generate MultiQC report
    //

    // channel: [ meta, somatic_dragen_dir ]
    ch_input_dragen_somatic_dir = ch_inputs
        .map { meta -> [meta, meta.dragen_somatic_dir] }

    // channel: [ meta, germline_dragen_dir ]
    ch_input_dragen_germline_dir = ch_inputs
        .map { meta -> [meta, meta.dragen_germline_dir] }

    // channel: [ meta_multiqc, [somatic_dragen_dir, germline_dragen_dir, somatic_bcftools_stats, germline_bcftools_stats, somatic_counts_type, germline_counts_type, purple_dir] ]
    ch_multiqc_report_inputs = WorkflowSash.groupByMeta(
        ch_input_dragen_somatic_dir,
        ch_input_dragen_germline_dir,
        ch_smlv_somatic_report_stats_out,
        ch_smlv_germline_report_stats_out,
        ch_smlv_somatic_report_counts_type_out,
        ch_smlv_germline_report_counts_type_out,
        PURPLE_CALLING.out.purple_dir,
    )
        .map {
            def meta = it[0]
            def other = it[1..-1]

            def meta_multiqc = [
                key: meta.id,
                id: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
            ]

            return [meta_multiqc, other]
        }

    BOLT_OTHER_MULTIQC_REPORT(
        ch_multiqc_report_inputs,
    )

    ch_versions = ch_versions.mix(BOLT_OTHER_MULTIQC_REPORT.out.versions)




    //
    // Annotate post processed strucutral variant events
    //

    LINX_ANNOTATION(
        ch_inputs,
        PURPLE_CALLING.out.purple_dir,
        genome.version,
        umccr_data.ensembl_data_resources,
        umccr_data.known_fusion_data,
        umccr_data.driver_gene_panel,
    )
    ch_versions = ch_versions.mix(LINX_ANNOTATION.out.versions)

    LINX_PLOTTING(
        ch_inputs,
        LINX_ANNOTATION.out.somatic,
        genome.version,
        umccr_data.ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(LINX_PLOTTING.out.versions)




    //
    // TASK: Aggregate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true,
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
