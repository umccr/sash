#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    umccr/sash
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/umccr/sash----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SASH } from './workflows/sash'

//
// WORKFLOW: Run main umccr/sash analysis pipeline
//
workflow UMCCR_SASH {
    main:
    SASH()

    emit:
    // NOTE(QC): pass-through of SASH's publishable outputs so the entry `workflow {}`
    // block below can wire them into its `publish:` section / top-level `output {}` block.
    smlv_somatic_rescue_vcf         = SASH.out.smlv_somatic_rescue_vcf
    smlv_somatic_annotate_vcf       = SASH.out.smlv_somatic_annotate_vcf
    smlv_somatic_filter_vcf         = SASH.out.smlv_somatic_filter_vcf
    smlv_somatic_filter_vcf_filters = SASH.out.smlv_somatic_filter_vcf_filters
    smlv_germline_prepare_vcf       = SASH.out.smlv_germline_prepare_vcf
    sv_somatic_annotate_vcf         = SASH.out.sv_somatic_annotate_vcf
    sv_somatic_prioritise_sv_tsv    = SASH.out.sv_somatic_prioritise_sv_tsv
    sv_somatic_prioritise_sv_vcf    = SASH.out.sv_somatic_prioritise_sv_vcf
    sv_somatic_prioritise_cnv_tsv   = SASH.out.sv_somatic_prioritise_cnv_tsv
    purple_dir                      = SASH.out.purple_dir
    sigrap_hrdetect_json            = SASH.out.sigrap_hrdetect_json
    sigrap_mutpat_dir               = SASH.out.sigrap_mutpat_dir
    esvee_caller_dir                = SASH.out.esvee_caller_dir
    linx_germline_annotation_dir    = SASH.out.linx_germline_annotation_dir
    linx_somatic_annotation_dir     = SASH.out.linx_somatic_annotation_dir
    linx_somatic_plot_dir           = SASH.out.linx_somatic_plot_dir
    linxreport_html                 = SASH.out.linxreport_html
    vcf2maf_maf                     = SASH.out.vcf2maf_maf
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    main:
    // Validate parameters and print summary to screen
    WorkflowMain.initialise(workflow, params, log)

    UMCCR_SASH()

    // Completion email and summary
    workflow.onComplete {
        def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
        if (params.email || params.email_on_fail) {
            NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
        }
        NfcoreTemplate.summary(workflow, params, log)
        if (params.hook_url) {
            NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
        }
    }

    // NOTE(QC): the workflow output `path {}` closures below can only emit `>>` publish
    // statements from their own (non-nested) closure body, so directory outputs that need
    // per-file flattening (to reproduce the old publishDir/saveAs behaviour exactly) are
    // expanded into one channel item per file here, ahead of the `publish:` section.
    ch_linx_germline_annotation_files = UMCCR_SASH.out.linx_germline_annotation_dir
        .flatMap { meta, dir -> file(dir).listFiles().collect { f -> [meta, f] } }
    ch_linx_somatic_annotation_files = UMCCR_SASH.out.linx_somatic_annotation_dir
        .flatMap { meta, dir -> file(dir).listFiles().collect { f -> [meta, f] } }
    ch_linx_somatic_plot_files = UMCCR_SASH.out.linx_somatic_plot_dir
        .flatMap { meta, dir ->
            def base = file(dir)
            def files = []
            base.eachFileRecurse(groovy.io.FileType.FILES) { f -> files << f }
            files.collect { f -> [meta, f, base.relativize(f)] }
        }

    publish:
    // NOTE(QC): pilot migration of publishDir -> workflow output definition (see conf/modules.config
    // and the `output {}` block below) for the 15 "simple" processes; see git history for details.
    smlv_somatic_rescue_vcf         = UMCCR_SASH.out.smlv_somatic_rescue_vcf
    smlv_somatic_annotate_vcf       = UMCCR_SASH.out.smlv_somatic_annotate_vcf
    smlv_somatic_filter_vcf         = UMCCR_SASH.out.smlv_somatic_filter_vcf
    smlv_somatic_filter_vcf_filters = UMCCR_SASH.out.smlv_somatic_filter_vcf_filters
    smlv_germline_prepare_vcf       = UMCCR_SASH.out.smlv_germline_prepare_vcf
    sv_somatic_annotate_vcf         = UMCCR_SASH.out.sv_somatic_annotate_vcf
    sv_somatic_prioritise_sv_tsv    = UMCCR_SASH.out.sv_somatic_prioritise_sv_tsv
    sv_somatic_prioritise_sv_vcf    = UMCCR_SASH.out.sv_somatic_prioritise_sv_vcf
    sv_somatic_prioritise_cnv_tsv   = UMCCR_SASH.out.sv_somatic_prioritise_cnv_tsv
    purple_dir                      = UMCCR_SASH.out.purple_dir
    sigrap_hrdetect_json            = UMCCR_SASH.out.sigrap_hrdetect_json
    sigrap_mutpat_dir               = UMCCR_SASH.out.sigrap_mutpat_dir
    esvee_caller_dir                = UMCCR_SASH.out.esvee_caller_dir
    linx_germline_annotation_dir    = ch_linx_germline_annotation_files
    linx_somatic_annotation_dir     = ch_linx_somatic_annotation_files
    linx_somatic_plot_dir           = ch_linx_somatic_plot_files
    linxreport_html                 = UMCCR_SASH.out.linxreport_html
    vcf2maf_maf                     = UMCCR_SASH.out.vcf2maf_maf
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUT DEFINITION

    Configures how each named output assigned in the entry `workflow {}` block's `publish:`
    section above is written under `--outdir`. This replicates, path-for-path, what the
    corresponding `publishDir`/`saveAs` blocks previously did in conf/modules.config for the
    15 "simple" (non filename-routing) processes. See conf/modules.config for the processes
    still using the legacy publishDir mechanism (filename-pattern-based routing/renaming is
    out of scope for this migration).

    `f.name` strips a process's own "output/" staging subdirectory to reproduce the old
    saveAs `fp.replaceFirst(/output\//, '')` behaviour. Where the process emits a directory
    (e.g. PURPLE's `purple/`), the directory is published as-is; the old dynamic (per-file)
    saveAs closures did not strip the directory's own name either, so this is unchanged
    behaviour, not a gap.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    smlv_somatic_rescue_vcf {
        path { meta, f -> f >> "${meta.key}/smlv_somatic/rescue/${f.name}" }
    }
    smlv_somatic_annotate_vcf {
        path { meta, f -> f >> "${meta.key}/smlv_somatic/annotate/${f.name}" }
    }
    smlv_somatic_filter_vcf {
        path { meta, vcf, tbi ->
            vcf >> "${meta.key}/smlv_somatic/filter/${vcf.name}"
            tbi >> "${meta.key}/smlv_somatic/filter/${tbi.name}"
        }
    }
    smlv_somatic_filter_vcf_filters {
        path { meta, f -> f >> "${meta.key}/smlv_somatic/filter/${f.name}" }
    }
    smlv_germline_prepare_vcf {
        path { meta, _f -> "${meta.key}/smlv_germline/prepare" }
    }
    sv_somatic_annotate_vcf {
        path { meta, f -> f >> "${meta.key}/sv_somatic/annotate/${f.name}" }
    }
    sv_somatic_prioritise_sv_tsv {
        path { meta, f -> f >> "${meta.key}/sv_somatic/prioritise/${f.name}" }
    }
    sv_somatic_prioritise_sv_vcf {
        path { meta, f -> f >> "${meta.key}/sv_somatic/prioritise/${f.name}" }
    }
    sv_somatic_prioritise_cnv_tsv {
        path { meta, f -> f >> "${meta.key}/sv_somatic/prioritise/${f.name}" }
    }
    purple_dir {
        path { meta, _dir -> "${meta.id}" }
    }
    sigrap_hrdetect_json {
        path { meta, _f -> "${meta.key}/sigrap/hrdetect" }
    }
    sigrap_mutpat_dir {
        path { meta, _dir -> "${meta.key}/sigrap" }
    }
    esvee_caller_dir {
        path { meta, _dir -> "${meta.key}/esvee" }
    }
    linx_germline_annotation_dir {
        path { meta, f -> f >> "${meta.id}/linx/germline_annotations/${f.name}" }
    }
    linx_somatic_annotation_dir {
        path { meta, f -> f >> "${meta.id}/linx/somatic_annotations/${f.name}" }
    }
    linx_somatic_plot_dir {
        path { meta, f, rel -> f >> "${meta.id}/linx/somatic_plots/${rel}" }
    }
    linxreport_html {
        path { meta, _f -> "${meta.key}" }
    }
    vcf2maf_maf {
        path { meta, _f -> "${meta.key}/vcf2maf" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
