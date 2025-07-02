//
// SMLV somatic processing performs rescue, annotation, filtering, and PAVE annotation
//

include { BOLT_SMLV_SOMATIC_ANNOTATE } from '../../modules/local/bolt/smlv_somatic/annotate/main'
include { BOLT_SMLV_SOMATIC_FILTER   } from '../../modules/local/bolt/smlv_somatic/filter/main'
include { BOLT_SMLV_SOMATIC_RESCUE   } from '../../modules/local/bolt/smlv_somatic/rescue/main'
include { PAVE_SOMATIC               } from '../../modules/local/pave/somatic/main'

workflow SMLV_SOMATIC_PROCESSING {
    take:
        // Sample data
        ch_inputs                    // channel: [mandatory] [ meta ]
        ch_input_vcf_somatic         // channel: [mandatory] [ meta, dragen_somatic_vcf, dragen_somatic_tbi ]
        ch_sage_somatic              // channel: [mandatory] [ meta, sage_somatic_vcf, sage_somatic_tbi ]

        // Reference data
        genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
        genome_version               // channel: [mandatory] genome version
        genome_fai                   // channel: [mandatory] /path/to/genome_fai
        umccr_hotspots               // channel: [mandatory] /path/to/umccr_hotspots
        somatic_panel_regions_gene   // channel: [mandatory] /path/to/somatic_panel_regions_gene
        annotations_dir              // channel: [mandatory] /path/to/annotations_dir
        pon_dir                      // channel: [mandatory] /path/to/pon_dir
        pcgr_dir                     // channel: [mandatory] /path/to/pcgr_dir
        vep_dir                      // channel: [mandatory] /path/to/vep_dir
        clinvar_annotations          // channel: [mandatory] /path/to/clinvar_annotations
        segment_mappability          // channel: [mandatory] /path/to/segment_mappability
        driver_gene_panel            // channel: [mandatory] /path/to/driver_gene_panel
        ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Prepare rescue inputs with meta transformation
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
            umccr_hotspots,
        )

        ch_versions = ch_versions.mix(BOLT_SMLV_SOMATIC_RESCUE.out.versions)

        BOLT_SMLV_SOMATIC_ANNOTATE(
            BOLT_SMLV_SOMATIC_RESCUE.out.vcf,
            somatic_panel_regions_gene,
            annotations_dir,
            pon_dir,
            pcgr_dir,
            vep_dir
        )

        ch_versions = ch_versions.mix(BOLT_SMLV_SOMATIC_ANNOTATE.out.versions)

        BOLT_SMLV_SOMATIC_FILTER(
            BOLT_SMLV_SOMATIC_ANNOTATE.out.vcf,
        )

        ch_versions = ch_versions.mix(BOLT_SMLV_SOMATIC_FILTER.out.versions)

        // Restore meta and create clean outputs
        // channel: [ meta, smlv_somatic_vcf ]
        ch_smlv_somatic_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_FILTER.out.vcf, ch_inputs)
            .map { meta, vcf, tbi -> [meta, vcf] }
        // channel: [ meta, smlv_somatic_filters_vcf ]
        ch_smlv_somatic_filters_out = WorkflowSash.restoreMeta(BOLT_SMLV_SOMATIC_FILTER.out.vcf_filters, ch_inputs)

        // NOTE(SW): PAVE is run so that we obtain a complete PURPLE driver catalogue
        PAVE_SOMATIC(
            BOLT_SMLV_SOMATIC_FILTER.out.vcf,
            genome_fasta,
            genome_version,
            genome_fai,
            clinvar_annotations,
            segment_mappability,
            driver_gene_panel,
            ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(PAVE_SOMATIC.out.versions)

        // channel: [ meta, pave_somatic_vcf ]
        ch_pave_somatic_out = WorkflowSash.restoreMeta(PAVE_SOMATIC.out.vcf, ch_inputs)

    emit:
        smlv_somatic         = ch_smlv_somatic_out         // channel: [ meta, smlv_somatic_vcf ]
        smlv_somatic_filters = ch_smlv_somatic_filters_out // channel: [ meta, smlv_somatic_filters_vcf ]
        pave_somatic         = ch_pave_somatic_out         // channel: [ meta, pave_somatic_vcf ]

        versions             = ch_versions                 // channel: [ versions.yml ]
}
