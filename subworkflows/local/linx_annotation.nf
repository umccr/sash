//
// LINX annotates and interprets structural variants
//

include { LINX_GERMLINE as GERMLINE } from '../../modules/local/linx/germline/main'
include { LINX_SOMATIC as SOMATIC   } from '../../modules/local/linx/somatic/main'

workflow LINX_ANNOTATION {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_purple              // channel: [mandatory] [ meta, purple_dir ]

        // Reference data
        genome_version         // channel: [mandatory] genome version
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
        known_fusion_data      // channel: [mandatory] /path/to/known_fusion_data
        driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Create inputs and create process-specific meta
        // channel: [ meta, sv_vcf ]
        ch_linx_inputs_germline = ch_purple
            .map { meta, purple_dir ->

                def sv_vcf = file(purple_dir).resolve("${meta.tumor_id}.purple.sv.germline.vcf.gz")

                def meta_linx = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.tumor_id,
                ]

                return [meta_linx, sv_vcf]
            }

        // channel: [ meta_linx, purple_dir ]
        ch_linx_inputs_somatic = ch_purple
            .map { meta, purple_dir ->
                def meta_linx = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.tumor_id,
                ]
                return [meta_linx, purple_dir]
            }

        GERMLINE(
            ch_linx_inputs_germline,
            genome_version,
            ensembl_data_resources,
            driver_gene_panel,
        )

        ch_versions = ch_versions.mix(GERMLINE.out.versions)

        SOMATIC(
            ch_linx_inputs_somatic,
            genome_version,
            ensembl_data_resources,
            known_fusion_data,
            driver_gene_panel,
            [],
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)

        ch_linx_somatic_out = WorkflowSash.restoreMeta(SOMATIC.out.annotation_dir, ch_inputs)
        ch_linx_germline_out = WorkflowSash.restoreMeta(GERMLINE.out.annotation_dir, ch_inputs)

    emit:
        somatic       = ch_linx_somatic_out  // channel: [ meta, linx_annotation_dir ]
        germline      = ch_linx_germline_out // channel: [ meta, linx_annotation_dir ]

        versions      = ch_versions          // channel: [ versions.yml ]
}
