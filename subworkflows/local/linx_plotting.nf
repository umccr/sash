//
// LINX plotting visualises clusters structural variants
//

include { LINXREPORT as REPORT          } from '../../modules/local/linxreport/main'
include { LINX_VISUALISER as VISUALISER } from '../../modules/local/linx/visualiser/main'

workflow LINX_PLOTTING {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_annotations         // channel: [mandatory] [ meta, linx_annotation_dir ]

        // Reference data
        genome_version         // channel: [mandatory] genome version
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta_linx, linx_annotation_dir ]
        ch_linx_visualiser_inputs = ch_annotations
            .map { meta, anno_dir ->
                def meta_linx = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.tumor_id,
                ]
                return [meta_linx, anno_dir]
            }

        VISUALISER(
            ch_linx_visualiser_inputs,
            genome_version,
            ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(VISUALISER.out.versions)

        ch_visualiser_out = WorkflowSash.restoreMeta(VISUALISER.out.plots, ch_inputs)

        // Create inputs and create process-specific meta
        // channel: [ meta_linxreport, linx_annotation_dir, linx_plot_dir ]
        ch_linxreport_inputs = WorkflowSash.groupByMeta(
            ch_annotations,
            ch_visualiser_out,
        )
            .map { meta, anno_dir, plot_dir ->
                def meta_linxreport = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.tumor_id,
                ]
                return [meta_linxreport, anno_dir, plot_dir]
            }

        REPORT(
            ch_linxreport_inputs,
        )

        ch_versions = ch_versions.mix(REPORT.out.versions)

    emit:
        plot_dir = ch_visualiser_out // channel: [ meta, plot_dir ]

        versions = ch_versions       // channel: [ versions.yml ]
}
