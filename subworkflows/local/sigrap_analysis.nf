//
// SIGRAP performs mutational signature analysis
//

include { SIGRAP_CHORD    } from '../../modules/local/sigrap/chord/main'
include { SIGRAP_HRDETECT } from '../../modules/local/sigrap/hrdetect/main'
include { SIGRAP_MUTPAT   } from '../../modules/local/sigrap/mutpat/main'

workflow SIGRAP_ANALYSIS {
    take:
        // Sample data
        ch_inputs                    // channel: [mandatory] [ meta ]
        ch_smlv_somatic              // channel: [mandatory] [ meta, smlv_somatic_vcf ]
        ch_sv_somatic_vcf            // channel: [mandatory] [ meta, sv_vcf ]
        ch_sv_somatic_cnv_tsv        // channel: [mandatory] [ meta, cnv_tsv ]

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        SIGRAP_CHORD(
            ch_smlv_somatic,
            ch_sv_somatic_vcf
        )

        // channel: [ meta, chord_json ]
        ch_chord_out = WorkflowSash.restoreMeta(SIGRAP_CHORD.out.chord_json, ch_inputs)
        ch_versions = ch_versions.mix(SIGRAP_CHORD.out.versions)

        SIGRAP_HRDETECT(
            ch_smlv_somatic,
            ch_sv_somatic_vcf,
            ch_sv_somatic_cnv_tsv,
        )

        // channel: [ meta, hrdetect_json ]
        ch_hrdetect_out = WorkflowSash.restoreMeta(SIGRAP_HRDETECT.out.hrdetect_json, ch_inputs)
        ch_versions = ch_versions.mix(SIGRAP_HRDETECT.out.versions)
        
        SIGRAP_MUTPAT(
            ch_smlv_somatic
        )

        // channel: [ meta, mutpat_output ]
        ch_mutpat_out = WorkflowSash.restoreMeta(SIGRAP_MUTPAT.out.mutpat_output, ch_inputs)
        ch_versions = ch_versions.mix(SIGRAP_MUTPAT.out.versions)

    emit:
        chord     = ch_chord_out     // channel: [ meta, chord_json ]
        hrdetect  = ch_hrdetect_out  // channel: [ meta, hrdetect_json ]
        mutpat    = ch_mutpat_out    // channel: [ meta, mutpat_output ]

        versions  = ch_versions      // channel: [ versions.yml ]
}
