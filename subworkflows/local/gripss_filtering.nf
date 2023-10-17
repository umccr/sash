//
// GRIPSS performs SV filtering.
//

include { GRIPSS_GERMLINE as GERMLINE } from '../../modules/local/gripss/germline/main'
include { GRIPSS_SOMATIC as SOMATIC   } from '../../modules/local/gripss/somatic/main'

workflow GRIPSS_FILTERING {
    take:
        // Sample inputs
        ch_inputs                // channel: [mandatory] [ meta ]
        ch_gridss                // channel: [mandatory] [ meta, gridss_vcf ]

        // Reference data
        genome_fasta             // channel: [mandatory] /path/to/genome_fasta
        genome_version           // channel: [mandatory] genome version
        genome_fai               // channel: [mandatory] /path/to/genome_fai
        breakend_pon             // channel: [mandatory] /path/to/breakend_pon
        breakpoint_pon           // channel: [mandatory] /path/to/breakpoint_pon
        known_fusions            // channel: [mandatory] /path/to/known_fusions
        repeatmasker_annotations // channel: [mandatory] /path/to/repeatmasker_annotations

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Create process input channel
        // channel: [ meta_gripss, gridss_vcf ]
        ch_gripss_inputs = ch_gridss
            .map { meta, gridss_vcf ->

                def meta_gripss = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.tumor_id,
                    normal_id: meta.normal_id,
                ]

                return [meta_gripss, gridss_vcf]
            }

        //
        // MODULE: GRIPSS germline
        //
        GERMLINE(
            ch_gripss_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            breakend_pon,
            breakpoint_pon,
            known_fusions,
            repeatmasker_annotations,
        )

        ch_versions = ch_versions.mix(GERMLINE.out.versions)

        //
        // MODULE: GRIPSS somatic
        //
        SOMATIC(
            ch_gripss_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            breakend_pon,
            breakpoint_pon,
            known_fusions,
            repeatmasker_annotations,
            [],
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, gripss_vcf, gripss_tbi ]
        ch_somatic_out = WorkflowSash.restoreMeta(SOMATIC.out.vcf, ch_inputs)
        ch_somatic_unfiltered_out = WorkflowSash.restoreMeta(SOMATIC.out.vcf_unfiltered, ch_inputs)
        ch_germline_out = WorkflowSash.restoreMeta(GERMLINE.out.vcf, ch_inputs)
        ch_germline_unfiltered_out = WorkflowSash.restoreMeta(GERMLINE.out.vcf_unfiltered, ch_inputs)

    emit:
        somatic             = ch_somatic_out             // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline            = ch_germline_out            // channel: [ meta, gripss_vcf, gripss_tbi ]
        somatic_unfiltered  = ch_somatic_unfiltered_out  // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline_unfiltered = ch_germline_unfiltered_out // channel: [ meta, gripss_vcf, gripss_tbi ]

        versions = ch_versions                           // channel: [ versions.yml ]
}
