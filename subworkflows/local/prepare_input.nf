import nextflow.Nextflow

workflow PREPARE_INPUT {
    take:
        ch_samplesheet

    main:
        // Parse samplesheet and group entries by sample ID
        // channel: [ meta ]
        ch_metas = Channel.of(ch_samplesheet)
            .splitCsv(header: true)
            .map { [it.id, it] }
            .groupTuple()
            .map { key, entries ->
                def meta = [id: key]
                entries.each {
                    switch (it.filetype) {
                        case 'oncoanalyser_dir':
                            break;
                        default:
                            log.error "got bad filetype: ${it.filetype}"
                            Nextflow.exit(1)
                    }

                    meta[it.filetype] = it.filepath

                    if (it.tumor_id && !meta.containsKey('tumor_id')) meta.tumor_id = it.tumor_id
                    if (it.normal_id && !meta.containsKey('normal_id')) meta.normal_id = it.normal_id

                    // Sample name
                    if (! meta.containsKey('subject_id')) {
                        meta.subject_id = it.subject_name
                    } else if (meta.subject_id != it.subject_name) {
                        log.error "expected ${meta.subject_id} as subject name but got ${it.subject_name}"
                        Nextflow.exit(1)
                    }

                }

                return meta
            }

        // Map oncoanalyser assets into channels

        // AMBER: copy number segmentation data
        // channel: [ meta, amber_dir ]
        ch_amber = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def amber_dir = "${base}/amber/"
            if (!file(amber_dir).exists()) {
                log.error "AMBER directory not found for ${meta.id}: ${amber_dir}"
                Nextflow.exit(1)
            }
            return [meta, amber_dir]
        }

        // COBALT: read depth ratio data
        // channel: [ meta, cobalt_dir ]
        ch_cobalt = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def cobalt_dir = "${base}/cobalt/"
            if (!file(cobalt_dir).exists()) {
                log.error "COBALT directory not found for ${meta.id}: ${cobalt_dir}"
                Nextflow.exit(1)
            }
            return [meta, cobalt_dir]
        }

        // eSVee: structural variant calling inputs
        // channel: [ meta_esvee, esvee_ref_depth_vcf, esvee_prep_dir ]
        ch_call_inputs = ch_metas.map { meta ->
            def meta_esvee = [
                key: meta.id,
                id: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
                sample_id: meta.tumor_id
            ]

            def base = file(meta.oncoanalyser_dir).toUriString()
            def esvee_ref_depth_vcf = "${base}/esvee/${meta.tumor_id}.esvee.ref_depth.vcf.gz"
            def esvee_prep_dir = "${base}/esvee/"
            if (!file(esvee_ref_depth_vcf).exists()) {
                log.error "eSVee depth VCF not found for ${meta.id}: ${esvee_ref_depth_vcf}"
                Nextflow.exit(1)
            }
            if (!file(esvee_prep_dir).exists()) {
                log.error "eSVee directory not found for ${meta.id}: ${esvee_prep_dir}"
                Nextflow.exit(1)
            }
            return [meta_esvee, esvee_ref_depth_vcf, esvee_prep_dir]
        }

        // SAGE: somatic small variant calls
        // channel: [ meta, sage_somatic_vcf, sage_somatic_tbi ]
        // NOTE(QC): OA v2.3.0+ uses sage/ instead of sage_calling/ — try both paths
        ch_sage_somatic = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def sage_somatic_vcf = "${base}/sage_calling/somatic/${meta.tumor_id}.sage.somatic.vcf.gz"
            if (!file(sage_somatic_vcf).exists()) {
                sage_somatic_vcf = "${base}/sage/somatic/${meta.tumor_id}.sage.somatic.vcf.gz"
            }
            def sage_somatic_tbi = "${sage_somatic_vcf}.tbi"
            if (!file(sage_somatic_vcf).exists()) {
                log.error "SAGE somatic VCF not found for ${meta.id} (tried sage_calling/ and sage/)"
                Nextflow.exit(1)
            }
            if (!file(sage_somatic_tbi).exists()) {
                log.error "SAGE somatic TBI not found for ${meta.id}: ${sage_somatic_tbi}"
                Nextflow.exit(1)
            }
            return [meta, sage_somatic_vcf, sage_somatic_tbi]
        }

        // VirusBreakend: viral integration detection
        // channel: [ meta, virusbreakend_dir ]
        ch_virusbreakend = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def virusbreakend_dir = "${base}/virusbreakend/"
            if (!file(virusbreakend_dir).exists()) {
                log.error "VirusBreakend directory not found for ${meta.id}: ${virusbreakend_dir}"
                Nextflow.exit(1)
            }
            return [meta, virusbreakend_dir]
        }

        // CHORD: homologous recombination deficiency prediction
        // channel: [ meta, chord_prediction_tsv ]
        ch_chord = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def chord_prediction_tsv = "${base}/chord/${meta.tumor_id}.chord.prediction.tsv"
            if (!file(chord_prediction_tsv).exists()) {
                log.error "CHORD prediction file not found for ${meta.id}: ${chord_prediction_tsv}"
                Nextflow.exit(1)
            }
            return [meta, chord_prediction_tsv]
        }

        // Germline variants — use PAVE germline from oncoanalyser output if available
        // channel: [ meta, germline_vcf ]
        ch_input_vcf_germline = ch_metas.map { meta ->
            if (!meta.normal_id) {
                return [meta, []]
            }
            def base = file(meta.oncoanalyser_dir).toUriString()
            def pave_germline_vcf = "${base}/pave/${meta.normal_id}.pave.germline.vcf.gz"
            if (file(pave_germline_vcf).exists()) {
                return [meta, pave_germline_vcf]
            }
            log.warn "No PAVE germline VCF found for ${meta.id} — germline report will be skipped"
            return [meta, []]
        }
    emit:
        // Sample metadata
        metas            = ch_metas                   // channel: [ meta ]

        // oncoanalyser channels
        amber            = ch_amber                   // channel: [ meta, amber_dir ]
        cobalt           = ch_cobalt                  // channel: [ meta, cobalt_dir ]
        sage_somatic     = ch_sage_somatic            // channel: [ meta, sage_somatic_vcf, sage_somatic_tbi ]
        virusbreakend    = ch_virusbreakend           // channel: [ meta, virusbreakend_dir ]
        chord            = ch_chord                   // channel: [ meta, chord_prediction_tsv ]
        call_inputs      = ch_call_inputs             // channel: [ meta_esvee, esvee_ref_depth_vcf, esvee_prep_dir ]

        // OA germline
        vcf_germline     = ch_input_vcf_germline      // channel: [ meta, pave_germline_vcf ]
}
