import java.nio.file.NoSuchFileException

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
                    // Filetype
                    switch (it.filetype) {
                        case 'dragen_somatic_dir':
                            meta.tumor_id = it.sample_name;
                            break;
                        case 'dragen_germline_dir':
                            meta.normal_id = it.sample_name;
                            break;
                        case 'oncoanalyser_dir':
                            break;
                        default:
                            log.error "\nERROR: got bad filetype: ${it.filetype}"
                            System.exit(1)
                    }

                    meta[it.filetype] = it.filepath

                    // Sample name
                    if (! meta.containsKey('subject_id')) {
                        meta.subject_id = it.subject_name
                    } else if (meta.subject_id != it.subject_name) {
                        log.error "\nERROR: expected ${meta.subject_id} as subject name but got ${it.subject_name}"
                        System.exit(1)
                    }

                }

                return meta
            }

        // Map oncoanalyser assets and DRAGEN outputs into channels

        // AMBER: copy number segmentation data
        // channel: [ meta, amber_dir ]
        ch_amber = ch_metas.map { meta ->
            def amber_dir = Utils.resolve_input_path(log, meta, meta.oncoanalyser_dir, "amber/", "AMBER directory")
            return [meta, amber_dir]
        }

        // COBALT: read depth ratio data
        // channel: [ meta, cobalt_dir ]
        ch_cobalt = ch_metas.map { meta ->
            def cobalt_dir = Utils.resolve_input_path(log, meta, meta.oncoanalyser_dir, "cobalt/", "COBALT directory")
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

            def esvee_vcf = Utils.resolve_input_path(log, meta, meta.oncoanalyser_dir, "esvee/${meta.tumor_id}.esvee.ref_depth.vcf.gz", "eSVee depth VCF")
            def esvee_dir = Utils.resolve_input_path(log, meta, meta.oncoanalyser_dir, "esvee/", "eSVee directory")
            return [meta_esvee, esvee_vcf, esvee_dir]
        }

        // SAGE: somatic small variant calls
        // channel: [ meta, sage_somatic_vcf, sage_somatic_tbi ]
        ch_sage_somatic = ch_metas.map { meta ->
            def sage_vcf = Utils.resolve_input_path(log, meta, meta.oncoanalyser_dir, "sage_calling/somatic/${meta.tumor_id}.sage.somatic.vcf.gz", "SAGE somatic VCF")
            def sage_tbi = Utils.resolve_input_path(log, meta, meta.oncoanalyser_dir, "sage_calling/somatic/${meta.tumor_id}.sage.somatic.vcf.gz.tbi", "SAGE somatic TBI")
            return [meta, sage_vcf, sage_tbi]
        }

        // VirusBreakend: viral integration detection
        // channel: [ meta, virusbreakend_dir ]
        ch_virusbreakend = ch_metas.map { meta ->
            def virusbreakend_dir = Utils.resolve_input_path(log, meta, meta.oncoanalyser_dir, "virusbreakend/", "VirusBreakend directory")
            return [meta, virusbreakend_dir]
        }

                // HRD: homologous recombination deficiency scores
        // channel: [ meta, hrdscore_csv ]
        ch_input_hrd = ch_metas.map { meta ->
            def hrdscore_csv = Utils.resolve_input_path(log, meta, meta.dragen_somatic_dir, "${meta.tumor_id}.hrdscore.csv", "HRD score CSV", true)
            return [meta, hrdscore_csv ?: []]  // Convert null to empty list
        }

        // DRAGEN germline variants
        // channel: [ meta, dragen_germline_vcf ]
        ch_input_vcf_germline = ch_metas.map { meta ->
            def dragen_germline_vcf = Utils.resolve_input_path(log, meta, meta.dragen_germline_dir, "${meta.normal_id}.hard-filtered.vcf.gz", "DRAGEN germline VCF")
            return [meta, dragen_germline_vcf]
        }

        // DRAGEN somatic variants
        // channel: [ meta, dragen_somatic_vcf, dragen_somatic_tbi ]
        ch_input_vcf_somatic = ch_metas.map { meta ->
            def dragen_somatic_vcf = Utils.resolve_input_path(log, meta, meta.dragen_somatic_dir, "${meta.tumor_id}.hard-filtered.vcf.gz", "DRAGEN somatic VCF")
            def dragen_somatic_tbi = Utils.resolve_input_path(log, meta, meta.dragen_somatic_dir, "${meta.tumor_id}.hard-filtered.vcf.gz.tbi", "DRAGEN somatic TBI")
            return [meta, dragen_somatic_vcf, dragen_somatic_tbi]
        }
    emit:
        // Sample metadata
        metas            = ch_metas                   // channel: [ meta ]

        // oncoanalyser channels
        amber            = ch_amber                   // channel: [ meta, amber_dir ]
        cobalt           = ch_cobalt                  // channel: [ meta, cobalt_dir ]
        sage_somatic     = ch_sage_somatic            // channel: [ meta, sage_somatic_vcf, sage_somatic_tbi ]
        virusbreakend    = ch_virusbreakend           // channel: [ meta, virusbreakend_dir ]
        call_inputs      = ch_call_inputs             // channel: [ meta_esvee, esvee_ref_depth_vcf, esvee_prep_dir ]

        // DRAGEN channels
        hrd              = ch_input_hrd               // channel: [ meta, hrdscore_csv ]
        vcf_germline     = ch_input_vcf_germline      // channel: [ meta, dragen_germline_vcf ]
        vcf_somatic      = ch_input_vcf_somatic       // channel: [ meta, dragen_somatic_vcf, dragen_somatic_tbi ]
}
