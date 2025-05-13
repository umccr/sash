workflow PREPARE_INPUT {
    take:
        ch_samplesheet

    main:
        def ch_metas = Channel.of(ch_samplesheet)
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

        // map oncoanalyser assets and DRAGEN outputs into channels
        def ch_amber = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            return [meta, "${base}/amber/"]
        }

        def ch_cobalt = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            return [meta, "${base}/cobalt/"]
        }

        def ch_call_inputs = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def dir = "${base}/esvee/prep/"
            def vcf = "${base}/esvee/depth_annotation/${meta.tumor_id}.esvee.ref_depth.vcf.gz"
            return [meta, dir, vcf]
        }

        def ch_sage_somatic = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def vcf = "${base}/sage/somatic/${meta.tumor_id}.sage.somatic.vcf.gz"
            return [meta, vcf, "${vcf}.tbi"]
        }

        def ch_virusbreakend = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            return [meta, "${base}/virusbreakend/"]
        }

        def ch_input_hrd = ch_metas.map { meta ->
            def base = file(meta.dragen_somatic_dir).toUriString()
            return [meta, "${base}/${meta.tumor_id}.hrdscore.csv"]
        }

        def ch_input_vcf_germline = ch_metas.map { meta ->
            def base = file(meta.dragen_germline_dir).toUriString()
            return [meta, "${base}/${meta.normal_id}.hard-filtered.vcf.gz"]
        }

        def ch_input_vcf_somatic = ch_metas.map { meta ->
            def base = file(meta.dragen_somatic_dir).toUriString()
            def vcf = "${base}/${meta.tumor_id}.hard-filtered.vcf.gz"
            return [meta, vcf, "${vcf}.tbi"]
        }
    emit:
        // Meta information for each sample
        metas = ch_metas                   // channel: [ meta ]

        // Oncoanalyser inputs
        amber = ch_amber                   // channel: [ meta, amber_dir ]
        cobalt = ch_cobalt                 // channel: [ meta, cobalt_dir ]
        sage_somatic = ch_sage_somatic     // channel: [ meta, vcf, tbi ]
        virusbreakend = ch_virusbreakend   // channel: [ meta, virusbreakend_dir ]

        // DRAGEN inputs
        hrd = ch_input_hrd                 // channel: [ meta, hrdscore_csv ]
        vcf_germline = ch_input_vcf_germline // channel: [ meta, germline_vcf ]
        vcf_somatic = ch_input_vcf_somatic   // channel: [ meta, somatic_vcf, tbi ]

        // eSVee inputs
        call_inputs = ch_call_inputs             // channel: [ meta, esvee_prep_dir, esvee_ref_depth_vcf ]
}
