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
        def ch_esvee_somatic = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def vcf = "${base}/esvee/caller/${meta.tumor_id}.esvee.somatic.vcf.gz"
            return [meta, vcf]
        }
        def ch_esvee_germline = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def vcf = "${base}/esvee/caller/${meta.tumor_id}.esvee.germline.vcf.gz"
            return [meta, vcf]
        }
        def ch_esvee_somatic_unfiltered = ch_metas.map { meta ->
            def base = file(meta.oncoanalyser_dir).toUriString()
            def vcf = "${base}/esvee/caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz"
            return [meta, vcf]
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
        metas = ch_metas
        amber = ch_amber
        cobalt = ch_cobalt
        esvee_somatic = ch_esvee_somatic
        esvee_germline = ch_esvee_germline
        esvee_somatic_unfiltered = ch_esvee_somatic_unfiltered
        sage_somatic = ch_sage_somatic
        virusbreakend = ch_virusbreakend
        hrd = ch_input_hrd
        vcf_germline = ch_input_vcf_germline
        vcf_somatic = ch_input_vcf_somatic
}
