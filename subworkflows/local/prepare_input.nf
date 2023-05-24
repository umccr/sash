workflow PREPARE_INPUT {
    take:
        ch_samplesheet

    main:
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

    emit:
        metas = ch_metas
}
