//
// Prepare reference data as required
//

include { CUSTOM_EXTRACTTARBALL as DECOMP_MISC_DATA } from '../../modules/local/custom/extract_tarball/main'

workflow PREPARE_REFERENCE {
    take:

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        //
        // Set UMCCR and HMF reference data paths
        //
        umccr_reference_data_path = joinPath(params.ref_data_path, "umccr_reference_data/v${params.data_versions.umccr_reference_data}/")

        ch_hmf_data = createDataMap(params.hmfdata_paths, params.ref_data_path)
        ch_umccr_data = createDataMap(params.umccrdata_paths, umccr_reference_data_path)
        ch_misc_data = createDataMap(params.miscdata_paths, params.ref_data_path)

        //
        // Extract tarball resources (e.g. PCGR data, VEP cache) when provided as .tar.gz/.tgz
        //
        misc_tarball_inputs = getTarballInputs(params.miscdata_paths, params.ref_data_path)
        if (misc_tarball_inputs) {
            ch_misc_data_inputs = Channel.fromList(misc_tarball_inputs)
            
            DECOMP_MISC_DATA(ch_misc_data_inputs)

            ch_misc_data_extracted = DECOMP_MISC_DATA.out.extracted_dir
                .collect()
                .map { dir_list ->
                    // Convert list of directories to a map
                    // For vep_dir extraction (which becomes 'homo_sapiens' dir), map back to 'vep_dir' key
                    def extracted_map = dir_list.collectEntries { dir ->
                        def dirname = dir.getFileName().toString()
                        def key = dirname == 'homo_sapiens' ? 'vep_dir' : dirname
                        [(key): dir]
                    }
                    // Merge extracted data with existing misc_data map
                    return createDataMap(params.miscdata_paths, params.ref_data_path) + extracted_map
                }
            
            ch_misc_data = ch_misc_data_extracted
        }

        //
        // Prepare genome paths and info
        //
        genome = [
            fasta: joinPath(params.ref_data_path, params.genome.fasta),
            fai: joinPath(params.ref_data_path, params.genome.fai),
            dict: joinPath(params.ref_data_path, params.genome.dict),
            version: '38',
        ]

    emit:
        genome                 = genome                         // map: Genome paths and info
        hmf_data               = ch_hmf_data                    // map: HMF data paths
        umccr_data             = ch_umccr_data                  // map: UMCCR data paths
        misc_data              = ch_misc_data                   // map: Misc data paths

        versions               = ch_versions                    // channel: [versions.yml]
}

def createDataMap(entries, ref_data_base_path) {
    return entries
        .collectEntries { name, path ->
            def ref_data_file = joinPath(ref_data_base_path, path)
            return [name, ref_data_file]
        }
}

def getTarballInputs(entries, ref_data_base_path) {
    return entries
        .findAll { name, relpath ->
            if (!relpath) {
                return false
            }
            def rel = relpath.toString()
            return rel.endsWith('.tar.gz') || rel.endsWith('.tgz')
        }
        .collect { name, relpath ->
            def tarball = joinPath(ref_data_base_path, relpath)
            // For VEP cache, use 'homo_sapiens' as the directory name (parent of version dir)
            // PCGR expects --vep_dir to point to parent containing homo_sapiens/113_GRCh38/
            def id = name == 'vep_dir' ? 'homo_sapiens' : name
            def meta = [id: id, strip_components: 1]
            return [meta, tarball]
        }
}

def joinPath(a, b) {
    def a_noslash = file(a).toUriString().replaceAll('/$', '')
    return file("${a_noslash}/${b}", checkIfExists: true)
}
