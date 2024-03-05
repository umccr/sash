//
// Prepare reference data as required
//

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

def joinPath(a, b) {
    def a_noslash = file(a).toUriString().replaceAll('/$', '')
    return file("${a_noslash}/${b}", checkIfExists: true)
}
