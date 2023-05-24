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
        ch_hmf_data = createDataMap(params.hmfdata_paths, params.ref_data_path)
        ch_umccr_data = createDataMap(params.umccrdata_paths, params.ref_data_path)

        //
        // Prepare genome paths and info
        //
        genome = [
            fasta: getRefdataFile(params.genome.fasta, params.ref_data_path),
            fai: getRefdataFile(params.genome.fai, params.ref_data_path),
            dict: getRefdataFile(params.genome.dict, params.ref_data_path),
            version: '38',
        ]

    emit:
        genome                 = genome                         // map: Genome paths and info
        hmf_data               = ch_hmf_data                    // map: HMF data paths
        umccr_data             = ch_umccr_data                  // map: UMCCR data paths

        versions               = ch_versions                    // channel: [versions.yml]
}

def createDataMap(entries, ref_data_path) {
    return entries
        .collectEntries { name, path ->
            def ref_data_file = getRefdataFile(path, ref_data_path)
            return [name, ref_data_file]
        }
}

def getRefdataFile(filepath, ref_data_path) {
    def data_path_noslash = ref_data_path.toString().replaceAll('/$', '')
    return file("${data_path_noslash}/${filepath}", checkIfExists: true)
}
