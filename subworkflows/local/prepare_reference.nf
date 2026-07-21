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
        if (!misc_tarball_inputs) {
            error "No tarball entries found in miscdata_paths. At minimum pcgr_dir and vep_dir " +
                  "must be provided as .tar.gz/.tgz archives. Plain directory inputs are no longer supported."
        }

        ch_misc_data_inputs = Channel.fromList(misc_tarball_inputs)

        DECOMP_MISC_DATA(ch_misc_data_inputs)

        ch_misc_data_channel = DECOMP_MISC_DATA.out.extracted_dir
            .collect()
            .map { dir_list ->
                assert dir_list.size() == misc_tarball_inputs.size()
                def extracted_map = dir_list.collectEntries { dir ->
                    [(dir.getFileName().toString()): dir]
                }
                return addEnsemblFasta(createDataMap(params.miscdata_paths, params.ref_data_path) + extracted_map)
            }
            // .first() converts from a queue channel (collect().map{} emits once) to a value
            // channel, allowing it to broadcast to all per-sample processes in multi-sample runs.
            .first()

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
        misc_data_ch           = ch_misc_data_channel           // channel: Misc data with runtime-extracted dirs and fasta_ensembl

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
            def meta
            if (name == 'vep_dir') {
                // VEP cache: strip wrapper dir, extract into homo_sapiens subdir
                // Result: vep_dir/homo_sapiens/113_GRCh38/
                meta = [id: name, strip_components: 1, subdir: 'homo_sapiens']
            } else if (name == 'pcgr_dir') {
                // PCGR bundle: don't strip, tarball contains data/ directory
                // Result: pcgr_dir/data/grch38/
                meta = [id: name, strip_components: 0]
            } else {
                // Default: strip top-level wrapper directory
                meta = [id: name, strip_components: 1]
            }
            return [meta, tarball]
        }
}

def joinPath(a, b) {
    def a_noslash = file(a).toUriString().replaceAll('/$', '')
    return file("${a_noslash}/${b}", checkIfExists: true)
}

def addEnsemblFasta(map) {
    def fa_base = 'data/grch38/misc/fasta/assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    return map + [
        fasta_ensembl:     map['pcgr_dir'].resolve("${fa_base}.gz"),
        fasta_ensembl_fai: map['pcgr_dir'].resolve("${fa_base}.gz.fai"),
        fasta_ensembl_gzi: map['pcgr_dir'].resolve("${fa_base}.gz.gzi"),
    ]
}
