process CUSTOM_EXTRACTTARBALL {
    label 'process_single'

    conda "conda-forge::tar=1.34"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(tarball)

    output:
    path "${meta.id}/", emit: extracted_dir
    path '.command.*' , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def strip = meta.strip_components != null ? meta.strip_components : 1
    def target = meta.subdir ? "${meta.id}/${meta.subdir}" : "${meta.id}"

    """
    mkdir -p ${target}
    tar ${args} -xzvf ${tarball} --strip-components ${strip} -C ${target}/
    """

    stub:
    """
    mkdir -p ${meta.id}/
    """
}
