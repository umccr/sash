process LINXREPORT {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/scwatts/r-linxreport:1.0.0--0'

    input:
    tuple val(meta), path(linx_annotation_dir), path(linx_visualiser_dir)

    output:
    tuple val(meta), path('*_linx.html'), emit: html
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Create input directory if it doesn't exist for linxreport
    if [[ ! -e ${linx_visualiser_dir} ]]; then
        mkdir -p ${linx_visualiser_dir};
    fi;

    linxreport.R \\
        ${args} \\
        --sample ${meta.sample_id} \\
        --plot ${linx_visualiser_dir} \\
        --table ${linx_annotation_dir} \\
        --out ${meta.sample_id}_linx.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/^R version \\([0-9.]\\+\\).\\+/\\1/')
        linxreport: \$(linxreport.R --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}_linx.html
    echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
    """
}
