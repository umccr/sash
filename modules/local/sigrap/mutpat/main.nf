process SIGRAP_MUTPAT {
    tag "${meta.id}"
    label 'process_medium_memory'
    label 'process_long'

    container 'docker.io/qclayssen/sigrap:0.2.0-dev-7'

    input:
    tuple val(meta), path(smlv_somatic_vcf)

    output:
    tuple val(meta), path('mutpat/')                   , emit: mutpat_output
    path 'versions.yml'                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    sigrap.R mutpat \\
        --sample ${meta.tumor_id} \\
        --snv ${smlv_somatic_vcf} \\
        --rainfall \\
        --strand-bias \\
        --predefined-dbs-mbs \\
        --out mutpat/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigrap: \$(sigrap.R --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p mutpat/
    touch mutpat/stub_output
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
