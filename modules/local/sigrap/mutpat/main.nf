process SIGRAP_MUTPAT {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/qclayssen/sigrap:0.2.0-dev-6'

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
        --sample ${meta.id} \\
        --snv ${smlv_somatic_vcf} \\
        --rainfall \\
        --strand-bias \\
        --out mutpat/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigrap: \$(sigrap.R --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sigrap/mutpat/
    touch sigrap/mutpat/stub_output
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
