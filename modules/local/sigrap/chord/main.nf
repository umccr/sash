process SIGRAP_CHORD {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/umccr/sigrap:0.2.0'

    input:
    tuple val(meta), path(smlv_somatic_vcf)
    tuple val(meta), path(sv_somatic_vcf)

    output:
    tuple val(meta), path('output/chord.json.gz')                           , emit: chord_json
    path 'versions.yml'                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    sigrap.R chord \\
        --sample ${meta.subject_id} \\
        --snv ${smlv_somatic_vcf} \\
        --sv ${sv_somatic_vcf} \\
        --out  output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigrap: \$(sigrap --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/
    touch output/chord.json.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
