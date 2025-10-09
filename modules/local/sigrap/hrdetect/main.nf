process SIGRAP_HRDETECT {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/qclayssen/sigrap:0.2.0-dev-5'

    input:
    tuple val(meta), path(smlv_somatic_vcf), path(sv_somatic_vcf), path(cnv_somatic_tsv)

    output:
    tuple val(meta), path('hrdetect.json.gz')                 , emit: hrdetect_json
    path 'versions.yml'                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    sigrap.R hrdetect \\
        --sample ${meta.id} \\
        --snv ${smlv_somatic_vcf} \\
        --sv ${sv_somatic_vcf} \\
        --cnv ${cnv_somatic_tsv} \\
        --out  hrdetect.json.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigrap: \$(sigrap.R --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch hrdetect.json.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
