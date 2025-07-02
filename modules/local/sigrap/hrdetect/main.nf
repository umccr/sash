process SIGRAP_HRDETECT {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/umccr/sigrap:0.2.0'

    input:
    tuple val(meta), path(smlv_somatic_vcf)
    tuple val(meta), path(sv_somatic_vcf)
    tuple val(meta), path(cnv_somatic_tsv)

    output:
    tuple val(meta), path('output/hrdetect.json.gz')                           , emit: hrdetect_json
    path 'versions.yml'                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    sigrap.R hrdetect \\
        --sample ${meta.subject_id} \\
        --snv ${smlv_somatic_vcf} \\
        --sv ${sv_somatic_vcf} \\
        --cnv ${cnv_somatic_tsv} \\
        --out  output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigrap: \$(sigrap.R --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/
    touch output/hrdetect.json.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
