process BOLT_SMLV_SOMATIC_FILTER {
    tag "${meta.key}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.1.0'

    input:
    tuple val(meta), path(smlv_vcf)

    output:
    tuple val(meta), path("${meta.tumor_id}*pass.vcf.gz")       , emit: vcf
    tuple val(meta), path("${meta.tumor_id}*filters_set.vcf.gz"), emit: vcf_unfiltered
    path 'versions.yml'                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt smlv_somatic filter \\
        --tumor_name ${meta.tumor_id} \\
        --vcf_fp ${smlv_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.filters_set.vcf.gz
    touch ${meta.tumor_id}.filters_set.pass.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
