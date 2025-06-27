process BOLT_SMLV_SOMATIC_FILTER {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/umccr/bolt:0.3.0-dev-1'

    input:
    tuple val(meta), path(smlv_vcf)

    output:
    tuple val(meta), path("output/${meta.tumor_id}*pass.vcf.gz"), path("output/${meta.tumor_id}*pass.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("output/${meta.tumor_id}*filters_set.vcf.gz"),                                           emit: vcf_filters
    path 'versions.yml',                                                                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt smlv_somatic filter \\
        --tumor_name ${meta.tumor_id} \\
        --vcf_fp ${smlv_vcf} \\
        --output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/
    touch output/${meta.tumor_id}.filters_set.pass.vcf.gz
    touch output/${meta.tumor_id}.filters_set.pass.vcf.gz.tbi
    touch output/${meta.tumor_id}.filters_set.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
