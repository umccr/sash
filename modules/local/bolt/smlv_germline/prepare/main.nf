process BOLT_SMLV_GERMLINE_PREPARE {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/umccr/bolt:0.2.15'

    input:
    tuple val(meta), path(smlv_vcf)
    path germline_predisposition_panel_regions_transcript

    output:
    tuple val(meta), path("${meta.normal_id}.prepared.vcf.gz"), emit: vcf
    path 'versions.yml'                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt smlv_germline prepare \\
        --vcf_fp ${smlv_vcf} \\
        --transcript_regions_fp ${germline_predisposition_panel_regions_transcript} \\
        --output_fp ${meta.normal_id}.prepared.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.normal_id}.prepared.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
