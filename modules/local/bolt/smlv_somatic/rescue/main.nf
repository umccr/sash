process BOLT_SMLV_SOMATIC_RESCUE {
    tag "${meta.key}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.1.8'

    input:
    tuple val(meta), path(smlv_vcf), path(smlv_tbi), path(sage_smlv_vcf), path(sage_smlv_tbi)
    path umccr_hotspots

    output:
    tuple val(meta), path("${meta.tumor_id}.rescued.vcf.gz"), emit: vcf
    path 'versions.yml'                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcftools view -o ${meta.tumor_id}.main.dragen.vcf.gz ${smlv_vcf} chr{1..22} chr{X,Y,M}
    bcftools view -o ${meta.tumor_id}.main.sage.vcf.gz ${sage_smlv_vcf} chr{1..22} chr{X,Y,M}

    bcftools index -t ${meta.tumor_id}.main.dragen.vcf.gz
    bcftools index -t ${meta.tumor_id}.main.sage.vcf.gz

    bolt smlv_somatic rescue \\
        --vcf_fp ${meta.tumor_id}.main.dragen.vcf.gz \\
        --sage_vcf_fp ${meta.tumor_id}.main.sage.vcf.gz \\
        --hotspots_fp ${umccr_hotspots} \\
        --output_fp ${meta.tumor_id}.rescued.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.rescued.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
