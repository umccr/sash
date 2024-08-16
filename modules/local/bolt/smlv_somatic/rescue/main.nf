process BOLT_SMLV_SOMATIC_RESCUE {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/scwatts/bolt:0.2.12'

    input:
    tuple val(meta), path(smlv_vcf), path(smlv_tbi), path(sage_smlv_vcf), path(sage_smlv_tbi)
    path umccr_hotspots

    output:
    tuple val(meta), path("output/${meta.tumor_id}.rescued.vcf.gz"), emit: vcf
    path 'versions.yml'                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcftools view -o ${meta.tumor_id}.main.dragen.vcf.gz ${smlv_vcf} chr{1..22} chr{X,Y,M}
    bcftools view -f PASS,. -o ${meta.tumor_id}.main.sage.filtered.vcf.gz ${sage_smlv_vcf} chr{1..22} chr{X,Y,M}

    bcftools index -t ${meta.tumor_id}.main.dragen.vcf.gz
    bcftools index -t ${meta.tumor_id}.main.sage.filtered.vcf.gz

    bolt smlv_somatic rescue \\
        --tumor_name ${meta.tumor_id} \\
        --vcf_fp ${meta.tumor_id}.main.dragen.vcf.gz \\
        --sage_vcf_fp ${meta.tumor_id}.main.sage.filtered.vcf.gz \\
        --hotspots_fp ${umccr_hotspots} \\
        --output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/
    touch output/${meta.tumor_id}.rescued.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
