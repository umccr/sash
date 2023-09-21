process BOLT_SMLV_GERMLINE_REPORT {
    tag "${meta.key}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.1.7-pcgr'

    input:
    tuple val(meta), path(smlv_vcf), path(smlv_unfiltered_vcf)
    path germline_predisposition_panel_genes
    path pcgr_data_dir

    output:
    tuple val(meta), path('cpsr_output/')         , emit: cpsr_dir
    tuple val(meta), path("*.variant_counts.yaml"), emit: variant_counts
    tuple val(meta), path("*.bcftools_stats.txt") , emit: bcftools_stats
    path "*.annotations.vcf.gz"                   , emit: vcf
    path 'versions.yml'                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt smlv_germline report \\
        --normal_name ${meta.id} \\
        --vcf_fp ${smlv_vcf} \\
        --vcf_unfiltered_fp ${smlv_unfiltered_vcf} \\
        --pcgr_conda pcgr \\
        --pcgrr_conda pcgrr \\
        --germline_panel_list_fp ${germline_predisposition_panel_genes} \\
        --pcgr_data_dir ${pcgr_data_dir} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cpsr_output/
    touch ${meta.id}.germline.bcftools_stats.txt
    touch ${meta.id}.germline.variant_counts.yaml
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
