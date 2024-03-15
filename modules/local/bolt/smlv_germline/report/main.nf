process BOLT_SMLV_GERMLINE_REPORT {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/scwatts/bolt:0.2.10-pcgr'

    input:
    tuple val(meta), path(smlv_vcf), path(smlv_unfiltered_vcf)
    path germline_predisposition_panel_genes
    path pcgr_data_dir

    output:
    tuple val(meta), path("output/*.variant_counts_type.yaml"), emit: counts_type
    tuple val(meta), path("output/*.bcftools_stats.txt")      , emit: bcftools_stats
    path 'output/cpsr/'                                       , emit: cpsr_dir
    path 'output/*.cpsr.grch38.html'                          , emit: cpsr_report
    path "output/*.annotations.vcf.gz"                        , emit: vcf
    path 'versions.yml'                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt smlv_germline report \\
        --normal_name ${meta.normal_id} \\
        --vcf_fp ${smlv_vcf} \\
        --vcf_unfiltered_fp ${smlv_unfiltered_vcf} \\
        --pcgr_conda pcgr \\
        --pcgrr_conda pcgrr \\
        --germline_panel_list_fp ${germline_predisposition_panel_genes} \\
        --pcgr_data_dir ${pcgr_data_dir} \\
        --threads ${task.cpus} \\
        --output_dir output/

    mv output/cpsr/${meta.normal_id}.cpsr.grch38.html output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/cpsr/
    touch output/cpsr/${meta.normal_id}.cpsr.grch38.html
    touch output/cpsr/${meta.normal_id}.annotations.vcf.gz
    touch output/${meta.normal_id}.germline.variant_counts_type.yaml
    touch output/${meta.normal_id}.germline.bcftools_stats.txt
    touch output/${meta.normal_id}.cpsr.grch38.html
    touch output/${meta.normal_id}.annotations.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
