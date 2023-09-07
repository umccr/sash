process BOLT_SMLV_SOMATIC_REPORT {
    tag "${meta.key}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.1.0-pcgr'

    input:
    tuple val(meta), path(smlv_vcf), path(smlv_unfiltered_vcf), path(purple_purity)
    path pcgr_data_dir
    path somatic_driver_panel_regions_coding
    path giab_regions
    path genome_fasta

    output:
    tuple val(meta), path('pcgr_output/')         , emit: pcgr_dir
    tuple val(meta), path('af_tumor.txt')         , emit: af_global
    tuple val(meta), path('af_tumor_keygenes.txt'), emit: af_keygenes
    tuple val(meta), path("*.bcftools_stats.txt") , emit: bcftools_stats
    tuple val(meta), path("*.variant_counts.yaml"), emit: variant_counts
    path 'versions.yml'                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt smlv_somatic report \\
        --tumor_name ${meta.tumor_id} \\
        --normal_name ${meta.normal_id} \\
        \\
        --vcf_fp ${smlv_vcf} \\
        --vcf_unfiltered_fp ${smlv_unfiltered_vcf} \\
        \\
        --pcgr_conda pcgr \\
        --pcgrr_conda pcgrr \\
        --pcgr_data_dir ${pcgr_data_dir} \\
        --purple_purity_fp ${purple_purity} \\
        \\
        --cancer_genes_fp ${somatic_driver_panel_regions_coding} \\
        --giab_regions_fp ${giab_regions} \\
        --genome_fp ${genome_fasta} \\
        \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pcgr_output/
    touch af_tumor.txt
    touch af_tumor_keygenes.txt
    touch ${meta.tumor_id}.somatic.variant_counts.yaml
    touch ${meta.tumor_id}.somatic.bcftools_stats.txt
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

