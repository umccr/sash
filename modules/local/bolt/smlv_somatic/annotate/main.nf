process BOLT_SMLV_SOMATIC_ANNOTATE {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/umccr/bolt:0.2.14-pcgr'

    input:
    tuple val(meta), path(smlv_vcf)
    path somatic_driver_panel_regions_gene
    path annotations_dir
    path pon_dir
    path pcgr_data_dir

    output:
    tuple val(meta), path("output/${meta.tumor_id}.annotations.vcf.gz"), emit: vcf
    path 'versions.yml'                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt smlv_somatic annotate \\
        --tumor_name ${meta.tumor_id} \\
        --normal_name ${meta.normal_id} \\
        --vcf_fp ${smlv_vcf} \\
        --cancer_genes_fp ${somatic_driver_panel_regions_gene} \\
        --annotations_dir ${annotations_dir} \\
        --pon_dir ${pon_dir} \\
        --pcgr_data_dir ${pcgr_data_dir} \\
        --pcgr_conda pcgr \\
        --pcgrr_conda pcgrr \\
        --threads ${task.cpus} \\
        --output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/
    touch output/${meta.tumor_id}.annotations.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
