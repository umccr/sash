process BOLT_SV_SOMATIC_PRIORITISE {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.2.4'

    input:
    tuple val(meta), path(sv_vcf)
    path known_fusion_pairs
    path known_fusion_heads
    path known_fusion_tails
    path fusioncatcher_pairs
    path somatic_driver_panel_genes
    path somatic_driver_panel_genes_ts
    path appris

    output:
    tuple val(meta), path("output/*.sv.*.tsv")   , emit: sv_tsv
    tuple val(meta), path("output/*.sv.*.vcf.gz"), emit: sv_vcf
    tuple val(meta), path("output/*.cnv.*.tsv")  , emit: cnv_tsv
    path 'versions.yml'                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt sv_somatic prioritise \\
        --tumor_name ${meta.tumor_id} \\
        --sv_vcf ${sv_vcf} \\
        --refdata_known_fusion_pairs ${known_fusion_pairs} \\
        --refdata_known_fusion_heads ${known_fusion_heads} \\
        --refdata_known_fusion_tails ${known_fusion_tails} \\
        --refdata_fusioncatcher_pairs ${fusioncatcher_pairs} \\
        --refdata_key_genes ${somatic_driver_panel_genes} \\
        --refdata_key_tsgenes ${somatic_driver_panel_genes_ts} \\
        --appris_fp ${appris} \\
        --output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sv_prioritise/
    touch sv_prioritise/${meta.tumor_id}.prioritised.sv.placeholder.tsv
    touch sv_prioritise/${meta.tumor_id}.prioritised.sv.placeholder.vcf.gz
    touch sv_prioritise/${meta.tumor_id}.prioritised.cnv.placeholder.tsv
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
