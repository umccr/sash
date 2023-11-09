process BOLT_OTHER_PURPLE_BAF_PLOT {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.2.5-circos'

    input:
    tuple val(meta), path(purple_dir)
    path circos_config
    path circos_gaps

    output:
    tuple val(meta), path("output/*.png"), emit: plot
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    bolt other purple_baf_plot \
        --tumor_name ${meta.tumor_id} \
        --purple_dir ${purple_dir} \
        --circos_conf_fp ${circos_config} \
        --circos_gaps_fp ${circos_gaps} \
        --output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p purple_baf/
    touch purple_baf/${meta.tumor_id}.circos_baf.png
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
