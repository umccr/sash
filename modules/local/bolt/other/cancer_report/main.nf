process BOLT_OTHER_CANCER_REPORT {
    tag "${meta.key}"
    label 'process_low'

    container 'docker.io/scwatts/bolt_gpgr:0.1.0'

    input:
    tuple val(meta), path(smlv_somatic_vcf), path(sv_somatic_tsv), path(sv_somatic_vcf), path(cnv_somatic_tsv), path(af_global), path(af_keygenes), path(purple_baf_plot), path(purple_dir), path(virusbreakend_dir)
    path somatic_driver_panel

    output:
    tuple val(meta), path('cancer_report/'), emit: cancer_report_dir
    path 'versions.yml'                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    # NOTE(SW): gpgr requires aboslute paths

    bolt other cancer_report \\
        --subject_name ${meta.subject_id} \\
        --tumor_name ${meta.tumor_id} \\
        \\
        --af_global_fp \$(pwd)/${af_global} \\
        --af_keygenes_fp \$(pwd)/${af_keygenes} \\
        \\
        --smlv_somatic_vcf_fp \$(pwd)/${smlv_somatic_vcf} \\
        \\
        --sv_somatic_tsv_fp \$(pwd)/${sv_somatic_tsv} \\
        --sv_somatic_vcf_fp \$(pwd)/${sv_somatic_vcf} \\
        --cnv_somatic_tsv_fp \$(pwd)/${cnv_somatic_tsv} \\
        \\
        --purple_baf_plot_fp \$(pwd)/${purple_baf_plot} \\
        \\
        --purple_dir \$(pwd)/${purple_dir} \\
        --virusbreakend_dir \$(pwd)/${virusbreakend_dir} \\
        \\
        --cancer_genes_fp \$(pwd)/${somatic_driver_panel} \\
        \\
        --output_dir \$(pwd)/cancer_report/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cancer_report/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
