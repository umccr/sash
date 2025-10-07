process BOLT_OTHER_CANCER_REPORT {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/umccr/bolt:0.2.17-gpgr'

    input:
    tuple val(meta), path(smlv_somatic_vcf), path(smlv_somatic_bcftools_stats), path(smlv_somatic_counts_process), path(sv_somatic_tsv), path(sv_somatic_vcf), path(cnv_somatic_tsv), path(af_global), path(af_keygenes), path(purple_baf_plot), path(purple_dir), path(virusbreakend_dir), path(dragen_hrd, stageAs: 'dragen_hrd/*')
    path somatic_driver_panel
    path oncokb_genes

    output:
    path 'output/'             , emit: cancer_report_dir
    path '*.cancer_report.html', emit: cancer_report
    path 'versions.yml'        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def dragen_hrd_arg = dragen_hrd.name != 'dragen_hrd' ? "--dragen_hrd_fp \$(pwd)/dragen_hrd/${dragen_hrd.name}" : ''

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
        --smlv_somatic_bcftools_stats_fp \$(pwd)/${smlv_somatic_bcftools_stats} \\
        --smlv_somatic_counts_process_fp \$(pwd)/${smlv_somatic_counts_process} \\
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
        ${dragen_hrd_arg} \\
        \\
        --cancer_genes_fp \$(pwd)/${somatic_driver_panel} \\
        --oncokb_genes_fp \$(pwd)/${oncokb_genes} \\
        \\
        --output_dir \$(pwd)/output/

    mv output/${meta.tumor_id}.cancer_report.html ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/
    touch output/cancer_report.stub
    touch ${meta.tumor_id}.cancer_report.html
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
