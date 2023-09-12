process BOLT_OTHER_MULTIQC_REPORT {
    tag "${meta.key}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.1.2-multiqc'

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path('multiqc_report/'), emit: multiqc_report_dir
    path 'versions.yml'                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    bolt other multiqc_report \\
        --tumor_name ${meta.tumor_id} \\
        --normal_name ${meta.normal_id} \\
        --input_dir ./ \\
        --output_dir multiqc_report/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p multiqc_report/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
