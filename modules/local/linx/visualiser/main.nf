process LINX_VISUALISER {
    tag "${meta.id}"
    label 'process_medium'

    container 'ghcr.io/umccr/linx:2.1'

    input:
    tuple val(meta), path(linx_annotation_dir)
    val genome_ver
    path ensembl_data_resources

    output:
    tuple val(meta), path('plots/'), emit: plots
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    # NOTE(SW): the output plot directories are always required for ORANGE, which is straightfoward to handle with POSIX
    # fs but more involved with FusionFS since it will not write empty directories to S3. A placeholder file can't be
    # used in the plot directory to force FusionFS to create the directory as ORANGE will treat the placeholder as a PNG
    # and fail. Optional outputs are possible but requires further channel logic and output to detect when complete.
    # Instead I place the two plot output directories under a parent directory, only operating on that to allow use of a
    # placeholder and support empty outputs when using FusionFS. Handling missing/non-existent directories are deferred
    # to downstream processes, bypassing the need to implement further channel operations.

    mkdir -p plots/

    # NOTE(SW): LINX v1.24.1 require trailing slashes for the -plot_out and -data_out arguments since no filesystem
    # separator is used when constructing fusion plot output filepaths.

    # https://github.com/hartwigmedical/hmftools/blob/linx-v1.24.1/linx/src/main/java/com/hartwig/hmftools/linx/visualiser/circos/ChromosomeRangeExecution.java#L22-L29
    # https://github.com/hartwigmedical/hmftools/blob/linx-v1.24.1/linx/src/main/java/com/hartwig/hmftools/linx/visualiser/circos/FusionExecution.java#L18-L23

    # Generate all chromosome and cluster plots by default

    linx \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -vis_file_dir ${linx_annotation_dir} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -circos \$(which circos) \\
        -threads ${task.cpus} \\
        -plot_out plots/all/ \\
        -data_out data/all/

    # Create placeholders to force FusionFS to create parent plot directory on S3
    if [[ \$(ls plots/ | wc -l) -eq 0 ]]; then
        touch plots/.keep;
    fi;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: \$(linx -version | sed -n '/^Linx version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p plots/{all,reportable}/
    touch plots/{all,reportable}/placeholder

    echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
    """
}
