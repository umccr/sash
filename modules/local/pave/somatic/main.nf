process PAVE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/umccr/pave:1.8'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path genome_fasta
    val genome_ver
    path genome_fai
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    pave \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -sample ${meta.sample_id} \\
        -vcf_file ${vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -clinvar_vcf ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -read_pass_only \\
        -threads ${task.cpus} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.sage.pave_somatic.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
