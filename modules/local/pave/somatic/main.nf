process PAVE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/hmftools-pave:1.8--hdfd78af_0'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path genome_fasta
    val genome_ver
    path genome_fai
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    bcftools view --exclude 'INFO/MNVTAG!="."' --write-index=tbi --output ${meta.sample_id}.mnv_filtred.vcf.gz ${vcf}

    pave \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_vcf ${meta.sample_id}.mnv_filtred.vcf.gz \\
        -output_vcf ${meta.sample_id}.pave.somatic.vcf.gz \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -gnomad_freq_dir ${gnomad_resource} \\
        -clinvar_vcf ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -mappability_bed ${segment_mappability} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -write_pass_only \\
        -threads ${task.cpus} \\
        -output_dir ./

    rm ${meta.sample_id}.mnv_filtred.vcf.gz{,.tbi}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed -n '/^Pave version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.pave.somatic.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
