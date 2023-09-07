process BOLT_SV_SOMATIC_ANNOTATE {
    tag "${meta.key}"
    label 'process_low'

    container 'docker.io/scwatts/bolt:0.1.0-snpeff'

    input:
    tuple val(meta), path(sv_vcf), path(cnv_tsv)
    path genome_fasta
    path genome_fai
    path snpeff_database

    output:
    tuple val(meta), path("${meta.id}*annotated.vcf.gz"), emit: vcf
    path 'versions.yml'                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bolt sv_somatic annotate \\
        --tumor_name ${meta.id} \\
        --sv_fp ${sv_vcf} \\
        --cnv_fp ${cnv_tsv} \\
        --reference_fasta_fp ${genome_fasta} \\
        --snpeff_database_dir ${snpeff_database}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bolt: \$(bolt --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.annotated.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
