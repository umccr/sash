process VCF2MAF {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/vcf2maf:1.6.22--hdfd78af_0'

    input:
    tuple val(meta), path(vcf)
    path genome_fasta


    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gunzip -c ${vcf} > ${meta.id}-temp.vcf

    vcf2maf.pl \\
        --inhibit-vep \\
        --input-vcf ${meta.id}-temp.vcf \\
        --output-maf ${meta.id}.maf \\
        --ref-fasta ${genome_fasta} \\
        --tumor-id ${meta.tumor_id} \\
        --normal-id ${meta.normal_id} \\
        --ncbi-build "GRCh38" \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: \$(vcf2maf.pl --help | grep -o 'vcf2maf [0-9.]*' | sed 's/vcf2maf //' || echo "1.6.22")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.maf

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    END_VERSIONS
    """
}
