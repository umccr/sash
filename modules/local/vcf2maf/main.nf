process VCF2MAF {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/qclayssen/vcf2maf:debian_v1.6.22'

    input:
    tuple val(meta), path(vcf)
    path genome_fasta
    val genome_build


    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def uncompressed_vcf = "${meta.id}-temp.vcf"

    """
    gunzip -c ${vcf} > ${uncompressed_vcf}

    vcf2maf \\
        --inhibit-vep \\
        --input-vcf ${uncompressed_vcf} \\
        --output-maf ${meta.id}.maf \\
        --ref-fasta ${genome_fasta} \\
        --tumor-id ${meta.tumor_id} \\
        --normal-id ${meta.normal_id} \\
        --ncbi-build ${genome_build} \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: 1.6.22
    END_VERSIONS
    """
}
