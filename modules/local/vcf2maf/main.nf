process VCF2MAF {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/qclayssen/vcf2maf:v1.6.22'

    input:
    tuple val(meta), path(vcf)
    path genome_fasta
    path vep_dir


    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ncbi_build = params.genome.build == 'hg38' ? 'GRCh38' : params.genome.build == 'GRCh37' ? 'GRCh37' : 'GRCh38'
    def vep_dir_cmd = vep_dir ? "--vep-data ${vep_dir}" : ""
    def uncompressed_vcf = "${meta.id}-temp.vcf"

    """
    # Handle VEP settings like nf-core module
    if [ "${params.vep_dir}" ]; then
        VEP_CMD="--vep-path \$(dirname \$(type -p vep))"
    else
        VEP_CMD=""
    fi

    # Uncompress VCF if needed (custom enhancement for compressed files)
    if [[ ${vcf} == *.gz ]]; then
        gunzip -c ${vcf} > ${uncompressed_vcf}
    else
        cp ${vcf} ${uncompressed_vcf}
    fi

    vcf2maf.pl \\
        ${args} \\
        \$VEP_CMD \\
        ${vep_dir_cmd} \\
        --input-vcf ${uncompressed_vcf} \\
        --output-maf ${prefix}.maf \\
        --ref-fasta ${genome_fasta} \\
        --tumor-id ${meta.tumor_id} \\
        --normal-id ${meta.normal_id} \\
        --ncbi-build ${ncbi_build} \\
        2> >(grep -v "Use of uninitialized value" >&2)

    # Clean up temporary file
    rm -f ${uncompressed_vcf}

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
