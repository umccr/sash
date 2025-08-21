process VCF2MAF {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf2maf:1.6.22--hdfd78af_0':
        'biocontainers/vcf2maf:1.6.22--hdfd78af_0' }"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor_id = meta.tumor_id ?: meta.id
    def normal_id = meta.normal_id ?: "${meta.id}_normal"
    def ncbi_build = params.genome.build == 'hg38' ? 'GRCh38' : params.genome.build == 'GRCh37' ? 'GRCh37' : 'GRCh38'
    def vep_cache_cmd = params.vep_cache ? "--vep-data ${params.vep_cache}" : ""
    def uncompressed_vcf = "${prefix}-temp.vcf"
    
    """
    # Handle VEP settings like nf-core module
    if [ "${params.vep_cache}" ]; then
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
        ${vep_cache_cmd} \\
        --input-vcf ${uncompressed_vcf} \\
        --output-maf ${prefix}.maf \\
        --ref-fasta ${genome_fasta} \\
        --tumor-id ${tumor_id} \\
        --normal-id ${normal_id} \\
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
