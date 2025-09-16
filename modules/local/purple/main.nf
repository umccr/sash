process PURPLE {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/hmftools-purple:4.2--hdfd78af_0'

    input:
    tuple val(meta), path(amber), path(cobalt), path(sv_tumor_vcf), path(sv_tumor_tbi), path(sv_normal_vcf), path(sv_normal_tbi), path(smlv_tumor_vcf), path(smlv_normal_vcf)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path gc_profile
    path sage_known_hotspots_somatic
    path sage_known_hotspots_germline
    path driver_gene_panel
    path ensembl_data_resources
    path germline_del
    path target_region_bed
    path target_region_ratios
    path target_region_msi_indels

    output:
    tuple val(meta), path('purple/'), emit: purple_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // allow custom heap fraction
    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''

    def sv_tumor_vcf_arg = sv_tumor_vcf ? "-somatic_sv_vcf ${sv_tumor_vcf}" : ''
    def sv_normal_vcf_arg = sv_normal_vcf ? "-germline_sv_vcf ${sv_normal_vcf}" : ''

    // NOTE(SW): use of 'smlv_tumor.vcf.gz' is intended here; see comment below in script block
    def smlv_tumor_vcf_arg = smlv_tumor_vcf ? "-somatic_vcf smlv_tumor.vcf.gz" : ''
    def smlv_normal_vcf_arg = smlv_normal_vcf ? "-germline_vcf ${smlv_normal_vcf}" : ''

    def sage_known_hotspots_germline_arg = sage_known_hotspots_germline ? "-germline_hotspots ${sage_known_hotspots_germline}" : ''
    def germline_del_arg = germline_del ? "-germline_del_freq_file ${germline_del}" : ''

    def target_region_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''
    def target_region_ratios_arg = target_region_ratios ? "-target_regions_ratios ${target_region_ratios}" : ''
    def target_region_msi_indels_arg = target_region_msi_indels ? "-target_regions_msi_indels ${target_region_msi_indels}" : ''

    """
    # NOTE(SW): input somatic SNVs must have `INFO/TIER` set to a valid value other than `LOW_CONFIDENCE` in order for
    # variants to be included in purity fitting, which is critical in some cases for an accurate estimate. The
    # `INFO/TIER` field is usually set by SAGE but for DRAGEN SNVs we must set it manually here. For sash SNVs it's
    # acceptable at this point in the workflow to assume all variants are high quality and so they are simply assigned
    # `TIER=HIGH_CONFIDENCE` to reach the minimal required change. For more information on SAGE tiering:
    #   - https://github.com/hartwigmedical/hmftools/tree/master/sage
    if [[ -n "${smlv_tumor_vcf}" ]]; then
        line_info=\$(bcftools view -h ${smlv_tumor_vcf} | grep -n '^##INFO' | cut -f1 -d: | tail -n1)
        header_entry='##INFO=<ID=TIER,Number=1,Type=String,Description="Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]">'
        {
            bcftools view -h ${smlv_tumor_vcf} | sed "\${line_info}a \${header_entry}";
            bcftools view -H ${smlv_tumor_vcf} | awk -F\$'\t' 'BEGIN { OFS="\t" } { \$8 = \$8 ";TIER=HIGH_CONFIDENCE"; print \$0 }';
        } | bcftools view -o smlv_tumor.vcf.gz
    fi;

    purple \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        ${reference_arg} \\
        -amber ${amber} \\
        -cobalt ${cobalt} \\
        ${sv_tumor_vcf_arg} \\
        ${sv_normal_vcf_arg} \\
        ${smlv_tumor_vcf_arg} \\
        ${smlv_normal_vcf_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -somatic_hotspots ${sage_known_hotspots_somatic} \\
        ${sage_known_hotspots_germline_arg} \\
        ${target_region_bed_arg} \\
        ${target_region_ratios_arg} \\
        ${target_region_msi_indels_arg} \\
        ${germline_del_arg} \\
        -gc_profile ${gc_profile} \\
        -circos \$(which circos) \\
        -threads ${task.cpus} \\
        -output_dir purple/

    # Remove the artificial `INFO/TIER` field in output PURPLE SNV file
    if [[ -n "${smlv_tumor_vcf}" ]]; then
        bcftools annotate -x INFO/TIER -o smlv_tumor.tmp.vcf.gz purple/${meta.tumor_id}.purple.somatic.vcf.gz
        mv smlv_tumor.tmp.vcf.gz purple/${meta.tumor_id}.purple.somatic.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: \$(purple -version | sed -n '/^Purple version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir purple/
    touch purple/${meta.tumor_id}.purple.cnv.gene.tsv
    touch purple/${meta.tumor_id}.purple.cnv.somatic.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.germline.tsv
    touch purple/${meta.tumor_id}.purple.driver.catalog.somatic.tsv
    touch purple/${meta.tumor_id}.purple.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.purity.tsv
    touch purple/${meta.tumor_id}.purple.qc
    touch purple/${meta.tumor_id}.purple.somatic.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.germline.vcf.gz
    touch purple/${meta.tumor_id}.purple.sv.vcf.gz
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
