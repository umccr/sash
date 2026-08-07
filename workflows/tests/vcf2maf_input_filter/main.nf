// Test-only workflow: replicates the VCF2MAF input channel chain from
// workflows/sash.nf (the restoreMeta + .map + .filter block downstream of
// BOLT_SMLV_SOMATIC_REPORT.out.pcgr_pass_vcf) so it can be exercised in
// isolation, since this logic is inline in sash.nf and not otherwise
// importable. See modules/local/bolt/smlv_somatic/report/tests/main.nf.test
// for the empirical confirmation that an optional:true output that wasn't
// produced yields zero channel items for that sample (not a null entry).

workflow VCF2MAF_INPUT_FILTER {
    take:
    ch_pcgr_pass_vcf // channel: [ meta_bolt, vcf ] — only emitted for samples where PCGR ran
    ch_inputs        // channel: [ meta ] — one entry per sample, regardless of PCGR skip

    main:
    ch_out = WorkflowSash.restoreMeta(ch_pcgr_pass_vcf, ch_inputs)
        .map { meta, vcf ->
            def meta_vcf2maf = [
                key: meta.id,
                id: meta.id,
                tumor_id: meta.tumor_id,
                normal_id: meta.normal_id,
            ]
            return [meta_vcf2maf, vcf]
        }
        .filter { _meta, vcf -> vcf }

    emit:
    out = ch_out
}
