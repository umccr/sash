#!/usr/bin/env bash
# Regression test for the MNVTAG header-detection guard in ../main.nf (lines ~32-38).
# Runs the guard logic for real (not stubbed) against the nf-test fixtures, and against a
# synthetic large header, to catch the SIGPIPE/pipefail race fixed for umccr/sash#66:
# `bcftools view -h | grep -q` can make bcftools receive SIGPIPE mid-write on large real
# headers, and under `set -o pipefail` that flips the `if` to the wrong branch even though
# grep genuinely found (or didn't find) the pattern. nf-test only exercises `-stub`, which
# skips the script body entirely, so this guard logic needs its own coverage.
#
# NOTE: the guard snippet below must be kept identical to main.nf's script block.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "${tmp_dir}"' EXIT

run_guard() {
    local vcf="$1"
    local workdir="$2"
    (
        cd "${workdir}"
        bcftools view -h "${vcf}" > header.txt
        if grep -q '##INFO=<ID=MNVTAG,' header.txt; then
            echo "filtered"
        else
            echo "unfiltered"
        fi
    )
}

fail=0

with_mnvtag_vcf="${script_dir}/fixtures/with_mnvtag/TUMOR.hard-filtered.vcf.gz"
without_mnvtag_vcf="${script_dir}/fixtures/without_mnvtag/TUMOR.hard-filtered.vcf.gz"

work1="${tmp_dir}/with_mnvtag"
mkdir -p "${work1}"
result=$(run_guard "${with_mnvtag_vcf}" "${work1}")
if [[ "${result}" != "filtered" ]]; then
    echo "FAIL: with_mnvtag fixture expected 'filtered', got '${result}'" >&2
    fail=1
else
    echo "PASS: with_mnvtag fixture -> filtered"
fi

work2="${tmp_dir}/without_mnvtag"
mkdir -p "${work2}"
result=$(run_guard "${without_mnvtag_vcf}" "${work2}")
if [[ "${result}" != "unfiltered" ]]; then
    echo "FAIL: without_mnvtag fixture expected 'unfiltered', got '${result}'" >&2
    fail=1
else
    echo "PASS: without_mnvtag fixture -> unfiltered"
fi

# Reproduce the original SIGPIPE race: a large header (thousands of ##contig lines, like a
# real hg38 VCF) with MNVTAG declared early, well before the pipe buffer fills up.
large_header_vcf="${tmp_dir}/large_header.vcf"
{
    echo '##fileformat=VCFv4.2'
    echo '##INFO=<ID=MNVTAG,Number=1,Type=String,Description="MNV tag">'
    for i in $(seq 1 20000); do
        echo "##contig=<ID=chr1_synthetic_${i},length=1000000>"
    done
    echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR'
    echo -e 'chr1\t100\t.\tA\tT\t60\tPASS\tDP=30;MNVTAG=.\tGT\t0/1'
} > "${large_header_vcf}"
bgzip -f "${large_header_vcf}"
tabix -p vcf "${large_header_vcf}.gz"

work3="${tmp_dir}/large_header"
mkdir -p "${work3}"
result=$(run_guard "${large_header_vcf}.gz" "${work3}")
if [[ "${result}" != "filtered" ]]; then
    echo "FAIL: large synthetic header expected 'filtered' (SIGPIPE/pipefail race), got '${result}'" >&2
    fail=1
else
    echo "PASS: large synthetic header -> filtered (no SIGPIPE race)"
fi

exit "${fail}"
