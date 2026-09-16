#!/bin/bash
# Preflight for the cohort_gwas_v2 re-run.  Login node, read-only on inputs.
#
# For every cohort in validation/cohort_trillium_paths.tsv:
#   - config present; per-chromosome genotype input for chr1-22 resolves and is non-empty
#   - RNA count matrix, metadata, biospecimen, clinical metadata present
#   - samples_to_keep resolved (repo copy, else SCC copy, else none)
#   - SCC reference outputs present for the later comparison
#   - NABEC chr1 is the recovered file (size + BGZF EOF)
# Writes ${V2_PREFLIGHT}/manifest.tsv and md5s of the small inputs, and exits 1
# if anything required is missing.
#
#   scripts/preflight_cohort_gwas_v2.sh            # all 15
#   COHORTS="NABEC GVEX" scripts/preflight_cohort_gwas_v2.sh

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
source "${NF_DIR}/scripts/cohort_gwas_v2_lib.sh"

mkdir -p "${V2_PREFLIGHT}"
MANIFEST="${V2_PREFLIGHT}/manifest.tsv"
MD5S="${V2_PREFLIGHT}/small_inputs.md5"
: > "${MD5S}"
printf "cohort\tconfig\tqc_study\tgenotype_source\tn_chr_ok\tgenotype_GB\tcount_matrix\tmetadata\tbiospecimen\tclinical_metadata\tsamples_to_keep\tscc_raw_p\tscc_refs\tstatus\n" > "${MANIFEST}"

filter="${COHORTS:-}"; filter="${filter//,/ }"
fail=0

check_file() {  # check_file <path-or-dash> -> "OK:<bytes>" | "-" | "MISSING"
    local p="$1"
    [ "$p" = "-" ] || [ -z "$p" ] && { echo "-"; return; }
    [ -s "$p" ] && echo "OK:$(stat -c %s "$p")" || echo "MISSING"
}

while IFS=$'\t' read -r cohort config count_matrix metadata biospecimen clinical vcf_dir vcf_glob scc_results; do
    if [ -n "$filter" ]; then case " $filter " in *" $cohort "*) ;; *) continue ;; esac; fi
    cfg="${NF_DIR}/${config}"
    status=OK
    if [ ! -s "$cfg" ]; then
        echo "FAIL ${cohort}: config missing ${cfg}"; fail=1
        printf "%s\t%s\tNA\tNA\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tCONFIG_MISSING\n" "$cohort" "$config" >> "${MANIFEST}"
        continue
    fi
    study=$(qc_study_name "$cfg")

    # genotype inputs, chr1-22
    n_ok=0; bytes=0; first_missing=""
    for c in $(seq 1 22); do
        f=$(resolve_vcf "$cfg" "$c" "$cohort")
        if [ -s "$f" ]; then n_ok=$((n_ok+1)); bytes=$((bytes + $(stat -c %s "$f"))); else [ -z "$first_missing" ] && first_missing="$f"; fi
    done
    src_dir=$(dirname "$(resolve_vcf "$cfg" 1 "$cohort")")
    [ "$n_ok" -ne 22 ] && { status=FAIL; echo "FAIL ${cohort}: ${n_ok}/22 genotype files; first missing: ${first_missing}"; }

    # RNA / metadata
    cm=$(check_file "$count_matrix"); md=$(check_file "$metadata"); bio=$(check_file "$biospecimen"); cl=$(check_file "$clinical")
    for v in "$cm" "$md"; do [ "$v" = "MISSING" ] && status=FAIL; done
    [ "$bio" = "MISSING" ] && status=FAIL
    [ "$cl" = "MISSING" ] && status=FAIL
    for p in "$metadata" "$biospecimen" "$clinical"; do [ -s "$p" ] && md5sum "$p" >> "${MD5S}"; done

    # samples_to_keep
    stk=$(resolve_samples_to_keep "$cohort" "$cfg")
    if [ -n "$stk" ]; then stk_str="$stk ($(( $(wc -l < "$stk") - 1 )) samples)"; md5sum "$stk" >> "${MD5S}";
    elif [ -n "$(cfg_param "$cfg" samples_to_keep)" ]; then stk_str="CONFIG_REFERENCES_MISSING_FILE(all phenotyped samples kept)";
    else stk_str="-"; fi

    # SCC reference for comparison
    n_raw=$(ls "${SCC_RES}/${scc_results}"/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
    refs=""
    for r in cell_proportions.csv phenotypes_RINT.txt pca.csv; do [ -s "${SCC_RES}/${scc_results}/${r}" ] && refs="${refs}${r%%.*}," || refs="${refs}NO_${r%%.*},"; done
    ls "${SCC_RES}/${scc_results}"/*.QC.final.psam >/dev/null 2>&1 && refs="${refs}psam" || refs="${refs}NO_psam"
    [ "$n_raw" -ne 19 ] && echo "WARN ${cohort}: SCC has ${n_raw}/19 raw_p (comparison will be partial)"

    # NABEC recovered chr1
    if [ "$cohort" = "NABEC" ]; then
        f1=$(resolve_vcf "$cfg" 1 "$cohort"); sz=$(stat -c %s "$f1")
        if [ "$sz" -lt 400000000 ] || ! bgzf_eof_ok "$f1"; then status=FAIL; echo "FAIL NABEC: chr1 ${f1} is ${sz} bytes / BGZF EOF $(bgzf_eof_ok "$f1" && echo ok || echo MISSING) - not the recovered file"; else echo "OK   NABEC chr1 recovered file: ${sz} bytes, BGZF EOF present"; fi
    fi

    [ "$status" = FAIL ] && fail=1
    printf "%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$cohort" "$config" "$study" "$src_dir" "$n_ok" "$(echo "$bytes/1073741824" | bc -l)" \
        "$cm" "$md" "$bio" "$cl" "$stk_str" "$n_raw" "$refs" "$status" >> "${MANIFEST}"
    printf "%-5s %-17s geno=%2s/22 (%5.1f GB)  rna=%s meta=%s bio=%s clin=%s  keep=%s\n" "$status" "$cohort" "$n_ok" "$(echo "$bytes/1073741824" | bc -l)" "${cm%%:*}" "${md%%:*}" "${bio%%:*}" "${cl%%:*}" "${stk_str%% *}"
done < <(map_rows)

echo
echo "Manifest: ${MANIFEST}"
echo "md5s    : ${MD5S}"
[ "$fail" -eq 0 ] && echo "PREFLIGHT PASSED" || echo "PREFLIGHT FAILED (see FAIL lines above)"
exit "$fail"
