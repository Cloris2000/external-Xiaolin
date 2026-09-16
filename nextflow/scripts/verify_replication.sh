#!/usr/bin/env bash
# Compare the Trillium 15-cohort meta-analysis against the CAMH/SCC run.
#
# The SCC results directory holds TWO tables per cell type, and picking the right
# reference matters:
#
#   <ct>_meta_analysis_<suffix>1.tbl   what METAL actually produced on its last
#                                      run (post-liftover, ~946 MB)
#   <ct>_meta_analysis_<suffix>.tbl    what got published and fed downstream
#                                      (pre-liftover May 21 leftover, ~1.63 GB)
#
# METAL writes <prefix>1.tbl, and the old module pointed OUTFILE at the shared
# results directory and then picked the largest match with `ls -S`, so the stale
# larger table kept winning. The correct replication target is therefore the
# *1.tbl file, not the published one.
#
#   scripts/verify_replication.sh
#   TRI_DIR=/scratch/$USER/results/meta_analysis_15cohorts scripts/verify_replication.sh

set -uo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

TRI_DIR="${TRI_DIR:-${RESULTS_ROOT}/meta_analysis_15cohorts}"
SCC_META_DIR="${SCC_META_DIR:-${SCC_DIR}/results/meta_analysis_15cohorts}"
SUFFIX="${SUFFIX:-AMP_AD_Mayo_AMP_AD_Rush_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_v10_GVEX_MSBB_Mayo_NABEC_NIMH_HBCC_1M_NIMH_HBCC_Omni5M_NIMH_HBCC_h650_ROSMAP_ROSMAP_array}"
# metal     = SCC *1.tbl (what METAL last wrote: post-liftover, ~946 MB)
# published = SCC .tbl   (what downstream consumed: May 21 pre-liftover, ~1.63 GB)
REF_KIND="${REF_KIND:-metal}"
JOBS="${JOBS:-19}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

case "${REF_KIND}" in
    metal)     REF_SUFFIX="1.tbl" ; REF_LABEL="SCC *1.tbl (METAL last wrote, post-liftover)" ;;
    published) REF_SUFFIX=".tbl"  ; REF_LABEL="SCC .tbl (published / May 21 pre-liftover)" ;;
    *) echo "ERROR: REF_KIND must be metal or published, got '${REF_KIND}'" >&2; exit 2 ;;
esac

echo "Trillium : ${TRI_DIR}"
echo "Reference: ${SCC_META_DIR} (${REF_LABEL})"
echo ""

compare_one() {
    local ct="$1"
    local tri="${TRI_DIR}/${ct}_meta_analysis_${SUFFIX}.tbl"
    local ref="${SCC_META_DIR}/${ct}_meta_analysis_${SUFFIX}${REF_SUFFIX}"

    [ -s "$tri" ] || { printf '%-16s MISSING_TRILLIUM\n' "$ct"; return 2; }
    [ -s "$ref" ] || { printf '%-16s MISSING_REFERENCE\n' "$ct"; return 2; }

    local a b
    a=$(md5sum "$tri" | cut -d' ' -f1)
    b=$(md5sum "$ref" | cut -d' ' -f1)
    if [ "$a" = "$b" ]; then
        printf '%-16s MATCH    %s  %s bytes\n' "$ct" "${a:0:12}" "$(stat -c%s "$tri")"
    else
        printf '%-16s DIFFER   trillium=%s(%s B) scc=%s(%s B)\n' \
            "$ct" "${a:0:12}" "$(stat -c%s "$tri")" "${b:0:12}" "$(stat -c%s "$ref")"
        return 1
    fi
}
export -f compare_one
export TRI_DIR SCC_META_DIR SUFFIX REF_SUFFIX

printf '%s\n' "${CELL_TYPES[@]}" \
  | xargs -P "${JOBS}" -I{} bash -c 'compare_one "$@"' _ {} \
  | sort | tee /tmp/verify_replication.$$

echo ""
n_match=$(grep -c ' MATCH ' /tmp/verify_replication.$$ || true)
n_diff=$(grep -c ' DIFFER ' /tmp/verify_replication.$$ || true)
n_miss=$(grep -c ' MISSING' /tmp/verify_replication.$$ || true)
echo "match=${n_match}/${#CELL_TYPES[@]}  differ=${n_diff}  missing=${n_miss}"
rm -f /tmp/verify_replication.$$

[ "${n_match}" -eq "${#CELL_TYPES[@]}" ] || exit 1
echo "Replication verified: every cell type is byte-identical to the CAMH run."
