#!/bin/bash
#SBATCH --job-name=cmp_meta_v2_scc
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/cmp_meta_v2_scc_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/cmp_meta_v2_scc_%j.err

# Compare the corrected meta (v2) with the SCC manuscript meta for all 19 cell types:
# Miami plot (v2 above the axis, SCC mirrored below), per-locus lead table, summary row.
#
# SCC side reads METAL's <prefix>1.tbl -- the table behind the manuscript figures.
# The plain <prefix>.tbl in that tree has 14.8M rows vs 8.49M and is a pre-filter
# artifact of the old OUTFILE handling (fixed in c1efce9).
#
# ~2 min and ~6 GB per cell type; 8 run concurrently on the node.
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

NEW_DIR="${NEW_DIR:-${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2}"
SCC_META="${SCC_META:-${SCC_DIR}/results/meta_analysis_15cohorts}"
OUT_DIR="${OUT_DIR:-${NEW_DIR}/compare_manuscript}"
CTS=(Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP L5.ET L6.CT L6b LAMP5
     Microglia OPC Oligodendrocyte PAX6 PVALB Pericyte SST VIP VLMC)
mkdir -p "${OUT_DIR}" "${LOG_ROOT}/meta_15cohorts/cmp_workers"

echo "Comparing v2 vs SCC manuscript meta, $(date -Is)"
echo "  new: ${NEW_DIR}"
echo "  scc: ${SCC_META}  (using *1.tbl)"
echo "  out: ${OUT_DIR}"

pids=(); running=0
for ct in "${CTS[@]}"; do
    "${RSCRIPT}" "${NF_DIR}/scripts/compare_meta_v2_vs_scc.R" \
        --cell_type "${ct}" --new_dir "${NEW_DIR}" --scc_dir "${SCC_META}" --out_dir "${OUT_DIR}" \
        > "${LOG_ROOT}/meta_15cohorts/cmp_workers/${ct}.log" 2>&1 &
    pids+=($!); running=$((running+1))
    if [ "${running}" -ge 8 ]; then wait -n || true; running=$((running-1)); fi
done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done

# Merge the per-cell-type summaries into one table.
first=$(ls "${OUT_DIR}"/*_summary.tsv 2>/dev/null | head -1)
if [ -n "${first}" ]; then
    { head -1 "${first}"; for f in "${OUT_DIR}"/*_summary.tsv; do tail -n +2 "$f"; done; } > "${OUT_DIR}/ALL_summary.tsv"
    echo ""; echo "=== v2 vs SCC manuscript meta, all cell types ==="
    column -t "${OUT_DIR}/ALL_summary.tsv"
fi
echo ""
echo "Miami plots : ${OUT_DIR}/<cell_type>_miami.png  ($(ls "${OUT_DIR}"/*_miami.png 2>/dev/null | wc -l)/19)"
echo "Lead tables : ${OUT_DIR}/<cell_type>_leads.tsv"
exit "${rc}"
