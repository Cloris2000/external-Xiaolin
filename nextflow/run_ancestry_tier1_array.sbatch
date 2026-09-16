#!/bin/bash
#SBATCH --job-name=ancestry_ct
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_sensitivity/ancestry_tier1_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_sensitivity/ancestry_tier1_%j.err

# Heterogeneity-focused ancestry screen, all 19 cell types on one whole node.
# After this job finishes: sbatch run_ancestry_tier1_merge.sbatch

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

OUT_DIR="${OUT_DIR:-${RESULTS_ROOT}/meta_sensitivity/ancestry_tier1}"
HIT_OVERLAP="${HIT_OVERLAP:-${RESULTS_ROOT}/meta_sensitivity/stratified_comparison_tier1/stratified_hit_overlap.tsv}"
UNIQUE_HITS="${UNIQUE_HITS:-${RESULTS_ROOT}/meta_sensitivity/stratified_comparison_tier1/stratified_unique_hits.tsv}"
COLOC_TSV="${COLOC_TSV:-${RESULTS_ROOT}/coloc/coloc_all_results.tsv}"
WORKER_LOGS="${LOG_ROOT}/meta_sensitivity/tier1_workers_${SLURM_JOB_ID:-manual}"
mkdir -p "${OUT_DIR}" "${WORKER_LOGS}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

run_one() {
    local ct="$1"
    python3 "${NF_DIR}/scripts/ancestry_tier1_screen.py" \
        --project-dir "${NF_DIR}" \
        --cell-type "${ct}" \
        --hit-overlap-tsv "${HIT_OVERLAP}" \
        --unique-hits-tsv "${UNIQUE_HITS}" \
        --coloc-tsv "${COLOC_TSV}" \
        --output-dir "${OUT_DIR}" \
        --i2-thresh 50 \
        > "${WORKER_LOGS}/${ct}.log" 2>&1 \
        && echo "  OK   ${ct}" \
        || { echo "  FAIL ${ct} (see ${WORKER_LOGS}/${ct}.log)"; return 1; }
}
export -f run_one
export NF_DIR OUT_DIR HIT_OVERLAP UNIQUE_HITS COLOC_TSV WORKER_LOGS

printf '%s\n' "${CELL_TYPES[@]}" \
  | xargs -P "${#CELL_TYPES[@]}" -I{} bash -c 'run_one "$@"' _ {} \
  && rc=0 || rc=$?

echo ""
if [ "$rc" -ne 0 ]; then
    echo "ERROR: at least one cell type failed. Logs: ${WORKER_LOGS}"
    exit 1
fi
echo "Ancestry tier1 screen complete: ${OUT_DIR}"
