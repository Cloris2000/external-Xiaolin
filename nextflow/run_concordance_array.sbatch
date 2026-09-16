#!/bin/bash
#SBATCH --job-name=concordance
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/zhoux156/logs/concordance/concordance_%j.out
#SBATCH --error=/scratch/zhoux156/logs/concordance/concordance_%j.err

# Bulk vs sn direction concordance, all 19 cell types on one whole node.
# After the job finishes, plot from a login node:
#   Rscript scripts/plot_concordance_barplot.R

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

WORKER_LOGS="${LOG_ROOT}/concordance/workers_${SLURM_JOB_ID:-manual}"
mkdir -p "${WORKER_LOGS}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia Oligodendrocyte
  OPC PAX6 Pericyte PVALB SST VIP VLMC
)

run_one() {
    local ct="$1"
    "${RSCRIPT}" "${NF_DIR}/scripts/bulk_sn_concordance_per_ct.R" \
        --celltype "${ct}" \
        > "${WORKER_LOGS}/${ct}.log" 2>&1 \
        && echo "  OK   ${ct}" \
        || { echo "  FAIL ${ct} (see ${WORKER_LOGS}/${ct}.log)"; return 1; }
}
export -f run_one
export NF_DIR RSCRIPT WORKER_LOGS

printf '%s\n' "${CELL_TYPES[@]}" \
  | xargs -P "${#CELL_TYPES[@]}" -I{} bash -c 'run_one "$@"' _ {} \
  && rc=0 || rc=$?

echo ""
if [ "$rc" -ne 0 ]; then
    echo "ERROR: at least one cell type failed. Logs: ${WORKER_LOGS}"
    exit 1
fi
echo "Concordance complete."
