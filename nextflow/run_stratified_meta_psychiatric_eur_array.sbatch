#!/bin/bash
#SBATCH --job-name=meta_psych_EUR_ct
#SBATCH --partition=medium
#SBATCH --array=0-18
#SBATCH --time=0:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/meta_sensitivity/stratified_psychiatric_eur_array_%A_%a.out
#SBATCH --error=logs/meta_sensitivity/stratified_psychiatric_eur_array_%A_%a.err

# METAL-only: reuse pre-harmonized files from EUR_homogeneous meta (no Nextflow).
# After all 19 tasks complete, run: sbatch run_stratified_meta_psychiatric_eur_merge.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
METAL_PATH=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/METAL/generic-metal/executables
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

CELL_TYPE="${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}"
echo "psychiatric_EUR METAL — cell type: ${CELL_TYPE} (array task ${SLURM_ARRAY_TASK_ID})"

python3 "${NF_DIR}/scripts/run_metal_from_harmonized.py" \
    --cell-type "${CELL_TYPE}" \
    --cohorts CMC_MSSM CMC_PENN CMC_PITT \
    --harmonized-dir "${NF_DIR}/results/meta_sensitivity/stratified_eur_homogeneous_15cohorts/harmonized" \
    --output-dir "${NF_DIR}/results/meta_sensitivity/stratified_psychiatric_eur_15cohorts" \
    --metal-path "${METAL_PATH}"

echo "Done: ${CELL_TYPE}"
