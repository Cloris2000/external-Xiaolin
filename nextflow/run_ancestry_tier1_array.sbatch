#!/bin/bash
#SBATCH --job-name=ancestry_ct
#SBATCH --partition=medium
#SBATCH --array=0-18
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/meta_sensitivity/ancestry_tier1_array_%A_%a.out
#SBATCH --error=logs/meta_sensitivity/ancestry_tier1_array_%A_%a.err

# One cell type per array task — heterogeneity-focused ancestry screen.
# After all 19 tasks complete, run: sbatch run_ancestry_tier1_merge.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

CELL_TYPE="${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing cell type: ${CELL_TYPE} (array task ${SLURM_ARRAY_TASK_ID})"

python3 "${NF_DIR}/scripts/ancestry_tier1_screen.py" \
    --project-dir "${NF_DIR}" \
    --cell-type "${CELL_TYPE}" \
    --hit-overlap-tsv "results/meta_sensitivity/stratified_comparison_tier1/stratified_hit_overlap.tsv" \
    --unique-hits-tsv "results/meta_sensitivity/stratified_comparison_tier1/stratified_unique_hits.tsv" \
    --coloc-tsv "results/coloc/coloc_all_results.tsv" \
    --output-dir "results/meta_sensitivity/ancestry_tier1" \
    --i2-thresh 50

echo "Done: ${CELL_TYPE}"
