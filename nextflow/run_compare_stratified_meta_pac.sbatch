#!/bin/bash
#SBATCH --job-name=cmp_pac
#SBATCH --partition=medium
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/meta_sensitivity/compare_pac_%j.out
#SBATCH --error=logs/meta_sensitivity/compare_pac_%j.err

# Compare pooled meta vs AD / EUR / HBCC stratified metas (PAC analysis).
#
# Usage:
#   sbatch run_compare_stratified_meta_pac.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
OUT_DIR="${NF_DIR}/results/meta_sensitivity/stratified_comparison_pac"
mkdir -p "${NF_DIR}/logs/meta_sensitivity" "${OUT_DIR}/figures"
cd "${NF_DIR}"

python3 "${NF_DIR}/scripts/compare_stratified_meta.py" \
    --project-dir "${NF_DIR}" \
    --strata-tsv "${NF_DIR}/docs/stratified_comparison_pac_strata.tsv" \
    --output-dir "${OUT_DIR}" \
    --workers 4

Rscript "${NF_DIR}/scripts/plot_stratified_meta_comparison.R" \
    --hit-overlap "${OUT_DIR}/stratified_hit_overlap.tsv" \
    --coloc-tsv     "${NF_DIR}/results/coloc/coloc_all_results.tsv" \
    --output-dir    "${OUT_DIR}/figures"

echo ""
echo "Comparison: ${OUT_DIR}/"
echo "Figures:      ${OUT_DIR}/figures/"
