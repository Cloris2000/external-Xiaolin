#!/bin/bash
#SBATCH --job-name=cmp_diag_eur
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/scratch/zhoux156/logs/meta_sensitivity/compare_diagnosis_eur_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_sensitivity/compare_diagnosis_eur_%j.err
#SBATCH --account=rrg-shreejoy

# Compare pooled vs AD_EUR / psychiatric_EUR / EUR_homogeneous (parallel by cell type).
# Plots run in a separate quick job after tables are written.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
OUT_DIR="${NF_DIR}/results/meta_sensitivity/stratified_comparison_diagnosis_eur"
mkdir -p "${NF_DIR}/logs/meta_sensitivity" "${OUT_DIR}/figures"
cd "${NF_DIR}"

python3 "${NF_DIR}/scripts/compare_stratified_meta.py" \
    --project-dir "${NF_DIR}" \
    --strata-tsv "${NF_DIR}/docs/stratified_comparison_diagnosis_eur_strata.tsv" \
    --output-dir "${OUT_DIR}" \
    --workers 4

Rscript "${NF_DIR}/scripts/plot_diagnosis_eur_comparison.R" \
    --hit-overlap "${OUT_DIR}/stratified_hit_overlap.tsv" \
    --coloc-tsv     "${NF_DIR}/results/coloc/coloc_all_results.tsv" \
    --output-dir    "${OUT_DIR}/figures"

echo "Outputs: ${OUT_DIR}/"
