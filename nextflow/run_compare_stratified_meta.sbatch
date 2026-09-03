#!/bin/bash
#SBATCH --job-name=compare_strat_meta
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/meta_sensitivity/compare_stratified_meta_%j.out
#SBATCH --error=logs/meta_sensitivity/compare_stratified_meta_%j.err
#SBATCH --account=rrg-shreejoy

# Compare pooled 15-cohort meta vs stratified meta runs.
# Run after stratified meta jobs complete:
#   sbatch run_compare_stratified_meta.sbatch

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

python3 "${NF_DIR}/scripts/compare_stratified_meta.py" \
    --project-dir "${NF_DIR}" \
    --output-dir "${NF_DIR}/results/meta_sensitivity/stratified_comparison"

echo ""
echo "Comparison complete:"
echo "  ${NF_DIR}/results/meta_sensitivity/stratified_comparison/"
