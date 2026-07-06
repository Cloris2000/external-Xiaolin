#!/bin/bash
#SBATCH --job-name=het_plots
#SBATCH --partition=medium
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/het_breakdown/het_breakdown_plots_%j.out
#SBATCH --error=logs/het_breakdown/het_breakdown_plots_%j.err

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
OUT_DIR="${NF_DIR}/results/meta_analysis_15cohorts/heterogeneity_breakdown"
mkdir -p "${NF_DIR}/logs/het_breakdown" "${OUT_DIR}/figures" "${OUT_DIR}/forest_plots"
cd "${NF_DIR}"

Rscript "${NF_DIR}/scripts/plot_heterogeneity_breakdown.R" \
    --per-cohort-effects "${OUT_DIR}/per_cohort_effects.tsv" \
    --group-i2-summary   "${OUT_DIR}/group_i2_summary.tsv" \
    --overall-group-i2   "${OUT_DIR}/overall_group_i2.tsv" \
    --output-dir         "${OUT_DIR}/figures" \
    --forest-dir         "${OUT_DIR}/forest_plots" \
    --max-forest-plots   100

echo "Figures:  ${OUT_DIR}/figures/"
echo "Forests:  ${OUT_DIR}/forest_plots/"
