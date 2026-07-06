#!/bin/bash
#SBATCH --job-name=plot_anc_t1
#SBATCH --partition=medium
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/meta_sensitivity/plot_ancestry_tier1_%j.out
#SBATCH --error=logs/meta_sensitivity/plot_ancestry_tier1_%j.err

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

Rscript "${NF_DIR}/scripts/plot_ancestry_tier1.R" \
    --screened-hits "${NF_DIR}/results/meta_sensitivity/ancestry_tier1/ancestry_comparison_screened_hits.tsv" \
    --locus-detail   "${NF_DIR}/results/meta_sensitivity/ancestry_tier1/ancestry_locus_detail.tsv" \
    --output-dir     "${NF_DIR}/results/meta_sensitivity/ancestry_tier1/figures" \
    --forest-dir     "${NF_DIR}/results/meta_sensitivity/ancestry_tier1/forest_plots" \
    --max-forest-plots 30 \
    --max-labels 15

echo "Figures:  results/meta_sensitivity/ancestry_tier1/figures/"
echo "Forests:  results/meta_sensitivity/ancestry_tier1/forest_plots/"
