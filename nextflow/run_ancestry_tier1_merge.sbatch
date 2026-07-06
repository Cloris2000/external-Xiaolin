#!/bin/bash
#SBATCH --job-name=ancestry_merge
#SBATCH --partition=medium
#SBATCH --time=0:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/meta_sensitivity/ancestry_tier1_merge_%j.out
#SBATCH --error=logs/meta_sensitivity/ancestry_tier1_merge_%j.err

# Merge 19 per-cell-type outputs from run_ancestry_tier1_array.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

python3 "${NF_DIR}/scripts/ancestry_tier1_screen.py" \
    --project-dir "${NF_DIR}" \
    --merge-only \
    --output-dir "results/meta_sensitivity/ancestry_tier1"

echo ""
echo "Merged outputs:"
echo "  results/meta_sensitivity/ancestry_tier1/ancestry_comparison_screened_hits.tsv"
echo "  results/meta_sensitivity/ancestry_tier1/ancestry_masking_candidates.tsv"
echo "  results/meta_sensitivity/ancestry_tier1/ancestry_locus_detail.tsv"
