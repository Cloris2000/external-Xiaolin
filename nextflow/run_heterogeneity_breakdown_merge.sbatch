#!/bin/bash
#SBATCH --job-name=het_merge
#SBATCH --partition=medium
#SBATCH --time=0:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/het_breakdown/het_breakdown_merge_%j.out
#SBATCH --error=logs/het_breakdown/het_breakdown_merge_%j.err

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/het_breakdown"
cd "${NF_DIR}"

OUT_DIR="${NF_DIR}/results/meta_analysis_15cohorts/heterogeneity_breakdown"

python3 "${NF_DIR}/scripts/heterogeneity_breakdown.py" \
    --merge-only \
    --output-dir "${OUT_DIR}"

echo "Merged TSVs in ${OUT_DIR}/"
