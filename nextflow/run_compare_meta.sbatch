#!/bin/bash
#SBATCH --job-name=compare_meta
#SBATCH --partition=medium
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/compare_meta_%j.out
#SBATCH --error=logs/compare_meta_%j.err

set -euo pipefail
NF_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
mkdir -p "${NF_DIR}/logs"

echo "============================================"
echo "Meta-analysis comparison: design vs no-tech"
echo "Job ID : ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================"

source /nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh
conda activate test

Rscript "${NF_DIR}/scripts/compare_meta_analyses.R"

echo ""
echo "Done: $(date)"
echo "Results: ${NF_DIR}/results/meta_analysis_comparison/"
