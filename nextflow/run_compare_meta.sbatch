#!/bin/bash
#SBATCH --job-name=compare_meta
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/scratch/zhoux156/logs/compare_meta_%j.out
#SBATCH --error=/scratch/zhoux156/logs/compare_meta_%j.err
#SBATCH --account=rrg-shreejoy

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
mkdir -p "${NF_DIR}/logs"

echo "============================================"
echo "Meta-analysis comparison: design vs no-tech"
echo "Job ID : ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================"

conda activate test

Rscript "${NF_DIR}/scripts/compare_meta_analyses.R"

echo ""
echo "Done: $(date)"
echo "Results: ${NF_DIR}/results/meta_analysis_comparison/"
