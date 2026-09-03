#!/bin/bash
#SBATCH --job-name=compare_mgp_sn
#SBATCH --partition=compute
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/compare_mgp_sn_%j.out
#SBATCH --error=logs/compare_mgp_sn_%j.err
#SBATCH --account=rrg-shreejoy

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
mkdir -p "${NF_DIR}/logs"

echo "============================================"
echo "MGP vs snRNA-seq comparison (design= vs no design=)"
echo "Job ID : ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================"

conda activate test

Rscript "${NF_DIR}/scripts/compare_mgp_snrnaseq.R"

echo ""
echo "============================================"
echo "Done: $(date)"
echo "Outputs: ${NF_DIR}/results/ROSMAP/comparison_design_matrix/"
echo "============================================"
