#!/bin/bash
#SBATCH --job-name=compare_mgp_sn
#SBATCH --partition=medium
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/compare_mgp_sn_%j.out
#SBATCH --error=logs/compare_mgp_sn_%j.err

set -euo pipefail

NF_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
mkdir -p "${NF_DIR}/logs"

echo "============================================"
echo "MGP vs snRNA-seq comparison (design= vs no design=)"
echo "Job ID : ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================"

source /nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh
conda activate test

Rscript "${NF_DIR}/scripts/compare_mgp_snrnaseq.R"

echo ""
echo "============================================"
echo "Done: $(date)"
echo "Outputs: ${NF_DIR}/results/ROSMAP/comparison_design_matrix/"
echo "============================================"
