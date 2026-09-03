#!/bin/bash
#SBATCH --job-name=sn_vs_bulk
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/sn_vs_bulk_%j.out
#SBATCH --error=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/sn_vs_bulk_%j.err
#SBATCH --account=rrg-shreejoy

set -euo pipefail

cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow

SN_DIR="results/sn_ROSMAP/regenie_step2"

# Comparison 1: sn vs bulk with design matrix batch correction
echo "=== sn vs bulk (design matrix) ==="
mkdir -p results/sn_ROSMAP/sn_vs_bulk_design_matrix
Rscript scripts/compare_sn_vs_bulk_gwas.R \
    --sn_dir   "${SN_DIR}" \
    --bulk_dir results/ROSMAP/regenie_step2 \
    --output_dir results/sn_ROSMAP/sn_vs_bulk_design_matrix \
    --cohort ROSMAP

# Comparison 2: sn vs bulk with library-only batch correction
echo "=== sn vs bulk (library-only) ==="
mkdir -p results/sn_ROSMAP/sn_vs_bulk_libonly
Rscript scripts/compare_sn_vs_bulk_gwas.R \
    --sn_dir   "${SN_DIR}" \
    --bulk_dir results/tech_cov_policy/no_tech_except_rosmap_libprep/ROSMAP/regenie_step2 \
    --output_dir results/sn_ROSMAP/sn_vs_bulk_libonly \
    --cohort ROSMAP
