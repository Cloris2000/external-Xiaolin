#!/bin/bash
#SBATCH --job-name=sn_vs_bulk_hbcc
#SBATCH --partition=mediumtmp
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/sn_vs_bulk_hbcc_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/sn_vs_bulk_hbcc_%j.err

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

SN_DIR="results/sn_NIMH_HBCC/regenie_step2"

# ── Comparison: sn_HBCC vs bulk_HBCC (NIMH_HBCC_1M as representative) ────────
echo "=== sn_HBCC vs bulk_HBCC (1M platform) ==="
mkdir -p results/sn_NIMH_HBCC/sn_vs_bulk_comparison
Rscript scripts/compare_sn_vs_bulk_gwas.R \
    --sn_dir    "${SN_DIR}" \
    --bulk_dir  results/NIMH_HBCC_1M/regenie_step2 \
    --output_dir results/sn_NIMH_HBCC/sn_vs_bulk_comparison \
    --cohort     NIMH_HBCC \
    --sn_cohort  NIMH_HBCC \
    --bulk_cohort NIMH_HBCC_1M
