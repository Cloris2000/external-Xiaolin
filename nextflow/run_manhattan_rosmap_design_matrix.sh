#!/bin/bash
# Generate Manhattan plots for bulk ROSMAP (design matrix batch correction)
# Reads from results/ROSMAP/regenie_step2/*.regenie.raw_p
# Outputs to results/ROSMAP/manhattan_plots/

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

RESULTS_DIR="results/ROSMAP/regenie_step2"
OUT_DIR="results/ROSMAP/manhattan_plots"
SCRIPT="scripts/generate_cohort_manhattan.R"
COHORT="ROSMAP"

mkdir -p "$OUT_DIR"
mkdir -p logs/ROSMAP_manhattan

for raw_p in "${RESULTS_DIR}"/ROSMAP_*_step2.regenie.raw_p; do
    fname=$(basename "$raw_p")
    cell_type=$(echo "$fname" | sed 's/^ROSMAP_//; s/_step2\.regenie\.raw_p$//')

    echo "Submitting: ${cell_type}"
    sbatch \
        --job-name="mht_ros_dm_${cell_type}" \
        --partition=mediumtmp \
        --time=1:00:00 \
        --cpus-per-task=1 \
        --mem=16G \
        --output="logs/ROSMAP_manhattan/mht_${cell_type}_%j.out" \
        --error="logs/ROSMAP_manhattan/mht_${cell_type}_%j.err" \
        --wrap="Rscript ${SCRIPT} \
            --input ${raw_p} \
            --cell_type '${cell_type}' \
            --cohort '${COHORT}' \
            --output_prefix '${cell_type}' \
            --output_dir '${OUT_DIR}'"
done

echo "All Manhattan plot jobs submitted. Output: ${OUT_DIR}"
