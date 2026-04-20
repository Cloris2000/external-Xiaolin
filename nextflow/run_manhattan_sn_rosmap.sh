#!/bin/bash
# Generate Manhattan plots for all sn_ROSMAP cell types
# Reads from results/sn_ROSMAP/regenie_step2/*.regenie.raw_p
# Outputs to results/sn_ROSMAP/manhattan_plots/

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

RESULTS_DIR="results/sn_ROSMAP/regenie_step2"
OUT_DIR="results/sn_ROSMAP/manhattan_plots"
SCRIPT="scripts/generate_cohort_manhattan.R"
COHORT="sn_ROSMAP"

mkdir -p "$OUT_DIR"
mkdir -p logs/sn_ROSMAP_manhattan

for raw_p in "${RESULTS_DIR}"/ROSMAP_*_step2.regenie.raw_p; do
    # Extract cell type from filename: ROSMAP_Astrocyte_step2.regenie.raw_p -> Astrocyte
    fname=$(basename "$raw_p")
    cell_type=$(echo "$fname" | sed 's/^ROSMAP_//; s/_step2\.regenie\.raw_p$//')

    echo "Submitting: ${cell_type}"
    sbatch \
        --job-name="mht_sn_ros_${cell_type}" \
        --partition=mediumtmp \
        --time=1:00:00 \
        --cpus-per-task=1 \
        --mem=16G \
        --output="logs/sn_ROSMAP_manhattan/mht_${cell_type}_%j.out" \
        --error="logs/sn_ROSMAP_manhattan/mht_${cell_type}_%j.err" \
        --wrap="Rscript ${SCRIPT} \
            --input ${raw_p} \
            --cell_type '${cell_type}' \
            --cohort '${COHORT}' \
            --output_prefix '${cell_type}' \
            --output_dir '${OUT_DIR}'"
done

echo "All Manhattan plot jobs submitted. Output: ${OUT_DIR}"
