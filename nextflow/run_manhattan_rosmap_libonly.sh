#!/bin/bash
# Generate Manhattan plots for bulk ROSMAP (library-only batch correction)
# Reads from results/tech_cov_policy/no_tech_except_rosmap_libprep/ROSMAP/regenie_step2/
# Outputs to results/tech_cov_policy/no_tech_except_rosmap_libprep/ROSMAP/manhattan_plots/

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

RESULTS_DIR="results/tech_cov_policy/no_tech_except_rosmap_libprep/ROSMAP/regenie_step2"
OUT_DIR="results/tech_cov_policy/no_tech_except_rosmap_libprep/ROSMAP/manhattan_plots"
SCRIPT="scripts/generate_cohort_manhattan.R"
COHORT="ROSMAP_libonly"

mkdir -p "$OUT_DIR"
mkdir -p logs/ROSMAP_libonly_manhattan

for raw_p in "${RESULTS_DIR}"/ROSMAP_*_step2.regenie.raw_p; do
    fname=$(basename "$raw_p")
    cell_type=$(echo "$fname" | sed 's/^ROSMAP_//; s/_step2\.regenie\.raw_p$//')

    echo "Submitting: ${cell_type}"
    sbatch \
        --job-name="mht_ros_lib_${cell_type}" \
        --partition=mediumtmp \
        --time=1:00:00 \
        --cpus-per-task=1 \
        --mem=16G \
        --output="logs/ROSMAP_libonly_manhattan/mht_${cell_type}_%j.out" \
        --error="logs/ROSMAP_libonly_manhattan/mht_${cell_type}_%j.err" \
        --wrap="Rscript ${SCRIPT} \
            --input ${raw_p} \
            --cell_type '${cell_type}' \
            --cohort '${COHORT}' \
            --output_prefix '${cell_type}' \
            --output_dir '${OUT_DIR}'"
done

echo "All Manhattan plot jobs submitted. Output: ${OUT_DIR}"
