#!/bin/bash
# Generate Manhattan plots for all sn_NIMH_HBCC cell types
# Reads from results/sn_NIMH_HBCC/regenie_step2/*.regenie.raw_p
# Outputs to results/sn_NIMH_HBCC/manhattan_plots/

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

RESULTS_DIR="results/sn_NIMH_HBCC/regenie_step2"
OUT_DIR="results/sn_NIMH_HBCC/manhattan_plots"
SCRIPT="scripts/generate_cohort_manhattan.R"
COHORT="sn_NIMH_HBCC"

mkdir -p "$OUT_DIR"
mkdir -p logs/sn_NIMH_HBCC_manhattan

for raw_p in "${RESULTS_DIR}"/NIMH_HBCC_*_step2.regenie.raw_p; do
    # Extract cell type: NIMH_HBCC_Adaptive_step2.regenie.raw_p -> Adaptive
    fname=$(basename "$raw_p")
    cell_type=$(echo "$fname" | sed 's/^NIMH_HBCC_//; s/_step2\.regenie\.raw_p$//')

    echo "Submitting: ${cell_type}"
    sbatch \
        --job-name="mht_sn_hbcc_${cell_type}" \
        --partition=mediumtmp \
        --time=3:00:00 \
        --cpus-per-task=1 \
        --mem=32G \
        --output="logs/sn_NIMH_HBCC_manhattan/mht_${cell_type}_%j.out" \
        --error="logs/sn_NIMH_HBCC_manhattan/mht_${cell_type}_%j.err" \
        --wrap="Rscript ${SCRIPT} \
            --input ${raw_p} \
            --cell_type '${cell_type}' \
            --cohort '${COHORT}' \
            --output_prefix '${cell_type}' \
            --output_dir '${OUT_DIR}'"
done

echo "All Manhattan plot jobs submitted. Output: ${OUT_DIR}"
