#!/bin/bash
# Generate Manhattan plots for NIMH_HBCC_1M, NIMH_HBCC_Omni5M, NIMH_HBCC_h650
# Reads from results/<cohort>/regenie_step2/*.regenie.raw_p
# Outputs to results/<cohort>/manhattan_plots/

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

SCRIPT="scripts/generate_cohort_manhattan.R"

for COHORT in NIMH_HBCC_1M NIMH_HBCC_Omni5M NIMH_HBCC_h650; do
    RESULTS_DIR="results/${COHORT}/regenie_step2"
    OUT_DIR="results/${COHORT}/manhattan_plots"

    mkdir -p "$OUT_DIR"
    mkdir -p "logs/${COHORT}_manhattan"

    for raw_p in "${RESULTS_DIR}"/*.raw_p; do
        fname=$(basename "$raw_p")
        # e.g. NIMH_HBCC_1M_Astrocyte_step2.regenie.raw_p -> Astrocyte
        cell_type=$(echo "$fname" | sed "s/^${COHORT}_//; s/_step2\.regenie\.raw_p\$//")

        echo "Submitting: ${COHORT} / ${cell_type}"
        sbatch \
            --job-name="mht_${COHORT}_${cell_type}" \
            --partition=mediumtmp \
            --time=3:00:00 \
            --cpus-per-task=1 \
            --mem=32G \
            --output="logs/${COHORT}_manhattan/mht_${cell_type}_%j.out" \
            --error="logs/${COHORT}_manhattan/mht_${cell_type}_%j.err" \
            --wrap="Rscript ${SCRIPT} \
                --input ${raw_p} \
                --cell_type '${cell_type}' \
                --cohort '${COHORT}' \
                --output_prefix '${cell_type}' \
                --output_dir '${OUT_DIR}'"
    done
done

echo "All Manhattan plot jobs submitted."
