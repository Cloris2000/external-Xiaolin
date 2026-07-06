#!/bin/bash
# Submit Manhattan plot jobs for all cohorts whose .raw_p P column was wrong and has now been fixed:
#   - NIMH_HBCC_1M, NIMH_HBCC_Omni5M, NIMH_HBCC_h650 (19 cell types each)
#   - sn_NIMH_HBCC (27 cell types, regenerate existing plots)
#   - GVEX (19 cell types, generate for the first time)
#
# Stale plots are removed before submission so results are fresh.

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

SCRIPT="scripts/generate_cohort_manhattan.R"
RESULTS_BASE="results"

submit_cohort() {
    local COHORT="$1"        # e.g. NIMH_HBCC_1M, GVEX, sn_NIMH_HBCC
    local RAW_P_DIR="$2"     # directory containing *.raw_p files
    local FILE_PREFIX="$3"   # prefix to strip from filename to get cell_type (e.g. "NIMH_HBCC_1M_")

    local OUT_DIR="${RESULTS_BASE}/${COHORT}/manhattan_plots"
    local LOG_DIR="logs/${COHORT}_manhattan"

    mkdir -p "$OUT_DIR" "$LOG_DIR"

    # Remove stale plots (wrong P column or from cancelled runs)
    if ls "${OUT_DIR}"/*.png &>/dev/null; then
        echo "  Removing $(ls ${OUT_DIR}/*.png | wc -l) stale plots from ${OUT_DIR}"
        rm -f "${OUT_DIR}"/*.png
    fi

    local count=0
    for raw_p in "${RAW_P_DIR}"/*.raw_p; do
        fname=$(basename "$raw_p")
        cell_type=$(echo "$fname" | sed "s/^${FILE_PREFIX}//; s/_step2\.regenie\.raw_p\$//")

        sbatch \
            --job-name="mht_${COHORT}_${cell_type}" \
            --partition=mediumtmp \
            --time=3:00:00 \
            --cpus-per-task=1 \
            --mem=32G \
            --output="${LOG_DIR}/mht_${cell_type}_%j.out" \
            --error="${LOG_DIR}/mht_${cell_type}_%j.err" \
            --wrap="Rscript ${SCRIPT} \
                --input ${raw_p} \
                --cell_type '${cell_type}' \
                --cohort '${COHORT}' \
                --output_prefix '${cell_type}' \
                --output_dir '${OUT_DIR}'"
        count=$((count + 1))
    done
    echo "  Submitted ${count} jobs for ${COHORT}"
}

echo "=== Submitting Manhattan plot jobs ==="
echo ""

echo "[1/4] NIMH_HBCC_1M"
submit_cohort "NIMH_HBCC_1M" \
    "results/NIMH_HBCC_1M/regenie_step2" \
    "NIMH_HBCC_1M_"

echo ""
echo "[2/4] NIMH_HBCC_Omni5M"
submit_cohort "NIMH_HBCC_Omni5M" \
    "results/NIMH_HBCC_Omni5M/regenie_step2" \
    "NIMH_HBCC_Omni5M_"

echo ""
echo "[3/4] NIMH_HBCC_h650"
submit_cohort "NIMH_HBCC_h650" \
    "results/NIMH_HBCC_h650/regenie_step2" \
    "NIMH_HBCC_h650_"

echo ""
echo "[4/4] sn_NIMH_HBCC (regenerating existing plots)"
# sn_NIMH_HBCC files are named NIMH_HBCC_<CellType>_step2.regenie.raw_p
submit_cohort "sn_NIMH_HBCC" \
    "results/sn_NIMH_HBCC/regenie_step2" \
    "NIMH_HBCC_"

echo ""
echo "[5/5] GVEX (first-time generation)"
submit_cohort "GVEX" \
    "results/GVEX/regenie_step2" \
    "GVEX_"

echo ""
echo "All jobs submitted. Monitor with:"
echo "  squeue -u \$USER --format='%.10i %.30j %.8T %.10M' | grep mht_"
