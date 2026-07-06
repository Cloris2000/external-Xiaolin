#!/bin/bash
#SBATCH --job-name=meta_AD_EUR_merge
#SBATCH --partition=medium
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/meta_sensitivity/stratified_ad_eur_merge_%j.out
#SBATCH --error=logs/meta_sensitivity/stratified_ad_eur_merge_%j.err

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
OUTPUT_DIR="${NF_DIR}/results/meta_sensitivity/stratified_ad_eur_15cohorts"
HET_DIR="${OUTPUT_DIR}/heterogeneity"
AUDIT_DIR="${OUTPUT_DIR}/audit"

mkdir -p "${NF_DIR}/logs/meta_sensitivity" "${HET_DIR}" "${AUDIT_DIR}"
cd "${NF_DIR}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

TBLS=()
for ct in "${CELL_TYPES[@]}"; do
    tbl=$(ls -1t "${OUTPUT_DIR}/${ct}"_meta_analysis_*.tbl 2>/dev/null | head -1 || true)
    if [[ -z "${tbl}" ]]; then
        echo "ERROR: missing METAL .tbl for cell type ${ct}"
        exit 1
    fi
    TBLS+=("${tbl}")
done

echo "Found ${#TBLS[@]} METAL .tbl files — running heterogeneity summary..."
python3 "${NF_DIR}/scripts/summarize_meta_heterogeneity.py" \
    --meta-tbls "${TBLS[@]}" \
    --output-summary-tsv "${HET_DIR}/meta_heterogeneity_summary.tsv" \
    --output-top-hits-tsv "${HET_DIR}/meta_heterogeneity_top_hits.tsv" \
    --gw-p-thresh 5e-8 \
    --suggestive-p-thresh 1e-5 \
    --i2-thresh 50 \
    --het-p-thresh 0.05

echo "Running cohort QC audit (once per stratum)..."
python3 "${NF_DIR}/scripts/audit_meta_cohorts.py" \
    --metadata "${NF_DIR}/docs/meta_cohort_metadata.tsv" \
    --base-dir "${NF_DIR}" \
    --cohorts-json '["ROSMAP","ROSMAP_array","Mayo","MSBB"]' \
    --output-tsv "${AUDIT_DIR}/meta_cohort_qc_matrix.tsv" \
    --output-md "${AUDIT_DIR}/meta_cohort_qc_summary.md"

echo "AD_EUR merge complete: ${OUTPUT_DIR}/"
