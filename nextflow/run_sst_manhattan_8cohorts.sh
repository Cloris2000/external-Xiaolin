#!/bin/bash
# Generate Manhattan plots for SST in each of the 8 cohorts separately
# Cohorts: CMC_MSSM, CMC_PENN, CMC_PITT, GTEx, MSBB, Mayo, NABEC, ROSMAP

set -e
cd "$(dirname "$0")"
RESULTS_DIR="results"
OUT_DIR="manhattan_plots/SST_per_cohort"
mkdir -p "$OUT_DIR"

echo "Generating Manhattan plots for SST (8 cohorts, separately)..."
echo "Output: $OUT_DIR"
echo ""

for cohort in CMC_MSSM CMC_PENN CMC_PITT GTEx MSBB Mayo NABEC ROSMAP; do
  inp="${RESULTS_DIR}/${cohort}/regenie_step2/${cohort}_SST_step2.regenie.raw_p"
  if [ ! -f "$inp" ]; then
    echo "WARNING: Input not found: $inp (skipping)"
    continue
  fi
  echo "=== $cohort ==="
  Rscript scripts/generate_cohort_manhattan.R \
    --input "$inp" \
    --output_prefix "${cohort}_SST" \
    --cell_type "SST" \
    --cohort "$cohort" \
    --output_dir "$OUT_DIR"
  echo ""
done

echo "Done! Plots saved in $OUT_DIR"
ls -la "$OUT_DIR"/*.png 2>/dev/null || true
