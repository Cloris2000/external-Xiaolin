#!/bin/bash
# Generate Manhattan plots for SST, VIP, and L5.ET meta-analysis results
# 8 cohorts: CMC_MSSM, CMC_PENN, CMC_PITT, GTEx, MSBB, Mayo, NABEC, ROSMAP

set -e
cd "$(dirname "$0")"
mkdir -p manhattan_plots/meta_8cohorts
META_DIR="results/meta_analysis"
OUT_DIR="manhattan_plots/meta_8cohorts"
SUFFIX="CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP"

echo "Generating Manhattan plots for meta-analysis (8 cohorts)..."
echo "Output: $OUT_DIR"
echo ""

for ct in SST VIP L5.ET; do
  inp="${META_DIR}/${ct}_meta_analysis_${SUFFIX}.tbl"
  if [ ! -f "$inp" ]; then
    echo "ERROR: Input not found: $inp"
    exit 1
  fi
  echo "=== $ct ==="
  Rscript scripts/generate_meta_manhattan.R \
    --input "$inp" \
    --output_prefix "${ct}_meta_8cohorts" \
    --cell_type "${ct} (8 cohorts meta)" \
    --output_dir "$OUT_DIR"
  echo ""
done

echo "Done! Plots saved in $OUT_DIR"
ls -la "$OUT_DIR"/*.png
