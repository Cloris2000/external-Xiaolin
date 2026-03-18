#!/bin/bash
# Generate Manhattan plots for SST in each cohort from individual script (WGS) outputs
# Same cohorts as in Joint_AMP_AD_meta_with_CMC_NABEC_GTEx/SST_metal_script.txt

set -e
cd "$(dirname "$0")"
OUT_DIR="manhattan_plots/SST_ind_scripts"
mkdir -p "$OUT_DIR"

echo "Generating Manhattan plots for SST (individual scripts / WGS outputs)..."
echo "Output: $OUT_DIR"
echo ""

declare -A inputs
inputs[ROSMAP]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2/ROSMAP_WGS_step2_update_SST.regenie.raw_p"
inputs[Mayo]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/Mayo_joint_wgs_step2/Mayo_WGS_step2_update_SST.regenie.raw_p"
inputs[MSBB]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/MSBB_joint_wgs_step2/MSBB_WGS_step2_update_SST.regenie.raw_p"
inputs[NABEC]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/NABEC_wgs_step2/NABEC_WGS_step2_update_hg19_SST.regenie.raw_p"
inputs[CMC_MSSM]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/CMC_MSSM_SNP_array_step2/CMC_MSSM_SNP_array_reimputed_step2_SST.regenie.raw_p"
inputs[CMC_PENN]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/CMC_PENN_SNP_array_step2/CMC_PENN_SNP_array_reimputed_step2_SST.regenie.raw_p"
inputs[CMC_PITT]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/CMC_PITT_SNP_array_step2/CMC_PITT_SNP_array_reimputed_step2_SST.regenie.raw_p"
inputs[GTEx]="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/GTEx_WGS/step2/GTEx_WGS_step2_SST.regenie.raw_p"

for cohort in ROSMAP Mayo MSBB NABEC CMC_MSSM CMC_PENN CMC_PITT GTEx; do
  inp="${inputs[$cohort]}"
  if [ ! -f "$inp" ]; then
    echo "WARNING: Input not found: $inp (skipping)"
    continue
  fi
  echo "=== $cohort ==="
  Rscript scripts/generate_cohort_manhattan.R \
    --input "$inp" \
    --output_prefix "${cohort}_SST_ind" \
    --cell_type "SST" \
    --cohort "${cohort} (ind)" \
    --output_dir "$OUT_DIR"
  echo ""
done

echo "Done! Plots saved in $OUT_DIR"
ls -la "$OUT_DIR"/*_manhattan_annotated.png 2>/dev/null || true
