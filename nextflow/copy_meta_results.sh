#!/bin/bash
# Copy METAL meta-analysis .tbl files from work dirs to results/meta_analysis/
# for NABEC and NIMH HBCC cohorts
set -e
NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

for cohort in NABEC NIMH_HBCC_1M NIMH_HBCC_Omni5M NIMH_HBCC_h650; do
  meta_dir=${NF_DIR}/results/${cohort}/meta_analysis
  mkdir -p ${meta_dir}
  nf_log=${NF_DIR}/run_dirs/${cohort}/.nextflow.log

  # Extract work dir hashes for METAL_META_ANALYSIS steps from nextflow log
  grep "Submitted process.*METAL_META_ANALYSIS" ${nf_log} 2>/dev/null | \
    grep -oP '\[([0-9a-f]{2}/[0-9a-f]+)\]' | tr -d '[]' | while read hash; do
    prefix=${hash%%/*}
    suffix=${hash##*/}
    wdir="${NF_DIR}/work/${cohort}/${prefix}"
    full=$(ls -d ${wdir}/${suffix}* 2>/dev/null | head -1)
    if [ -n "$full" ]; then
      cp ${full}/*.tbl ${meta_dir}/ 2>/dev/null || true
      cp ${full}/*.tbl.info ${meta_dir}/ 2>/dev/null || true
    fi
  done

  n=$(ls ${meta_dir}/*.tbl 2>/dev/null | wc -l)
  echo "${cohort}: ${n} .tbl files copied to ${meta_dir}"
done

# Also fix NABEC stale covariates.txt
echo ""
echo "Fixing NABEC stale covariates.txt..."
nabec_cov_work=$(grep "Submitted process.*REGENIE_STEP1" ${NF_DIR}/run_dirs/NABEC/.nextflow.log 2>/dev/null | \
  grep -oP 'workDir: \K\S+' | tail -1)
if [ -z "$nabec_cov_work" ]; then
  # Try finding it via COV_PREP
  nabec_cov_work=$(grep "Submitted process.*COV_PREP\|PHENO_PREP" ${NF_DIR}/run_dirs/NABEC/.nextflow.log 2>/dev/null | \
    grep -oP '\[([0-9a-f]{2}/[0-9a-f]+)\]' | tr -d '[]' | head -1)
fi
# Find covariates.txt in work dirs newer than the stale one
find ${NF_DIR}/work/NABEC -maxdepth 4 -name "covariates.txt" \
  -newer ${NF_DIR}/results/NABEC/top_tech_covariates.txt 2>/dev/null | \
  xargs -I{} sh -c 'lines=$(wc -l < {}); echo "  {} : $lines lines"' | head -5
