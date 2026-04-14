#!/bin/bash
#SBATCH --job-name=meta_13cohorts
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_13cohorts_%j.out
#SBATCH --error=logs/meta_13cohorts_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/meta_13cohorts

# Verify all 13 cohorts have complete raw_p files (19 cell types each)
echo "Checking GWAS results for all 13 cohorts..."
all_good=true
for cohort in ROSMAP Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx NABEC NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M GVEX ROSMAP_array; do
    count=$(ls ${NF_DIR}/results/${cohort}/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
    if [ "$count" -eq 19 ]; then
        echo "  OK  ${cohort}: ${count}/19"
    else
        echo "  FAIL ${cohort}: ${count}/19 (INCOMPLETE)"
        all_good=false
    fi
done

if [ "$all_good" = false ]; then
    echo ""
    echo "ERROR: Some cohorts have incomplete GWAS results. Aborting."
    exit 1
fi

echo ""
echo "All 13 cohorts complete. Starting 13-cohort meta-analysis..."

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_13cohorts \
    -work-dir ${NF_DIR}/work/meta_13cohorts

echo ""
echo "13-cohort meta-analysis complete."
echo "Results written to: ${NF_DIR}/results/meta_analysis_13cohorts_design_matrix/"
