#!/bin/bash
#SBATCH --job-name=meta_15cohorts
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=12:00:00
#SBATCH --output=logs/meta_15cohorts/meta_15cohorts_%j.out
#SBATCH --error=logs/meta_15cohorts/meta_15cohorts_%j.err

# Trillium allocates whole nodes (SelectType=select/linear), so this job owns all
# 192 cores and 767 GB and Nextflow schedules tasks locally inside the allocation.
# Do not add --mem; node memory is allocated in full.

set -euo pipefail

NF_DIR=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow
SCC_DIR=/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow
WORK_DIR=/scratch/zhoux156/nf_work/meta_15cohorts

# Pin 25.10.2 to match the 25.10.0 used for the original run: Nextflow 26's config
# parser is stricter and rejects the cross-reference style these configs rely on
# (e.g. output_dir = "${base_dir}/...").
module load nextflow/25.10.2

# R (plot stages call `Rscript` bare) and python live in the imported `test` env.
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate test

mkdir -p "${NF_DIR}/run_dirs/meta_15cohorts/logs/meta_15cohorts" "${WORK_DIR}"
cd "${NF_DIR}/run_dirs/meta_15cohorts"

# Stage 0 (cohort audit) reads the per-cohort REGENIE output from the SCC tree.
echo "Checking GWAS results for all 15 cohorts..."
all_good=true
for cohort in ROSMAP ROSMAP_array Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx_v10 NABEC NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M GVEX AMP_AD_Rush AMP_AD_Mayo; do
    count=$(ls ${SCC_DIR}/results/${cohort}/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
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

# Liftover + harmonization are skipped: pre_harmonized_input_dir points at the
# already-harmonized, already-lifted sumstats from the CAMH run.
PRE_HARM=/project/rrg-shreejoy/zhoux156/harmonized_15cohorts
harm_count=$(ls ${PRE_HARM}/*/harmonized/*_harmonized.raw_p 2>/dev/null | wc -l)
echo ""
echo "Pre-harmonized inputs: ${harm_count}/285"
if [ "$harm_count" -ne 285 ]; then
    echo "ERROR: expected 285 pre-harmonized files (19 cell types x 15 cohorts). Aborting."
    exit 1
fi

echo ""
echo "All inputs present. Starting 15-cohort meta-analysis..."

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_15cohorts \
    -work-dir ${WORK_DIR}

echo ""
echo "15-cohort meta-analysis complete."
echo "Results written to: ${NF_DIR}/results/meta_analysis_15cohorts/"
