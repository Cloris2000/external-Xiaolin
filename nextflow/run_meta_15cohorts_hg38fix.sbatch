#!/bin/bash
#SBATCH --job-name=meta_15cohorts_hg38fix
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_hg38fix_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_hg38fix_%j.err

# Track B: 15-cohort meta-analysis with ROSMAP_array correctly declared hg38.
#
# Unlike run_meta_15cohorts.sbatch this does NOT reuse the pre-harmonized CAMH
# sumstats (they bake in the wrong build for ROSMAP_array), so it runs the full
# LIFTOVER_SUMSTATS -> HARMONIZE_META_SUMSTATS -> METAL path from the raw REGENIE
# output. That is why the walltime is 24 h rather than 12 h.
#
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
WORK_DIR="${WORK_ROOT}/meta_15cohorts_hg38fix"
RUN_DIR="${SCRATCH_ROOT}/run_dirs/meta_15cohorts_hg38fix"
OUT_DIR="${RESULTS_ROOT}/meta_analysis_15cohorts_hg38fix"

module load "${NEXTFLOW_MODULE}"

mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${OUT_DIR}"
cd "${RUN_DIR}"

echo "Checking raw REGENIE output for all 15 cohorts..."
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

CHAIN=/project/rrg-shreejoy/pipeline_refs/reference_data/hg38ToHg19.over.chain.gz
if [ ! -s "${CHAIN}" ]; then
    echo "ERROR: liftover chain file missing: ${CHAIN}" >&2
    exit 1
fi

echo ""
echo "Liftover chain present. 4 cohorts will be lifted hg38 -> hg19:"
echo "  GTEx_v10, AMP_AD_Rush, AMP_AD_Mayo, ROSMAP_array"
echo "Starting corrected 15-cohort meta-analysis..."

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_15cohorts_hg38fix \
    -work-dir ${WORK_DIR} \
    -resume

echo ""
echo "Corrected 15-cohort meta-analysis complete."
echo "Results written to: ${OUT_DIR}/"
echo "Sync them to ${PROJECT_RESULTS}/ from a login node with:"
echo "  ${NF_DIR}/scripts/sync_results_to_project.sh meta_analysis_15cohorts_hg38fix"
