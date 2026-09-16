#!/bin/bash
#SBATCH --job-name=meta_15cohorts_replicate
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_replicate_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_replicate_%j.err

# Track A: reproduce the May 21 SCC meta (no liftover) from raw REGENIE output.
# Trillium allocates whole nodes; do not add --mem.
#
# Geometry is 1 task x 192 cpus so Nextflow's local executor sees the whole node.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
WORK_DIR="${WORK_ROOT}/meta_15cohorts_replicate"
RUN_DIR="${SCRATCH_ROOT}/run_dirs/meta_15cohorts_replicate"
OUT_DIR="${RESULTS_ROOT}/meta_replicate_15cohorts"

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

echo ""
echo "No cohort will be lifted (cohort_genome_builds is empty)."
echo "Starting Track A replication meta-analysis..."

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_15cohorts_replicate \
    -work-dir ${WORK_DIR} \
    -resume

echo ""
echo "Track A replication complete."
echo "Results written to: ${OUT_DIR}/"
echo "Verify against the SCC 1.63 GB tables with:"
echo "  TRI_DIR=${OUT_DIR} REF_KIND=published ${NF_DIR}/scripts/verify_replication.sh"
