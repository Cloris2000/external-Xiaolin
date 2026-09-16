#!/bin/bash
#SBATCH --job-name=meta_nobio_allct
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/zhoux156/logs/validation/meta_nobio_allct_%j.out
#SBATCH --error=/scratch/zhoux156/logs/validation/meta_nobio_allct_%j.err

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"

GWAS_ROOT="${RESULTS_ROOT}/validation/gwas_nobio_allct"
OUT_DIR="${RESULTS_ROOT}/validation/gwas_nobio_allct_meta"
WORK_DIR="${WORK_ROOT}/validation_gwas_nobio_allct_meta"
RUN_DIR="${SCRATCH_ROOT}/run_dirs/validation_gwas_nobio_allct_meta"
mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${OUT_DIR}"
cd "${RUN_DIR}"

all_good=true
for cohort in ROSMAP ROSMAP_array Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx_v10 NABEC NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M GVEX AMP_AD_Rush AMP_AD_Mayo; do
    n=$(ls "${GWAS_ROOT}/${cohort}"/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
    if [ "${n}" -eq 19 ]; then
        echo "  OK  ${cohort}: ${n}/19"
    else
        echo "  FAIL ${cohort}: ${n}/19"
        all_good=false
    fi
done
[ "${all_good}" = true ] || { echo "ERROR: incomplete GWAS"; exit 1; }

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/nextflow.config.standalone_meta_15cohorts_nobio_allct" \
    -work-dir "${WORK_DIR}"

echo "Meta written to ${OUT_DIR}"
