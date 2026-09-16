#!/bin/bash
#SBATCH --job-name=meta_astro_june9
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/zhoux156/logs/validation/meta_astro_june9_%j.out
#SBATCH --error=/scratch/zhoux156/logs/validation/meta_astro_june9_%j.err

# Astrocyte-only meta from validation GWAS, June 9 lift settings.
# Compare the resulting .tbl to SCC Astrocyte *1.tbl.

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"

GWAS_ROOT="${RESULTS_ROOT}/validation/gwas_Astrocyte"
OUT_DIR="${RESULTS_ROOT}/validation/meta_astrocyte_june9"
WORK_DIR="${WORK_ROOT}/validation_meta_astrocyte_june9"
RUN_DIR="${SCRATCH_ROOT}/run_dirs/validation_meta_astrocyte_june9"
SCC_META="${SCC_DIR}/results/meta_analysis_15cohorts"
SUFFIX="${META_SUFFIX:-AMP_AD_Mayo_AMP_AD_Rush_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_v10_GVEX_MSBB_Mayo_NABEC_NIMH_HBCC_1M_NIMH_HBCC_Omni5M_NIMH_HBCC_h650_ROSMAP_ROSMAP_array}"

mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${OUT_DIR}" "${LOG_ROOT}/validation"
cd "${RUN_DIR}"

echo "Checking Astrocyte .raw_p in ${GWAS_ROOT}..."
all_good=true
for cohort in ROSMAP ROSMAP_array Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx_v10 NABEC NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M GVEX AMP_AD_Rush AMP_AD_Mayo; do
    n=$(ls "${GWAS_ROOT}/${cohort}"/regenie_step2/*Astrocyte*.regenie.raw_p 2>/dev/null | wc -l)
    if [ "${n}" -ge 1 ]; then
        echo "  OK  ${cohort}: ${n}"
    else
        echo "  FAIL ${cohort}: no Astrocyte .raw_p"
        all_good=false
    fi
done
if [ "${all_good}" = false ]; then
    echo "ERROR: missing Astrocyte GWAS outputs. Aborting."
    exit 1
fi

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/nextflow.config.standalone_meta_15cohorts_astrocyte_june9" \
    -work-dir "${WORK_DIR}" \
    -resume

TRI="${OUT_DIR}/Astrocyte_meta_analysis_${SUFFIX}.tbl"
REF="${SCC_META}/Astrocyte_meta_analysis_${SUFFIX}1.tbl"
CMP="${OUT_DIR}/compare_vs_june9.txt"
{
    echo "Trillium: ${TRI}"
    echo "June 9  : ${REF} (SCC *1.tbl)"
    if [ -s "${TRI}" ] && [ -s "${REF}" ]; then
        a=$(md5sum "${TRI}" | awk '{print $1}')
        b=$(md5sum "${REF}" | awk '{print $1}')
        echo "md5 tri : ${a}"
        echo "md5 ref : ${b}"
        if [ "${a}" = "${b}" ]; then
            echo "RESULT  : MATCH"
        else
            echo "RESULT  : DIFF"
            echo "bytes tri=$(stat -c%s "${TRI}") ref=$(stat -c%s "${REF}")"
        fi
    else
        echo "RESULT  : MISSING"
    fi
} | tee "${CMP}"

REPORT="${RESULTS_ROOT}/validation/MORNING_REPORT.txt"
{
    echo "=== validation morning report $(date -Is) ==="
    echo ""
    echo "--- deconv vs SCC ---"
    cat "${RESULTS_ROOT}/validation/deconv_compare/summary.tsv" 2>/dev/null || echo "(no deconv summary)"
    echo ""
    echo "--- Astrocyte .raw_p vs SCC ---"
    cat "${RESULTS_ROOT}/validation/gwas_Astrocyte_compare/summary.tsv" 2>/dev/null || echo "(no gwas summary)"
    echo ""
    echo "--- Astrocyte meta vs June 9 *1.tbl ---"
    cat "${CMP}"
} > "${REPORT}"
echo "Wrote ${REPORT}"
