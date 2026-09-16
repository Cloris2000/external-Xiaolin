#!/bin/bash
#SBATCH --job-name=deconv_15cohort
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/zhoux156/logs/validation/deconv_15cohort_%j.out
#SBATCH --error=/scratch/zhoux156/logs/validation/deconv_15cohort_%j.err

# Stage 1+2 only (RNA → deconv → phenotypes_RINT) for all 15 bulk cohorts.
# Does not touch VCFs or REGENIE.
#
# Opt-in env (defaults keep the original deconv output path):
#   WORK_NAME=validation_deconv
#   OUT_ROOT / CMP_DIR
#   EXCLUDE_BIO=1  →  --exclude_bio_from_tech_cov true

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"

MAP="${NF_DIR}/validation/cohort_trillium_paths.tsv"
WORK_NAME="${WORK_NAME:-validation_deconv}"
OUT_ROOT="${OUT_ROOT:-${RESULTS_ROOT}/validation/deconv}"
CMP_DIR="${CMP_DIR:-${RESULTS_ROOT}/validation/deconv_compare}"
SCC_RES="${SCC_DIR}/results"
mkdir -p "${OUT_ROOT}" "${CMP_DIR}" "${LOG_ROOT}/validation/deconv_workers"

run_one() {
    local cohort="$1" config="$2" count_matrix="$3" metadata="$4"
    local biospecimen="$5" clinical_metadata="$6" scc_results="$7"
    local work run out
    echo "======== START ${cohort} $(date -Is) ========"
    work="${WORK_ROOT}/${WORK_NAME}/${cohort}"
    run="${SCRATCH_ROOT}/run_dirs/${WORK_NAME}/${cohort}"
    out="${OUT_ROOT}/${cohort}"
    mkdir -p "${work}" "${run}" "${out}"
    cd "${run}"
    export NXF_HOME="${SCRATCH_ROOT}/.nextflow_deconv/${cohort}"
    export TMPDIR="${SCRATCH_ROOT}/tmp/deconv_${cohort}"
    export NXF_ANSI_LOG=false
    mkdir -p "${NXF_HOME}" "${TMPDIR}"

    extra=()
    [ -n "${biospecimen}" ] && [ "${biospecimen}" != "-" ] && extra+=(--biospecimen_file "${biospecimen}")
    [ -n "${clinical_metadata}" ] && [ "${clinical_metadata}" != "-" ] && extra+=(--clinical_metadata_file "${clinical_metadata}")
    if [ "${EXCLUDE_BIO:-0}" = "1" ]; then
        extra+=(--exclude_bio_from_tech_cov true)
        echo "  exclude_bio_from_tech_cov=true"
    fi

    if nextflow run "${NF_DIR}/combined_pipeline_v2.nf" \
        -resume \
        -c "${NF_DIR}/${config}" \
        -c "${NF_DIR}/nextflow.config.trillium" \
        -c "${NF_DIR}/nextflow.config.e2e_parallel_cohorts" \
        -work-dir "${work}" \
        --stop_after pheno \
        --count_matrix_file "${count_matrix}" \
        --metadata_file "${metadata}" \
        --output_dir "${out}" \
        --marker_file "${REFS_ROOT}/markers/new_MTGnCgG_lfct2.5_Publication.csv" \
        --hgnc_mapping_file "${REFS_ROOT}/markers/hgnc_complete_set.txt" \
        "${extra[@]+"${extra[@]}"}" \
        > "${LOG_ROOT}/validation/deconv_workers/${cohort}.log" 2>&1
    then
        echo "  OK   ${cohort} $(date -Is)"
    else
        echo "  FAIL ${cohort} $(date -Is)"
        touch "${FAIL_FILE}"
        return 1
    fi

    new_prop=$(ls "${out}"/cell_proportions.csv 2>/dev/null | head -1 || true)
    new_pheno=$(ls "${out}"/phenotypes_RINT.txt 2>/dev/null | head -1 || true)
    "${RSCRIPT}" "${NF_DIR}/scripts/compare_deconv_to_scc.R" \
        --cohort "${cohort}" \
        --new_prop "${new_prop}" \
        --old_prop "${SCC_RES}/${scc_results}/cell_proportions.csv" \
        --new_pheno "${new_pheno}" \
        --old_pheno "${SCC_RES}/${scc_results}/phenotypes_RINT.txt" \
        --out "${CMP_DIR}/${cohort}.tsv" \
        >> "${LOG_ROOT}/validation/deconv_workers/${cohort}.log" 2>&1 || true
    return 0
}

FAIL_FILE="${OUT_ROOT}/.fail"
rm -f "${FAIL_FILE}"
MAP_FILLED="${OUT_ROOT}/.cohort_map_filled.tsv"
sed -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' "${MAP}" > "${MAP_FILLED}"

echo "Launching deconv in parallel $(date -Is) work=${WORK_NAME} exclude_bio=${EXCLUDE_BIO:-0}"
pids=()
while IFS=$'\t' read -r cohort config count_matrix metadata biospecimen clinical_metadata vcf_dir vcf_glob scc_results; do
    [ "${cohort}" = "cohort" ] && continue
    run_one "${cohort}" "${config}" "${count_matrix}" "${metadata}" \
        "${biospecimen}" "${clinical_metadata}" "${scc_results}" &
    pids+=($!)
    echo "  launched ${cohort} pid $!"
done < "${MAP_FILLED}"

echo "Waiting on ${#pids[@]} cohort pipelines..."
wait_rc=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        wait_rc=1
    fi
done
echo "All cohort pipelines finished $(date -Is) wait_rc=${wait_rc}"

echo ""
echo "Per-cohort comparison tables: ${CMP_DIR}"
if ls "${CMP_DIR}"/*.tsv >/dev/null 2>&1; then
    { head -1 "$(ls "${CMP_DIR}"/*.tsv | head -1)"; tail -n +2 "${CMP_DIR}"/*.tsv | grep -v '^cohort'; } \
      > "${CMP_DIR}/summary.tsv" || true
    cat "${CMP_DIR}/summary.tsv"
fi
[ ! -f "${FAIL_FILE}" ]
