#!/bin/bash
#SBATCH --job-name=gwas_nobio_allct
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/zhoux156/logs/validation/gwas_nobio_allct_%j.out
#SBATCH --error=/scratch/zhoux156/logs/validation/gwas_nobio_allct_%j.err

# PHENO_PREP + REGENIE for all 19 cell types, using the approved
# exclude_bio_from_tech_cov cell_proportions (skip Stage 1 RNA).
# 13 cohorts: SCC pgen. PENN/PITT: joint June 9 CMC.QC pgen from e2e.

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"

MAP="${NF_DIR}/validation/cohort_trillium_paths.tsv"
OUT_ROOT="${OUT_ROOT:-${RESULTS_ROOT}/validation/gwas_nobio_allct}"
PROP_ROOT="${PROP_ROOT:-${RESULTS_ROOT}/validation/deconv_nobio}"
META_ROOT="${META_ROOT:-${RESULTS_ROOT}/validation/e2e_june9}"
SCC_RES="${SCC_DIR}/results"
MAP_DIR="${NF_DIR}/validation/june9_fid_specimen_maps"
WORK_NAME="${WORK_NAME:-validation_gwas_nobio_allct}"
E2E_PGEN="${RESULTS_ROOT}/validation/e2e_june9"
mkdir -p "${OUT_ROOT}" "${LOG_ROOT}/validation/gwas_nobio_workers"

cfg_prefix() {
    case "$1" in
        NIMH_HBCC_*) echo "CMC_HBCC.QC.final" ;;
        *)           echo "${1}.QC.final" ;;
    esac
}
qc_study_name() {
    case "$1" in
        NIMH_HBCC_*) echo "CMC_HBCC" ;;
        *)           echo "$1" ;;
    esac
}
resolve_pgen_prefix() {
    local cohort="$1" scc_results="$2" out="$3" f
    if [ "${cohort}" = "CMC_PENN" ] || [ "${cohort}" = "CMC_PITT" ]; then
        for f in "${E2E_PGEN}/${cohort}"/*.QC.final.pgen; do
            [ -s "${f}" ] && { echo "${f%.pgen}"; return 0; }
        done
    fi
    for f in "${SCC_RES}/${scc_results}"/*.QC.final.pgen; do
        [ -s "${f}" ] && { echo "${f%.pgen}"; return 0; }
    done
    for f in "${out}"/*.QC.final.pgen; do
        [ -s "${f}" ] && { echo "${f%.pgen}"; return 0; }
    done
    return 1
}
config_tissue_filter() {
    python3 -c "
import re, sys
t=open(sys.argv[1]).read()
m=re.search(r\"tissue_filter\\s*=\\s*['\\\"]([^'\\\"]+)['\\\"]\", t)
print(m.group(1) if m else '')
" "$1"
}

run_one() {
    local cohort="$1" config="$2" count_matrix="$3" metadata="$4"
    local biospecimen="$5" clinical_metadata="$6" scc_results="$7"
    local work run out prefix cfgp qcs scc extra tissue
    echo "======== START ${cohort} $(date -Is) ========"
    work="${WORK_ROOT}/${WORK_NAME}/${cohort}"
    run="${SCRATCH_ROOT}/run_dirs/${WORK_NAME}/${cohort}"
    out="${OUT_ROOT}/${cohort}"
    mkdir -p "${work}" "${run}" "${out}/logs"
    cd "${run}"
    export NXF_HOME="${SCRATCH_ROOT}/.nextflow_gwas_nobio/${cohort}"
    export TMPDIR="${SCRATCH_ROOT}/tmp/gwas_nobio_${cohort}"
    export NXF_ANSI_LOG=false
    mkdir -p "${NXF_HOME}" "${TMPDIR}"

    if [ ! -s "${PROP_ROOT}/${cohort}/cell_proportions.csv" ]; then
        echo "  FAIL ${cohort}: missing ${PROP_ROOT}/${cohort}/cell_proportions.csv"
        touch "${FAIL_FILE}"; return 1
    fi
    if [ ! -s "${META_ROOT}/${cohort}/combined_metrics.csv" ]; then
        echo "  FAIL ${cohort}: missing combined_metrics.csv"
        touch "${FAIL_FILE}"; return 1
    fi

    qcs=$(qc_study_name "${cohort}")
    scc="${SCC_RES}/${scc_results}"
    extra=()
    prefix=""
    if ! prefix=$(resolve_pgen_prefix "${cohort}" "${scc_results}" "${out}"); then
        echo "  FAIL ${cohort}: no readable pgen"
        touch "${FAIL_FILE}"; return 1
    fi
    cfgp=$(cfg_prefix "${cohort}")
    if [ ! -s "${NF_DIR}/results/${cohort}/${cfgp}.pgen" ]; then
        echo "  FAIL ${cohort}: staged pgen missing ${NF_DIR}/results/${cohort}/${cfgp}.pgen"
        touch "${FAIL_FILE}"; return 1
    fi
    if [ -s "${scc}/pca.csv" ]; then
        [ "$(readlink -f "${scc}/pca.csv")" != "$(readlink -f "${out}/pca.csv" 2>/dev/null || true)" ] && \
            ln -sfn "${scc}/pca.csv" "${out}/pca.csv"
    elif [ -s "${E2E_PGEN}/${cohort}/pca.csv" ]; then
        ln -sfn "${E2E_PGEN}/${cohort}/pca.csv" "${out}/pca.csv"
    fi
    prune_src=""
    if [ -s "${scc}/${qcs}.QC.prune.in" ]; then
        prune_src="${scc}/${qcs}.QC.prune.in"
    elif [ -s "${E2E_PGEN}/${cohort}/${qcs}.QC.prune.in" ]; then
        prune_src="${E2E_PGEN}/${cohort}/${qcs}.QC.prune.in"
    elif [ -s "${E2E_PGEN}/${cohort}/${cohort}.QC.prune.in" ]; then
        prune_src="${E2E_PGEN}/${cohort}/${cohort}.QC.prune.in"
    fi
    if [ -n "${prune_src}" ]; then
        ln -sfn "${prune_src}" "${out}/${qcs}.QC.prune.in"
        ln -sfn "${prune_src}" "${out}/${cohort}.QC.prune.in" 2>/dev/null || true
    fi
    extra+=(--precomputed_pgen "${prefix}")
    extra+=(--precomputed_proportions "${PROP_ROOT}/${cohort}/cell_proportions.csv")
    extra+=(--precomputed_metadata "${META_ROOT}/${cohort}/combined_metrics.csv")
    [ -n "${biospecimen}" ] && [ "${biospecimen}" != "-" ] && extra+=(--biospecimen_file "${biospecimen}")
    [ -n "${clinical_metadata}" ] && [ "${clinical_metadata}" != "-" ] && extra+=(--clinical_metadata_file "${clinical_metadata}")
    extra+=(--work_dir "${work}/qc" --log_dir "${out}/logs")
    extra+=(--june9_replay true)
    extra+=(--regenie_step1_snp_filter mac1)
    tissue=$(config_tissue_filter "${NF_DIR}/${config}")
    [ -n "${tissue}" ] && extra+=(--pheno_tissue_filter "${tissue}")
    if [ -s "${MAP_DIR}/${cohort}.tsv" ]; then
        extra+=(--pheno_fid_specimen_map "${MAP_DIR}/${cohort}.tsv")
    fi

    if nextflow run "${NF_DIR}/combined_pipeline_v2.nf" \
        -resume \
        -c "${NF_DIR}/${config}" \
        -c "${NF_DIR}/nextflow.config.trillium" \
        -c "${NF_DIR}/nextflow.config.e2e_parallel_cohorts" \
        -work-dir "${work}" \
        --count_matrix_file "${count_matrix}" \
        --metadata_file "${metadata}" \
        --output_dir "${out}" \
        --skip_metal true \
        --marker_file "${REFS_ROOT}/markers/new_MTGnCgG_lfct2.5_Publication.csv" \
        --hgnc_mapping_file "${REFS_ROOT}/markers/hgnc_complete_set.txt" \
        "${extra[@]+"${extra[@]}"}" \
        > "${LOG_ROOT}/validation/gwas_nobio_workers/${cohort}.log" 2>&1
    then
        echo "  OK   ${cohort} $(date -Is)"
    else
        echo "  FAIL ${cohort} $(date -Is)"
        touch "${FAIL_FILE}"
        return 1
    fi
}

FAIL_FILE="${OUT_ROOT}/.fail"
rm -f "${FAIL_FILE}"
MAP_FILLED="${OUT_ROOT}/.cohort_map_filled.tsv"
sed -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' "${MAP}" > "${MAP_FILLED}"

echo "Launching all-cell-type GWAS $(date -Is) props=${PROP_ROOT}"
pids=()
while IFS=$'\t' read -r cohort config count_matrix metadata biospecimen clinical_metadata vcf_dir vcf_glob scc_results; do
    [ "${cohort}" = "cohort" ] && continue
    run_one "${cohort}" "${config}" "${count_matrix}" "${metadata}" \
        "${biospecimen}" "${clinical_metadata}" "${scc_results}" &
    pids+=($!)
    echo "  launched ${cohort} pid $!"
done < "${MAP_FILLED}"

wait_rc=0
for pid in "${pids[@]}"; do
    wait "${pid}" || wait_rc=1
done
echo "All cohort pipelines finished $(date -Is) wait_rc=${wait_rc}"
[ ! -f "${FAIL_FILE}" ]
