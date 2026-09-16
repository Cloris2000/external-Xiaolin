#!/bin/bash
#SBATCH --job-name=gwas_astro_15
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/zhoux156/logs/validation/gwas_astro_15_%j.out
#SBATCH --error=/scratch/zhoux156/logs/validation/gwas_astro_15_%j.err

# Astrocyte × 15: REGENIE only, using SCC (June 9) pgen + phenotypes.
# Skips VCF QC. CMC_PENN / CMC_PITT SCC pgen are dangling CAMH links, so
# those two reuse last night's Trillium-QC pgen (will not match June 9).
# REGENIE step 2 reads ${projectDir}/results/<cohort>/*.QC.final — those
# symlinks must be created on a login node before submit (see
# validation/stage_scc_pgen_for_regenie.sh).

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"

MAP="${NF_DIR}/validation/cohort_trillium_paths.tsv"
CT="${CELL_TYPE:-Astrocyte}"
OUT_ROOT="${RESULTS_ROOT}/validation/gwas_${CT}"
CMP="${RESULTS_ROOT}/validation/gwas_${CT}_compare"
SCC_RES="${SCC_DIR}/results"
# Fresh work dir so we do not -resume last night's broken QC / missing-pgen tasks.
WORK_NAME="${WORK_NAME:-validation_gwas_${CT}_pgen}"
mkdir -p "${OUT_ROOT}" "${CMP}" "${LOG_ROOT}/validation/gwas_workers"

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

# Prefer a readable SCC June 9 pgen; else last night's scratch QC pgen.
resolve_pgen_prefix() {
    local cohort="$1" scc_results="$2" out="$3"
    local f
    for f in "${SCC_RES}/${scc_results}"/*.QC.final.pgen; do
        if [ -s "${f}" ]; then
            echo "${f%.pgen}"
            return 0
        fi
    done
    for f in "${out}"/*.QC.final.pgen; do
        if [ -s "${f}" ]; then
            echo "${f%.pgen}"
            return 0
        fi
    done
    return 1
}

fail=0
FAIL_FILE="${OUT_ROOT}/.fail"
rm -f "${FAIL_FILE}"
MAP_FILLED="${OUT_ROOT}/.cohort_map_filled.tsv"
sed -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' "${MAP}" > "${MAP_FILLED}"
while IFS=$'\t' read -r cohort config count_matrix metadata biospecimen clinical_metadata vcf_dir vcf_glob scc_results; do
    [ "${cohort}" = "cohort" ] && continue
    echo "======== ${cohort} ${CT} (precomputed pgen) ========"
    work="${WORK_ROOT}/${WORK_NAME}/${cohort}"
    run="${SCRATCH_ROOT}/run_dirs/${WORK_NAME}/${cohort}"
    out="${OUT_ROOT}/${cohort}"
    mkdir -p "${work}" "${run}" "${out}"
    cd "${run}"

    if ls "${out}"/regenie_step2/*Astrocyte*.regenie.raw_p >/dev/null 2>&1; then
        echo "  SKIP ${cohort} (already have Astrocyte .raw_p)"
        new=$(ls "${out}"/regenie_step2/*Astrocyte*.regenie.raw_p | head -1)
        old=$(ls "${SCC_RES}/${scc_results}"/regenie_step2/*Astrocyte*.regenie.raw_p 2>/dev/null | head -1 || true)
        if [ -n "${old}" ] && [ ! -f "${CMP}/${cohort}.tsv" ]; then
            "${RSCRIPT}" "${NF_DIR}/scripts/compare_raw_p.R" \
                --cohort "${cohort}" --cell_type "${CT}" \
                --new "${new}" --old "${old}" \
                --out "${CMP}/${cohort}.tsv" || true
        fi
        continue
    fi

    prefix=""
    if ! prefix=$(resolve_pgen_prefix "${cohort}" "${scc_results}" "${out}"); then
        echo "  FAIL ${cohort}: no readable .QC.final.pgen (SCC or ${out})"
        touch "${FAIL_FILE}"
        continue
    fi
    src="${prefix}.pgen"
    case "${src}" in
        "${SCC_RES}"*) echo "  pgen SCC  ${src}" ;;
        *)             echo "  pgen LOCAL ${src} (not June 9; SCC link dangling)" ;;
    esac

    cfgp=$(cfg_prefix "${cohort}")
    staged="${NF_DIR}/results/${cohort}/${cfgp}.pgen"
    if [ ! -s "${staged}" ]; then
        echo "  FAIL ${cohort}: missing login-node symlink ${staged}"
        echo "       run: bash ${NF_DIR}/validation/stage_scc_pgen_for_regenie.sh"
        touch "${FAIL_FILE}"
        continue
    fi

    # genotyping_qc_pipeline precomputed emit reads these from output_dir.
    qcs=$(qc_study_name "${cohort}")
    scc="${SCC_RES}/${scc_results}"
    # ln -sfn to the same path is an error (CMC_PENN local fallback).
    if [ -s "${scc}/pca.csv" ]; then
        if [ "$(readlink -f "${scc}/pca.csv")" != "$(readlink -f "${out}/pca.csv" 2>/dev/null || true)" ]; then
            ln -sfn "${scc}/pca.csv" "${out}/pca.csv"
        fi
    elif [ ! -s "${out}/pca.csv" ]; then
        echo "  FAIL ${cohort}: no pca.csv"
        touch "${FAIL_FILE}"
        continue
    fi
    prune_src=""
    if [ -s "${scc}/${qcs}.QC.prune.in" ]; then
        prune_src="${scc}/${qcs}.QC.prune.in"
    elif [ -s "${out}/${qcs}.QC.prune.in" ]; then
        prune_src="${out}/${qcs}.QC.prune.in"
    fi
    if [ -z "${prune_src}" ]; then
        echo "  FAIL ${cohort}: no ${qcs}.QC.prune.in"
        touch "${FAIL_FILE}"
        continue
    fi
    if [ "$(readlink -f "${prune_src}")" != "$(readlink -f "${out}/${qcs}.QC.prune.in" 2>/dev/null || true)" ]; then
        ln -sfn "${prune_src}" "${out}/${qcs}.QC.prune.in"
    fi

    extra=()
    [ -n "${biospecimen}" ] && [ "${biospecimen}" != "-" ] && extra+=(--biospecimen_file "${biospecimen}")
    [ -n "${clinical_metadata}" ] && [ "${clinical_metadata}" != "-" ] && extra+=(--clinical_metadata_file "${clinical_metadata}")
    extra+=(--work_dir "${work}/qc" --log_dir "${out}/logs")
    mkdir -p "${work}/qc" "${out}/logs"
    extra+=(--precomputed_pgen "${prefix}")
    extra+=(--precomputed_proportions "${scc}/cell_proportions.csv")
    extra+=(--precomputed_metadata "${scc}/combined_metrics.csv")

    if nextflow run "${NF_DIR}/combined_pipeline_v2.nf" \
        -c "${NF_DIR}/${config}" \
        -c "${NF_DIR}/nextflow.config.trillium" \
        -work-dir "${work}" \
        --count_matrix_file "${count_matrix}" \
        --metadata_file "${metadata}" \
        --output_dir "${out}" \
        --cell_type_list "${CT}" \
        --skip_metal true \
        --marker_file "${REFS_ROOT}/markers/new_MTGnCgG_lfct2.5_Publication.csv" \
        --hgnc_mapping_file "${REFS_ROOT}/markers/hgnc_complete_set.txt" \
        "${extra[@]+"${extra[@]}"}" \
        > "${LOG_ROOT}/validation/gwas_workers/${cohort}.log" 2>&1
    then
        echo "  OK   ${cohort}"
    else
        echo "  FAIL ${cohort}"
        touch "${FAIL_FILE}"
        continue
    fi

    new=$(ls "${out}"/regenie_step2/*Astrocyte*.regenie.raw_p 2>/dev/null | head -1 || true)
    old=$(ls "${SCC_RES}/${scc_results}"/regenie_step2/*Astrocyte*.regenie.raw_p 2>/dev/null | head -1 || true)
    if [ -n "${new}" ] && [ -n "${old}" ]; then
        "${RSCRIPT}" "${NF_DIR}/scripts/compare_raw_p.R" \
            --cohort "${cohort}" --cell_type "${CT}" \
            --new "${new}" --old "${old}" \
            --out "${CMP}/${cohort}.tsv" \
            >> "${LOG_ROOT}/validation/gwas_workers/${cohort}.log" 2>&1 || true
    fi
done < "${MAP_FILLED}"
if [ -d "${CMP}" ] && ls "${CMP}"/*.tsv >/dev/null 2>&1; then
    { head -1 "$(ls "${CMP}"/*.tsv | head -1)"; tail -n +2 "${CMP}"/*.tsv | grep -v '^cohort'; } \
      > "${CMP}/summary.tsv" || true
    echo ""
    echo "Per-cohort .raw_p comparison: ${CMP}/summary.tsv"
    cat "${CMP}/summary.tsv"
fi
[ ! -f "${FAIL_FILE}" ]
