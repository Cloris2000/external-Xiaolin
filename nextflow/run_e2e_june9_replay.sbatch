#!/bin/bash
#SBATCH --job-name=e2e_june9
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/zhoux156/logs/validation/e2e_june9_%j.out
#SBATCH --error=/scratch/zhoux156/logs/validation/e2e_june9_%j.err

# End-to-end June 9 replay: RNA deconv → PHENO_PREP → Astrocyte REGENIE.
# Does NOT feed SCC phenotypes_RINT.txt.
# Opt-in pins only (defaults elsewhere unchanged):
#   --june9_replay true
#   --regenie_step1_snp_filter mac1
#   --pheno_tissue_filter from the cohort config
#   --pheno_fid_specimen_map if validation/june9_fid_specimen_maps/<cohort>.tsv exists
# Stage 1 tech cov stays at the pipeline default (auto_top_n). Do not pass
# fixed_list — that emptied the cov list and broke cohorts that already matched.
# Default: skip VCF QC and use SCC pgen (PENN/PITT SCC pgen is dangling).
# RUN_VCF_QC=1: Stage 3 from config vcf_pattern (June 9 joint CMC.QC for PENN/PITT).
#
# All 15 cohorts run as parallel Nextflow instances on this one node
# (nextflow.config.e2e_parallel_cohorts caps each instance).

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"

MAP="${NF_DIR}/validation/cohort_trillium_paths.tsv"
CT="${CELL_TYPE:-Astrocyte}"
OUT_ROOT="${RESULTS_ROOT}/validation/e2e_june9"
CMP_DECONV="${RESULTS_ROOT}/validation/e2e_june9_deconv_compare"
CMP_GWAS="${RESULTS_ROOT}/validation/e2e_june9_gwas_compare"
SCC_RES="${SCC_DIR}/results"
MAP_DIR="${NF_DIR}/validation/june9_fid_specimen_maps"
WORK_NAME="${WORK_NAME:-validation_e2e_june9_auto}"
mkdir -p "${OUT_ROOT}" "${CMP_DECONV}" "${CMP_GWAS}" "${LOG_ROOT}/validation/e2e_workers"

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
    for f in "${SCC_RES}/${scc_results}"/*.QC.final.pgen; do
        [ -s "${f}" ] && { echo "${f%.pgen}"; return 0; }
    done
    local fallback="/scratch/zhoux156/results/validation/gwas_Astrocyte/${cohort}"
    for f in "${fallback}"/*.QC.final.pgen "${out}"/*.QC.final.pgen; do
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
    local biospecimen="$5" clinical_metadata="$6" vcf_dir="$7" vcf_glob="$8" scc_results="$9"
    local work run out prefix cfgp qcs scc extra tissue
    echo "======== START ${cohort} $(date -Is) ========"
    work="${WORK_ROOT}/${WORK_NAME}/${cohort}"
    run="${SCRATCH_ROOT}/run_dirs/${WORK_NAME}/${cohort}"
    out="${OUT_ROOT}/${cohort}"
    mkdir -p "${work}" "${run}" "${out}/logs"
    cd "${run}"
    export NXF_HOME="${SCRATCH_ROOT}/.nextflow_e2e/${cohort}"
    export TMPDIR="${SCRATCH_ROOT}/tmp/e2e_${cohort}"
    export NXF_ANSI_LOG=false
    mkdir -p "${NXF_HOME}" "${TMPDIR}"

    qcs=$(qc_study_name "${cohort}")
    scc="${SCC_RES}/${scc_results}"
    extra=()
    # SCC PENN/PITT pgen is a dangling CAMH link; QC from the June 9 joint VCF.
    if [ "${RUN_VCF_QC:-0}" = "1" ] || [ "${cohort}" = "CMC_PENN" ] || [ "${cohort}" = "CMC_PITT" ]; then
        echo "  VCF QC from config vcf_pattern (no precomputed pgen)"
        rm -f "${out}/${qcs}.QC.prune.in" "${out}/pca.csv"
        rm -rf "${out}/regenie_step1" "${out}/regenie_step2"
        if [ -n "${vcf_dir}" ] && [ "${vcf_dir}" != "-" ] && [ -n "${vcf_glob}" ] && [ "${vcf_glob}" != "-" ]; then
            extra+=(--vcf_dir "${vcf_dir}")
            extra+=(--vcf_pattern "${vcf_dir}/${vcf_glob//\*/\{chrom\}}")
        fi
    else
        prefix=""
        if ! prefix=$(resolve_pgen_prefix "${cohort}" "${scc_results}" "${out}"); then
            echo "  FAIL ${cohort}: no readable pgen"
            touch "${FAIL_FILE}"
            return 1
        fi
        cfgp=$(cfg_prefix "${cohort}")
        if [ ! -s "${NF_DIR}/results/${cohort}/${cfgp}.pgen" ]; then
            echo "  FAIL ${cohort}: run validation/stage_scc_pgen_for_regenie.sh first"
            touch "${FAIL_FILE}"
            return 1
        fi
        [ -s "${scc}/pca.csv" ] && ln -sfn "${scc}/pca.csv" "${out}/pca.csv"
        prune_src=""
        if [ -s "${scc}/${qcs}.QC.prune.in" ]; then
            prune_src="${scc}/${qcs}.QC.prune.in"
        else
            local_prune="/scratch/zhoux156/results/validation/gwas_Astrocyte/${cohort}/${qcs}.QC.prune.in"
            if [ -s "${local_prune}" ]; then
                prune_src="${local_prune}"
            elif [ -s "${out}/${qcs}.QC.prune.in" ]; then
                prune_src="${out}/${qcs}.QC.prune.in"
            fi
        fi
        if [ -n "${prune_src}" ]; then
            [ "$(readlink -f "${prune_src}")" != "$(readlink -f "${out}/${qcs}.QC.prune.in" 2>/dev/null || true)" ] && \
                ln -sfn "${prune_src}" "${out}/${qcs}.QC.prune.in"
        fi
        extra+=(--precomputed_pgen "${prefix}")
    fi
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

    nf_extra_c=()
    if [ "${RUN_VCF_QC:-0}" = "1" ] || [ "${cohort}" = "CMC_PENN" ] || [ "${cohort}" = "CMC_PITT" ]; then
        nf_extra_c+=(-c "${NF_DIR}/nextflow.config.e2e_vcf_qc")
    fi
    if nextflow run "${NF_DIR}/combined_pipeline_v2.nf" \
        -resume \
        -c "${NF_DIR}/${config}" \
        -c "${NF_DIR}/nextflow.config.trillium" \
        -c "${NF_DIR}/nextflow.config.e2e_parallel_cohorts" \
        "${nf_extra_c[@]+"${nf_extra_c[@]}"}" \
        -work-dir "${work}" \
        --count_matrix_file "${count_matrix}" \
        --metadata_file "${metadata}" \
        --output_dir "${out}" \
        --cell_type_list "${CT}" \
        --skip_metal true \
        --marker_file "${REFS_ROOT}/markers/new_MTGnCgG_lfct2.5_Publication.csv" \
        --hgnc_mapping_file "${REFS_ROOT}/markers/hgnc_complete_set.txt" \
        "${extra[@]+"${extra[@]}"}" \
        > "${LOG_ROOT}/validation/e2e_workers/${cohort}.log" 2>&1
    then
        echo "  OK   ${cohort} $(date -Is)"
    else
        echo "  FAIL ${cohort} $(date -Is)"
        touch "${FAIL_FILE}"
        return 1
    fi

    if [ -s "${out}/cell_proportions.csv" ]; then
        "${RSCRIPT}" "${NF_DIR}/scripts/compare_deconv_to_scc.R" \
            --cohort "${cohort}" \
            --new_prop "${out}/cell_proportions.csv" \
            --old_prop "${scc}/cell_proportions.csv" \
            --new_pheno "${out}/phenotypes_RINT.txt" \
            --old_pheno "${scc}/phenotypes_RINT.txt" \
            --out "${CMP_DECONV}/${cohort}.tsv" \
            >> "${LOG_ROOT}/validation/e2e_workers/${cohort}.log" 2>&1 || true
    fi
    new=$(ls "${out}"/regenie_step2/*Astrocyte*.regenie.raw_p 2>/dev/null | head -1 || true)
    old=$(ls "${scc}"/regenie_step2/*Astrocyte*.regenie.raw_p 2>/dev/null | head -1 || true)
    if [ -n "${new}" ] && [ -n "${old}" ]; then
        "${RSCRIPT}" "${NF_DIR}/scripts/compare_raw_p.R" \
            --cohort "${cohort}" --cell_type "${CT}" \
            --new "${new}" --old "${old}" \
            --out "${CMP_GWAS}/${cohort}.tsv" \
            >> "${LOG_ROOT}/validation/e2e_workers/${cohort}.log" 2>&1 || true
    fi
    return 0
}

FAIL_FILE="${OUT_ROOT}/.fail"
rm -f "${FAIL_FILE}"
MAP_FILLED="${OUT_ROOT}/.cohort_map_filled.tsv"
sed -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' "${MAP}" > "${MAP_FILLED}"

COHORT_FILTER="${COHORTS:-}"
COHORT_FILTER="${COHORT_FILTER//,/ }"
echo "Launching cohorts in parallel $(date -Is)${COHORT_FILTER:+ (filter: ${COHORT_FILTER})}"
pids=()
while IFS=$'\t' read -r cohort config count_matrix metadata biospecimen clinical_metadata vcf_dir vcf_glob scc_results; do
    [ "${cohort}" = "cohort" ] && continue
    if [ -n "${COHORT_FILTER}" ]; then
        case " ${COHORT_FILTER} " in
            *" ${cohort} "*) ;;
            *) continue ;;
        esac
    fi
    run_one "${cohort}" "${config}" "${count_matrix}" "${metadata}" \
        "${biospecimen}" "${clinical_metadata}" "${vcf_dir}" "${vcf_glob}" "${scc_results}" &
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

summarize() {
    local d="$1"
    [ -d "${d}" ] || return 0
    ls "${d}"/*.tsv >/dev/null 2>&1 || return 0
    { head -1 "$(ls "${d}"/*.tsv | head -1)"; tail -n +2 "${d}"/*.tsv | grep -v '^cohort'; } > "${d}/summary.tsv" || true
    echo "Summary: ${d}/summary.tsv"
    cat "${d}/summary.tsv"
}
echo ""
summarize "${CMP_DECONV}"
summarize "${CMP_GWAS}"
[ ! -f "${FAIL_FILE}" ]
