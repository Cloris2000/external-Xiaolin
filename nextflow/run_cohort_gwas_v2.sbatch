#!/bin/bash
#SBATCH --job-name=cohort_gwas_v2
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/zhoux156/logs/cohort_gwas_v2/node_%j.out
#SBATCH --error=/scratch/zhoux156/logs/cohort_gwas_v2/node_%j.err

# cohort_gwas_v2: full per-cohort pipeline from the VCFs, on scratch.
#   Stage 1 RNA (PCA / tech-cov / MGP deconv, exclude_bio_from_tech_cov=true)
#   Stage 2 PHENO_PREP -> orientation gate -> Stage 3 VCF QC + PCA -> COV_PREP
#   -> REGENIE step 1 (hom_alt filter, the current default) + step 2, 19 cell types
#   -> comparison against the SCC run.
# No --precomputed_pgen / --precomputed_proportions / --june9_replay: everything
# is recomputed.  METAL is skipped here (run_meta_15cohorts_hg19_v2.sbatch).
#
# Cohorts are given as ARGUMENTS and run in parallel on this one node, each as
# its own Nextflow with 192/n CPUs and 700/n GB (nextflow.config.cohort_gwas_v2_node).
# Intended grouping (three concurrent jobs):
#   sbatch run_cohort_gwas_v2.sbatch GTEx_v10 GVEX ROSMAP                      # heavy WGS
#   sbatch run_cohort_gwas_v2.sbatch AMP_AD_Rush AMP_AD_Mayo Mayo MSBB CMC_MSSM
#   sbatch run_cohort_gwas_v2.sbatch NABEC ROSMAP_array CMC_PENN CMC_PITT NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M
# Re-submitting the same cohorts resumes (Nextflow -resume, same work dir).
#
# Walltime: 8 h default, sized from measured runs (REGENIE 19 CT on precomputed
# pgen took 155 min for GTEx, 122 for AMP_AD_Rush, <=90 for the rest; VCF->pgen
# QC scales from CMC_PENN's 25 GB / 20 min, giving ~3.6 h for GTEx's 622 GB).
# Slowest cohort per group: A ~6.5 h (GTEx), B ~4 h (AMP_AD_Rush), C ~1.5 h.
# Override at submit time for the lighter groups so they backfill sooner:
#   sbatch --time=6:00:00 run_cohort_gwas_v2.sbatch AMP_AD_Rush ...
# A job that hits the wall is resumable: resubmit the same cohorts (-resume).
#
# Arguments, not environment variables: ~/.bashrc wraps sbatch with
# --export=NONE, so `COHORTS=... sbatch ...` never reaches the job.
# For the same reason ORIENTATION_GATE is a flag: pass --orientation-gate=warn
# BEFORE the cohort names to let a cohort continue past a negative MGP
# orientation (default: that cohort stops before genotype QC and is reported).
# Trillium allocates whole nodes; do not add --mem.

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
source "${NF_DIR}/scripts/cohort_gwas_v2_lib.sh"
module load "${NEXTFLOW_MODULE}"

ORIENTATION_GATE=fail
while [ $# -gt 0 ]; do
    case "$1" in
        --orientation-gate=*) ORIENTATION_GATE="${1#*=}"; shift ;;
        --) shift; break ;;
        -*) echo "ERROR: unknown option $1" >&2; exit 1 ;;
        *)  break ;;
    esac
done
[ $# -gt 0 ] || { echo "ERROR: name the cohorts as arguments, e.g." >&2
                  echo "  sbatch run_cohort_gwas_v2.sbatch GTEx_v10 GVEX ROSMAP" >&2
                  echo "(not COHORTS=...: ~/.bashrc wraps sbatch with --export=NONE)" >&2; exit 1; }
COHORT_LIST=()
for a in "$@"; do IFS=', ' read -r -a _parts <<< "$a"; COHORT_LIST+=("${_parts[@]}"); done
n=${#COHORT_LIST[@]}

# Fail fast on a typo rather than silently running fewer cohorts than intended.
for c in "${COHORT_LIST[@]}"; do
    case " ${ALL_COHORTS[*]} " in *" ${c} "*) ;;
        *) echo "ERROR: unknown cohort '${c}'. Known: ${ALL_COHORTS[*]}" >&2; exit 1 ;;
    esac
done
export NF_SHARE_CPUS=$(( 192 / n ))
export NF_SHARE_MEM_GB=$(( 700 / n ))
# plink2 per-process cap: up to 22 chromosome conversions run at once per cohort.
PLINK_MEM_MB=$(( NF_SHARE_MEM_GB * 1024 / 26 )); [ "$PLINK_MEM_MB" -lt 4000 ] && PLINK_MEM_MB=4000; [ "$PLINK_MEM_MB" -gt 32000 ] && PLINK_MEM_MB=32000
OVERLAY="${NF_DIR}/nextflow.config.cohort_gwas_v2_node"

mkdir -p "${V2_RESULTS}" "${V2_COMPARE}" "${V2_WORK}" "${V2_LOGS}/workers"
echo "cohort_gwas_v2 node job ${SLURM_JOB_ID:-local} $(date -Is)"
echo "  cohorts (${n}): ${COHORT_LIST[*]}"
echo "  per-cohort share: ${NF_SHARE_CPUS} CPUs / ${NF_SHARE_MEM_GB} GB; plink2 cap ${PLINK_MEM_MB} MB"
echo "  results: ${V2_RESULTS}   compare: ${V2_COMPARE}   work: ${V2_WORK}"

FAIL_FILE="${V2_ROOT}/.fail_${SLURM_JOB_ID:-local}"; rm -f "${FAIL_FILE}"

run_one() {
    local cohort="$1" config="$2" count_matrix="$3" metadata="$4" biospecimen="$5" clinical="$6" scc_results="$7"
    local cfg="${NF_DIR}/${config}" out="${V2_RESULTS}/${cohort}" work="${V2_WORK}/${cohort}"
    local run="${SCRATCH_ROOT}/run_dirs/cohort_gwas_v2/${cohort}" wlog="${V2_LOGS}/workers/${cohort}.log"
    echo "======== START ${cohort} $(date -Is) ========"
    mkdir -p "${out}/logs" "${work}" "${run}"
    cd "${run}" || return 1
    export NXF_HOME="${SCRATCH_ROOT}/.nextflow_v2/${cohort}" TMPDIR="${SCRATCH_ROOT}/tmp/v2_${cohort}" NXF_ANSI_LOG=false
    mkdir -p "${NXF_HOME}" "${TMPDIR}"

    local args=(
        -c "${cfg}" -c "${NF_DIR}/nextflow.config.trillium" -c "${OVERLAY}"
        -work-dir "${work}" -resume
        --output_dir "${out}" --work_dir "${work}/qc" --log_dir "${out}/logs"
        --count_matrix_file "${count_matrix}" --metadata_file "${metadata}"
        --marker_file "${REFS_ROOT}/markers/new_MTGnCgG_lfct2.5_Publication.csv"
        --hgnc_mapping_file "${REFS_ROOT}/markers/hgnc_complete_set.txt"
        --exclude_bio_from_tech_cov true
        --plink_memory "${PLINK_MEM_MB}"
        --skip_metal true
    )
    [ "${biospecimen}" != "-" ] && args+=(--biospecimen_file "${biospecimen}")
    [ "${clinical}"    != "-" ] && args+=(--clinical_metadata_file "${clinical}")
    local stk; stk=$(resolve_samples_to_keep "${cohort}" "${cfg}")
    [ -n "${stk}" ] && args+=(--samples_to_keep "${stk}") && echo "  ${cohort}: samples_to_keep ${stk}"
    genotype_override_args "${cohort}" args
    echo "  ${cohort}: genotype input $(dirname "$(resolve_vcf "${cfg}" 1 "${cohort}")")"

    # Pass 1: RNA + phenotypes only, so the MGP orientation can be checked before
    # any genotype work starts.
    echo "  ${cohort}: pass 1 (deconv + phenotypes)"
    if ! nextflow run "${NF_DIR}/combined_pipeline_v2.nf" "${args[@]}" --stop_after pheno > "${wlog}" 2>&1; then
        echo "  FAIL ${cohort}: pass 1 (see ${wlog})"; touch "${FAIL_FILE}"; return 1
    fi
    "${RSCRIPT}" "${NF_DIR}/scripts/check_deconv_orientation.R" --cohort "${cohort}" \
        --proportions "${out}/cell_proportions.csv" --corrected "${out}/corrected_data.RData" \
        --marker_file "${REFS_ROOT}/markers/new_MTGnCgG_lfct2.5_Publication.csv" \
        --out "${out}/deconv_orientation.tsv" >> "${wlog}" 2>&1
    local rc=$?
    if [ "${rc}" -ne 0 ]; then
        grep "ORIENTATION" "${wlog}" | tail -1
        if [ "${ORIENTATION_GATE}" != "warn" ]; then
            echo "  STOP ${cohort}: orientation gate failed; GWAS not run (ORIENTATION_GATE=warn to override)"
            touch "${FAIL_FILE}"; return 1
        fi
        echo "  WARN ${cohort}: continuing past orientation failure (ORIENTATION_GATE=warn)"
    fi

    # Pass 2: full pipeline; Stage 1/2 tasks are cached by -resume.
    echo "  ${cohort}: pass 2 (genotype QC + PCA + REGENIE) $(date -Is)"
    if ! nextflow run "${NF_DIR}/combined_pipeline_v2.nf" "${args[@]}" >> "${wlog}" 2>&1; then
        echo "  FAIL ${cohort}: pass 2 (see ${wlog})"; touch "${FAIL_FILE}"; return 1
    fi
    local nraw; nraw=$(ls "${out}"/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
    echo "  OK   ${cohort} $(date -Is)  raw_p files: ${nraw}/19"
    [ "${nraw}" -eq 19 ] || { echo "  FAIL ${cohort}: expected 19 raw_p"; touch "${FAIL_FILE}"; }

    echo "  ${cohort}: comparison vs SCC"
    "${NF_DIR}/scripts/compare_cohort_to_scc.sh" "${cohort}" "${out}" "${SCC_RES}/${scc_results}" \
        "${V2_COMPARE}/${cohort}" >> "${wlog}" 2>&1
    echo "  ${cohort}: flags ($(wc -l < "${V2_COMPARE}/${cohort}/flags.txt")):"; sed 's/^/      /' "${V2_COMPARE}/${cohort}/flags.txt"
    return 0
}

pids=()
while IFS=$'\t' read -r cohort config count_matrix metadata biospecimen clinical vcf_dir vcf_glob scc_results; do
    case " ${COHORT_LIST[*]} " in *" ${cohort} "*) ;; *) continue ;; esac
    run_one "${cohort}" "${config}" "${count_matrix}" "${metadata}" "${biospecimen}" "${clinical}" "${scc_results}" &
    pids+=($!); echo "  launched ${cohort} pid $!"
done < <(map_rows)
[ ${#pids[@]} -eq "${n}" ] || echo "WARNING: ${#pids[@]} of ${n} requested cohorts found in ${V2_MAP}"

rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "All cohort pipelines finished $(date -Is) rc=${rc}"

# Node-level summary of the comparisons that ran.
SUM="${V2_COMPARE}/summary_gwas_all_ct.tsv"
first=$(ls "${V2_COMPARE}"/*/gwas_all_ct.tsv 2>/dev/null | head -1)
if [ -n "${first}" ]; then
    { head -1 "${first}"; for f in "${V2_COMPARE}"/*/gwas_all_ct.tsv; do tail -n +2 "$f"; done; } > "${SUM}"
    echo "Per-cohort x cell-type comparison: ${SUM}"
    awk -F'\t' 'NR>1{n++; if($6>=0.999) e++; s+=$6} END{printf "  %d cohort x cell-type pairs; %d exact (r_beta>=0.999); mean r_beta=%.3f\n", n, e, s/n}' "${SUM}"
fi
[ ! -f "${FAIL_FILE}" ] && [ "${rc}" -eq 0 ]
