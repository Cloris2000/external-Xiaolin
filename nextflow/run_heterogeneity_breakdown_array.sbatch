#!/bin/bash
#SBATCH --job-name=het_ct
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/zhoux156/logs/het_breakdown/het_breakdown_%j.out
#SBATCH --error=/scratch/zhoux156/logs/het_breakdown/het_breakdown_%j.err

# Heterogeneity breakdown for all 19 cell types.
#
# This was a 19-way SLURM array. Trillium allocates whole nodes, so an array
# would have burned 19 nodes on single-threaded jobs and queued badly. All 19
# now run as parallel workers inside one node, then the merge step runs once.
#
# Trillium allocates whole nodes; do not add --mem.
#
#   sbatch run_heterogeneity_breakdown_array.sbatch                     # Track A results
#   META_DIR=... OUT_DIR=... sbatch run_heterogeneity_breakdown_array.sbatch

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

# Defaults to the replication meta; point META_DIR at the hg38fix tree for Track B.
META_DIR="${META_DIR:-${RESULTS_ROOT}/meta_analysis_15cohorts}"
OUT_DIR="${OUT_DIR:-${META_DIR}/heterogeneity_breakdown}"
WORKER_LOGS="${LOG_ROOT}/het_breakdown/workers_${SLURM_JOB_ID:-manual}"
mkdir -p "${OUT_DIR}" "${WORKER_LOGS}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

COHORTS=(
  AMP_AD_Mayo AMP_AD_Rush
  CMC_MSSM CMC_PENN CMC_PITT
  GTEx_v10 GVEX MSBB Mayo NABEC
  NIMH_HBCC_1M NIMH_HBCC_Omni5M NIMH_HBCC_h650
  ROSMAP ROSMAP_array
)

# The per-cohort REGENIE output lives in the SCC tree copied to Trillium, not
# under NF_DIR, which only holds meta-analysis results.
# Per-cohort summary statistics.  Default: raw REGENIE output in the SCC tree.
# SUMSTAT_PATTERN overrides the whole path template ({cohort}, {cell_type}
# placeholders); the v2 run points it at the harmonized (hg19-lifted) files so
# the four hg38 cohorts match the meta markers instead of silently dropping out.
SUMSTAT_DIR="${SUMSTAT_DIR:-${SCC_DIR}/results}"
# Built in two steps on purpose: inside ${VAR:-default} bash ends the expansion at
# the first '}', so a default containing {cohort} would be truncated and the rest
# (with its trailing '}') appended literally -> Python format error "Single '}'".
DEFAULT_SUMSTAT_PATTERN="${SUMSTAT_DIR}/{cohort}/regenie_step2/{cohort}_{cell_type}_step2.regenie.raw_p"
SUMSTAT_PATTERN="${SUMSTAT_PATTERN:-$DEFAULT_SUMSTAT_PATTERN}"

echo "Meta tables : ${META_DIR}"
echo "Output      : ${OUT_DIR}"
echo "Cohort data : ${SUMSTAT_DIR}"
echo "Running ${#CELL_TYPES[@]} cell types across ${#CELL_TYPES[@]} parallel workers"

# Bash arrays do not survive `export -f` into the xargs subshells, so the cohort
# list crosses as a space-separated string and is split back apart here.
run_one() {
    local ct="$1"
    local -a cohorts
    read -r -a cohorts <<< "${COHORTS_STR}"
    python3 "${NF_DIR}/scripts/heterogeneity_breakdown.py" \
        --cell-type "${ct}" \
        --meta-tbl-dir "${META_DIR}" \
        --cohort-sumstat-pattern "${SUMSTAT_PATTERN}" \
        --cohort-metadata "${NF_DIR}/docs/meta_cohort_metadata.tsv" \
        --output-dir "${OUT_DIR}" \
        --cohorts "${cohorts[@]}" \
        --sumstat-results-dir "${SUMSTAT_DIR}" \
        > "${WORKER_LOGS}/${ct}.log" 2>&1 \
        && echo "  OK   ${ct}" \
        || { echo "  FAIL ${ct} (see ${WORKER_LOGS}/${ct}.log)"; return 1; }
}
export -f run_one
export NF_DIR META_DIR OUT_DIR SUMSTAT_DIR WORKER_LOGS
export COHORTS_STR="${COHORTS[*]}"

# xargs returns 123 if any worker failed; keep going so one bad cell type does
# not hide the rest, then report.
printf '%s\n' "${CELL_TYPES[@]}" \
  | xargs -P "${#CELL_TYPES[@]}" -I{} bash -c 'run_one "$@"' _ {} \
  && rc=0 || rc=$?

echo ""
if [ "$rc" -ne 0 ]; then
    echo "ERROR: at least one cell type failed; not merging. Logs: ${WORKER_LOGS}"
    exit 1
fi

echo "All cell types done. Merging..."
python3 "${NF_DIR}/scripts/heterogeneity_breakdown.py" \
    --merge-only \
    --output-dir "${OUT_DIR}"

echo ""
echo "Heterogeneity breakdown complete: ${OUT_DIR}"
