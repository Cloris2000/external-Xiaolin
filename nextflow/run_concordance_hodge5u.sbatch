#!/bin/bash
# Recompute bulk<->sn direction concordance using the 5-cohort sn Hodge meta,
# de-duplicated (Mathys restricted to donors NOT shared with ROSMAP_Green).
# Bulk: 15-cohort meta on scratch.
# sn:   SCC tree (ROSMAP_Green + PsychAD_HBCC + PsychAD_MSSM + Ruz_MSSM +
#       ROSMAP_Mathys_unique).
#
#SBATCH --job-name=conc_hodge5u
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/zhoux156/logs/concordance/concordance_hodge5u_%j.out
#SBATCH --error=/scratch/zhoux156/logs/concordance/concordance_hodge5u_%j.err

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

BULK_DIR="${BULK_DIR:-${RESULTS_ROOT}/meta_analysis_15cohorts}"
SN_DIR="${SN_DIR:-${SCC_DIR}/results/meta_analysis_sn_hodge5u}"
OUTDIR="${OUTDIR:-${RESULTS_ROOT}/sn_bulk_meta_similarity_hodge5u/top_hits}"
mkdir -p "${OUTDIR}"

echo "Bulk : ${BULK_DIR}"
echo "sn   : ${SN_DIR}"
echo "Out  : ${OUTDIR}"

"${RSCRIPT}" "${NF_DIR}/scripts/bulk_top_hits_sn_direction.R" \
  --bulk_dir "${BULK_DIR}" \
  --sn_dir   "${SN_DIR}" \
  --outdir   "${OUTDIR}" \
  --match_by_name TRUE \
  --sugg_thresh 1e-5 \
  --clump_window_kb 500

echo "Concordance tables written to ${OUTDIR}"
echo "Finished: $(date)"
