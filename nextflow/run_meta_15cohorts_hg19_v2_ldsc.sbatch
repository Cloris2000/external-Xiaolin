#!/bin/bash
#SBATCH --job-name=meta_v2_ldsc
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/meta_v2_ldsc_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/meta_v2_ldsc_%j.err

# LDSC h2 per cell type for the corrected 15-cohort meta (v2).
# Runs in the SAME run/work dirs as run_meta_15cohorts_hg19_v2.sbatch so -resume
# reuses the cached METAL + harmonized outputs and only the LDSC stages execute.
# Output: ${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2/ldsc/summary/ldsc_h2_summary.tsv
# (feeds the h2 panel of figure2_ctp_variation_merged.R via MF_H2_FILE).
#
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"

WORK_DIR="${WORK_ROOT}/meta_15cohorts_hg19_v2"          # must match the v2 meta run
RUN_DIR="${SCRATCH_ROOT}/run_dirs/meta_15cohorts_hg19_v2"
cd "${RUN_DIR}"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/nextflow.config.standalone_meta_15cohorts_hg19_v2_ldsc" \
    -work-dir "${WORK_DIR}" \
    -resume

echo ""
echo "v2 LDSC h2 complete: ${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2/ldsc/summary/ldsc_h2_summary.tsv"
