#!/bin/bash
#SBATCH --job-name=meta_15cohorts_ldsc
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_ldsc_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_ldsc_%j.err

# LDSC h2 follow-up for the 15-cohort meta-analysis.
# Runs in the same run/work dirs as run_meta_15cohorts.sbatch so that -resume
# reuses the cached METAL outputs and only the LDSC stages run
# (PREPARE_LDSC_SUMSTATS, RUN_LDSC_H2, SUMMARIZE_LDSC_RESULTS).
#
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
WORK_DIR="${WORK_ROOT}/meta_15cohorts"

# Must match run_meta_15cohorts.sbatch for -resume to hit the cache.
module load "${NEXTFLOW_MODULE}"

cd "${NF_DIR}/run_dirs/meta_15cohorts"

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_15cohorts_ldsc \
    -work-dir ${WORK_DIR} \
    -resume

echo ""
echo "15-cohort LDSC h2 complete."
echo "Summary: ${NF_DIR}/results/meta_analysis_15cohorts/ldsc/summary/ldsc_h2_summary.tsv"
