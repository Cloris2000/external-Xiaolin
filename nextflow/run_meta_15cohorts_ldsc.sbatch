#!/bin/bash
#SBATCH --job-name=meta_15cohorts_ldsc
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=12:00:00
#SBATCH --output=logs/meta_15cohorts/meta_15cohorts_ldsc_%j.out
#SBATCH --error=logs/meta_15cohorts/meta_15cohorts_ldsc_%j.err

# LDSC h2 follow-up for the 15-cohort meta-analysis.
# Runs in the same run/work dirs as run_meta_15cohorts.sbatch so that -resume
# reuses the cached METAL outputs and only the LDSC stages run
# (PREPARE_LDSC_SUMSTATS, RUN_LDSC_H2, SUMMARIZE_LDSC_RESULTS).
#
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail

NF_DIR=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow
WORK_DIR=/scratch/zhoux156/nf_work/meta_15cohorts

# Must match run_meta_15cohorts.sbatch for -resume to hit the cache.
module load nextflow/25.10.2

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate test

cd "${NF_DIR}/run_dirs/meta_15cohorts"

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_15cohorts_ldsc \
    -work-dir ${WORK_DIR} \
    -resume

echo ""
echo "15-cohort LDSC h2 complete."
echo "Summary: ${NF_DIR}/results/meta_analysis_15cohorts/ldsc/summary/ldsc_h2_summary.tsv"
