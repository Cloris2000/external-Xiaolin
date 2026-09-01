#!/bin/bash
#SBATCH --job-name=meta_15cohorts_ldsc
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_15cohorts/meta_15cohorts_ldsc_%j.out
#SBATCH --error=logs/meta_15cohorts/meta_15cohorts_ldsc_%j.err

# LDSC h2 follow-up for the 15-cohort meta-analysis.
# Runs in the same run/work dirs as run_meta_15cohorts.sbatch so that
# -resume reuses cached harmonization + METAL outputs and only the
# LDSC stages (PREPARE_LDSC_SUMSTATS, RUN_LDSC_H2, SUMMARIZE_LDSC_RESULTS) run.

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/meta_15cohorts

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_15cohorts_ldsc \
    -work-dir ${NF_DIR}/work/meta_15cohorts \
    -resume

echo ""
echo "15-cohort LDSC h2 complete."
echo "Summary: ${NF_DIR}/results/meta_analysis_15cohorts/ldsc/summary/ldsc_h2_summary.tsv"
