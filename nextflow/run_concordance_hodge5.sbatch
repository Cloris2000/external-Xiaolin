#!/bin/bash
# Recompute bulk<->sn direction concordance using the 5-cohort sn Hodge meta.
# Bulk: 15-cohort meta.
# sn:  ROSMAP_Green + PsychAD_HBCC + PsychAD_MSSM + PsychAD_RADC(Mathys) + Ruz_MSSM (Ruzicka) Hodge meta.
#
# This is the hodge5 counterpart of run_concordance_hodge3.sbatch. It leaves the
# 3-cohort hodge3 outputs untouched and writes to sn_bulk_meta_similarity_hodge5.
#
#SBATCH --job-name=conc_hodge5
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/concordance_hodge5_%j.out
#SBATCH --error=logs/concordance_hodge5_%j.err
#SBATCH --account=rrg-shreejoy

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd "$NF_DIR"
mkdir -p logs

OUTDIR=${NF_DIR}/results/sn_bulk_meta_similarity_hodge5/top_hits

conda run -n test Rscript ${NF_DIR}/scripts/bulk_top_hits_sn_direction.R \
  --bulk_dir ${NF_DIR}/results/meta_analysis_15cohorts \
  --sn_dir   ${NF_DIR}/results/meta_analysis_sn_hodge5 \
  --outdir   ${OUTDIR} \
  --match_by_name TRUE \
  --sugg_thresh 1e-5 \
  --clump_window_kb 500

echo "Concordance tables written to ${OUTDIR}"
echo "Finished: $(date)"
