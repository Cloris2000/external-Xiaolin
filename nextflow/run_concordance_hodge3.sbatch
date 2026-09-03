#!/bin/bash
# Recompute bulk↔sn direction concordance using the newest sn Hodge3 meta.
# Bulk: 15-cohort meta; sn: ROSMAP_Green + PsychAD_HBCC + PsychAD_MSSM Hodge meta.
#
#SBATCH --job-name=conc_hodge3
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/scratch/zhoux156/logs/concordance_hodge3_%j.out
#SBATCH --error=/scratch/zhoux156/logs/concordance_hodge3_%j.err
#SBATCH --account=rrg-shreejoy

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd "$NF_DIR"
mkdir -p logs

OUTDIR=${NF_DIR}/results/sn_bulk_meta_similarity_hodge3/top_hits

conda run -n test Rscript ${NF_DIR}/scripts/bulk_top_hits_sn_direction.R \
  --bulk_dir ${NF_DIR}/results/meta_analysis_15cohorts \
  --sn_dir   ${NF_DIR}/results/meta_analysis_sn_hodge3 \
  --outdir   ${OUTDIR} \
  --match_by_name TRUE \
  --sugg_thresh 1e-5 \
  --clump_window_kb 500

echo "Concordance tables written to ${OUTDIR}"
echo "Finished: $(date)"
