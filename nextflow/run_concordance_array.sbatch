#!/bin/bash
# SLURM array job: one task per bulk cell type (19 total).
# Each task reads the bulk tbl (~1.6 GB), filters to P < 1e-5, clumps
# into independent loci, looks each lead SNP up in the matched sn
# meta-analysis tbl, and writes a concordance TSV.
#
# Submit:  sbatch run_concordance_array.sbatch
# After completion, run:  conda run -n test Rscript scripts/plot_concordance_barplot.R
#
#SBATCH --job-name=concordance
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-18
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/concordance_%A_%a.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/concordance_%A_%a.err
#SBATCH --account=rrg-shreejoy

CELL_TYPES=(
    Astrocyte
    Endothelial
    IT
    L4.IT
    L5.6.IT.Car3
    L5.6.NP
    L5.ET
    L6.CT
    L6b
    LAMP5
    Microglia
    Oligodendrocyte
    OPC
    PAX6
    Pericyte
    PVALB
    SST
    VIP
    VLMC
)

CT=${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

echo "Array task ${SLURM_ARRAY_TASK_ID}: cell type = ${CT}"
echo "Started: $(date)"

conda run -n test Rscript ${NF_DIR}/scripts/bulk_sn_concordance_per_ct.R \
    --celltype "${CT}"

echo "Finished: $(date)"
