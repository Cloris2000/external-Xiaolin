#!/bin/bash
# SLURM array job: one task per cell type, all run in parallel.
# Each task loads a single 1.6 GB .tbl file and produces an annotated
# Manhattan + QQ plot using topr::manhattanExtra().
#
# Submit:  sbatch run_plot_annotated_array.sbatch
#
#SBATCH --job-name=plot_annotated
#SBATCH --partition=medium
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=0-18          # 19 cell types (indices 0..18)
#SBATCH --output=logs/meta_15cohorts/plot_annotated_%A_%a.out
#SBATCH --error=logs/meta_15cohorts/plot_annotated_%A_%a.err

# ---- cell type list (must match exactly the .tbl file prefixes) -------------
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

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
RESULTS_DIR=${NF_DIR}/results/meta_analysis_15cohorts
OUTPUT_DIR=${RESULTS_DIR}/plots/annotated

echo "Array task ${SLURM_ARRAY_TASK_ID}: cell type = ${CT}"
echo "Started: $(date)"

# Use the test conda env which has topr / manhattanExtra
conda run -n test Rscript ${NF_DIR}/scripts/plot_meta_annotated.R \
    --results_dir ${RESULTS_DIR} \
    --output_dir  ${OUTPUT_DIR} \
    --cell_type   ${CT} \
    --gw_thresh   5e-8 \
    --sugg_thresh 1e-5

echo "Finished: $(date)"
