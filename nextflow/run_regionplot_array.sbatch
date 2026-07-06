#!/bin/bash
# SLURM array job: one task per cell-type locus of interest.
# Each task extracts its region from the .tbl file and draws a region plot
# colored by direction of effect (positive β = red, negative β = blue).
#
# Submit:  sbatch run_regionplot_array.sbatch
#
#SBATCH --job-name=regionplot
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=0-1          # 2 tasks: Oligodendrocyte + VIP
#SBATCH --output=logs/meta_15cohorts/regionplot_%A_%a.out
#SBATCH --error=logs/meta_15cohorts/regionplot_%A_%a.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
RESULTS_DIR=${NF_DIR}/results/meta_analysis_15cohorts
OUTPUT_DIR=${RESULTS_DIR}/plots/regionplot

# ---- Per-task config: CELL_TYPE, LEAD_SNP -----------------------------------
CELL_TYPES=( "Oligodendrocyte"      "VIP"                  )
LEAD_SNPS=(  "chr7:12252119:G:GT"   "chr7:12284430:T:TA"   )

CT=${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}
LEAD=${LEAD_SNPS[$SLURM_ARRAY_TASK_ID]}

# Find the matching .tbl symlink (without trailing 1)
TBL=$(ls ${RESULTS_DIR}/${CT}_meta_analysis_*.tbl 2>/dev/null | grep -v '1\.tbl$' | head -1)

if [ -z "$TBL" ]; then
    echo "ERROR: No .tbl file found for ${CT}"
    exit 1
fi

mkdir -p ${OUTPUT_DIR}
echo "Task ${SLURM_ARRAY_TASK_ID}: ${CT}  lead=${LEAD}"
echo "TBL: ${TBL}"
echo "Started: $(date)"

conda run -n test Rscript ${NF_DIR}/scripts/plot_regionplot.R \
    --tbl_file    "${TBL}" \
    --cell_type   "${CT}" \
    --lead_snp    "${LEAD}" \
    --region_size 1000000 \
    --output_dir  "${OUTPUT_DIR}" \
    --gw_thresh   5e-8

echo "Finished: $(date)"
