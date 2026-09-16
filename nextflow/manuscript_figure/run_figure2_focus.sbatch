#!/bin/bash
#SBATCH --job-name=figure2_focus
#SBATCH --output=/scratch/zhoux156/logs/manuscript_figure/fig2focus_%j.log
#SBATCH --error=/scratch/zhoux156/logs/manuscript_figure/fig2focus_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute

set -euo pipefail
RSCRIPT="/home/zhoux156/miniforge3/envs/test/bin/Rscript"
SCRIPT_DIR="/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"
SCRATCH_OUT="/scratch/zhoux156/manuscript_figure"
FINAL_OUT="${SCRIPT_DIR}"

mkdir -p "${SCRATCH_OUT}"
# Merged figure needs the bulk long table in OUT_DIR
if [ -f "${FINAL_OUT}/bulk_celltype_proportions_long.tsv" ]; then
  cp -f "${FINAL_OUT}/bulk_celltype_proportions_long.tsv" "${SCRATCH_OUT}/"
elif [ -f "/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow/manuscript_figure/bulk_celltype_proportions_long.tsv" ]; then
  cp -f "/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow/manuscript_figure/bulk_celltype_proportions_long.tsv" "${SCRATCH_OUT}/"
fi

export FIG2_OUT_DIR="${SCRATCH_OUT}"
echo "Starting at $(date)"
echo "Rscript: ${RSCRIPT}"
echo "FIG2_OUT_DIR=${FIG2_OUT_DIR}"
"${RSCRIPT}" "${SCRIPT_DIR}/figure2_celltype_validation_focus.R"
"${RSCRIPT}" "${SCRIPT_DIR}/figure2_ctp_variation_merged.R"

echo "Copying outputs to ${FINAL_OUT}"
cp -f "${SCRATCH_OUT}/combined_bulk_snrna_paired.tsv" \
      "${SCRATCH_OUT}/figure2_celltype_accuracy.tsv" \
      "${SCRATCH_OUT}/figure2_celltype_validation_focus."{png,pdf,svg} \
      "${SCRATCH_OUT}/figure2_ctp_variation_merged."{png,pdf,svg} \
      "${FINAL_OUT}/"
echo "Finished at $(date)"
