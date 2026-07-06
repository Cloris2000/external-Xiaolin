#!/bin/bash
#SBATCH --job-name=regen_eigencor
#SBATCH --partition=medium
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/regen_eigencor_%j.out
#SBATCH --error=logs/regen_eigencor_%j.err

set -euo pipefail

NF_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
mkdir -p "${NF_DIR}/logs"

echo "============================================"
echo "Regenerating before/after eigencor plots"
echo "Job ID : ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "============================================"

source /nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh
conda activate test

Rscript "${NF_DIR}/scripts/regen_eigencor_plots.R" \
    --results_dir "${NF_DIR}/results" \
    --cohorts "ROSMAP"

echo ""
echo "============================================"
echo "Done: $(date)"
echo "Plots written to results/<cohort>/:"
echo "  eigencor_plot_before_correction.png"
echo "  eigencor_plot_removedBatchEff_cov.png  (overwritten)"
echo "============================================"
