#!/bin/bash
#SBATCH --job-name=figure2_focus
#SBATCH --output=/scratch/zhoux156/logs/manuscript_figure/fig2focus_%j.log
#SBATCH --error=/scratch/zhoux156/logs/manuscript_figure/fig2focus_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=rrg-shreejoy

set -euo pipefail
SCRIPT_DIR="/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"
echo "Starting at $(date)"
Rscript "${SCRIPT_DIR}/figure2_celltype_validation_focus.R"
echo "Finished at $(date)"
