#!/bin/bash
#SBATCH --job-name=figure2_focus
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/fig2focus_%j.log
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/fig2focus_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail
SCRIPT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure"
echo "Starting at $(date)"
Rscript "${SCRIPT_DIR}/figure2_celltype_validation_focus.R"
echo "Finished at $(date)"
