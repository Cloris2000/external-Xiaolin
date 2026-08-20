#!/bin/bash
#SBATCH --job-name=figure2_validation
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/fig2val_%j.log
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/fig2val_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail
SCRIPT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure"
echo "Starting at $(date)"
Rscript "${SCRIPT_DIR}/figure2_celltype_phenotypes_validation.R"
echo "Finished at $(date)"
