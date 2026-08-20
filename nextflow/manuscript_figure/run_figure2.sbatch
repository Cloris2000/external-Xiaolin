#!/bin/bash
#SBATCH --job-name=figure2_celltype_props
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/figure2_%j.log
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/figure2_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail

SCRIPT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure"

echo "Starting Figure 2 script at $(date)"
Rscript "${SCRIPT_DIR}/figure2_celltype_proportion_phenotypes.R"
echo "Finished at $(date)"
