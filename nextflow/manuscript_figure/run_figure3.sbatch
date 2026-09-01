#!/bin/bash
#SBATCH --job-name=figure3_gwas_discovery
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/figure3_%j.log
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure/figure3_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2

set -euo pipefail

SCRIPT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure"

echo "Starting Figure 3 script at $(date)"
Rscript "${SCRIPT_DIR}/figure3_gwas_discovery.R"
echo "Finished at $(date)"
