#!/bin/bash
#SBATCH --job-name=figure3_gwas_discovery
#SBATCH --output=/scratch/zhoux156/logs/manuscript_figure/figure3_%j.log
#SBATCH --error=/scratch/zhoux156/logs/manuscript_figure/figure3_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --account=rrg-shreejoy

set -euo pipefail

SCRIPT_DIR="/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"

echo "Starting Figure 3 script at $(date)"
Rscript "${SCRIPT_DIR}/figure3_gwas_discovery.R"
echo "Finished at $(date)"
