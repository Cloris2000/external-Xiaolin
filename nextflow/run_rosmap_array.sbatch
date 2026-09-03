#!/bin/bash
#SBATCH --job-name=rosmap_array
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/rosmap_array_%j.out
#SBATCH --error=logs/rosmap_array_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd ${NF_DIR}/run_dirs/ROSMAP_array

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.rosmap_array \
  -work-dir ${NF_DIR}/work/ROSMAP_array
