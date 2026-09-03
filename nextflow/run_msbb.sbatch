#!/bin/bash
#SBATCH --job-name=msbb
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/msbb_%j.out
#SBATCH --error=logs/msbb_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd ${NF_DIR}/run_dirs/MSBB

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.msbb \
  -work-dir ${NF_DIR}/work/MSBB
