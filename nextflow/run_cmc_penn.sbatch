#!/bin/bash
#SBATCH --job-name=cmc_penn
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/zhoux156/logs/cmc_penn_%j.out
#SBATCH --error=/scratch/zhoux156/logs/cmc_penn_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd ${NF_DIR}/run_dirs/CMC_PENN

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.cmc_penn \
  -work-dir ${NF_DIR}/work/CMC_PENN
