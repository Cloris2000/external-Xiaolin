#!/bin/bash
#SBATCH --job-name=cmc_mssm
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/zhoux156/logs/cmc_mssm_%j.out
#SBATCH --error=/scratch/zhoux156/logs/cmc_mssm_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"
cd ${NF_DIR}/run_dirs/CMC_MSSM

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.cmc_mssm \
  -work-dir ${WORK_ROOT}/CMC_MSSM
