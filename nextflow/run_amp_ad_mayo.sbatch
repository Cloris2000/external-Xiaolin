#!/bin/bash
#SBATCH --job-name=amp_ad_mayo
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/zhoux156/logs/AMP_AD_Mayo/amp_ad_mayo_%j.out
#SBATCH --error=/scratch/zhoux156/logs/AMP_AD_Mayo/amp_ad_mayo_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
module load "${NEXTFLOW_MODULE}"
cd ${NF_DIR}/run_dirs/AMP_AD_Mayo

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.amp_ad_mayo \
  -work-dir ${WORK_ROOT}/AMP_AD_Mayo
