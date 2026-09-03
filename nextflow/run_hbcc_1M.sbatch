#!/bin/bash
#SBATCH --job-name=hbcc_1M
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/hbcc_1M_%j.out
#SBATCH --error=logs/hbcc_1M_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd ${NF_DIR}/run_dirs/NIMH_HBCC_1M

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.nimh_hbcc_1M \
  -work-dir ${NF_DIR}/work/NIMH_HBCC_1M
