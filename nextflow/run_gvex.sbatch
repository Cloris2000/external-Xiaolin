#!/bin/bash
#SBATCH --job-name=gvex
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/gvex_%j.out
#SBATCH --error=logs/gvex_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd ${NF_DIR}/run_dirs/GVEX

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.gvex \
  -work-dir ${NF_DIR}/work/GVEX
