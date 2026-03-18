#!/bin/bash
#SBATCH --job-name=gvex
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/gvex_%j.out
#SBATCH --error=logs/gvex_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/GVEX

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.gvex \
  -work-dir ${NF_DIR}/work/GVEX \
  -resume
