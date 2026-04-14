#!/bin/bash
#SBATCH --job-name=gtex
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/gtex_%j.out
#SBATCH --error=logs/gtex_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/GTEx

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.gtex \
  -work-dir ${NF_DIR}/work/GTEx \
  -resume
