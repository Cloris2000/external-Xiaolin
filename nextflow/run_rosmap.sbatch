#!/bin/bash
#SBATCH --job-name=rosmap
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/rosmap_%j.out
#SBATCH --error=logs/rosmap_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/ROSMAP

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.rosmap \
  -work-dir ${NF_DIR}/work/ROSMAP
