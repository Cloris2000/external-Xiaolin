#!/bin/bash
#SBATCH --job-name=rosmap_array
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/rosmap_array_%j.out
#SBATCH --error=logs/rosmap_array_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/ROSMAP_array

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.rosmap_array \
  -work-dir ${NF_DIR}/work/ROSMAP_array
