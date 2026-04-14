#!/bin/bash
#SBATCH --job-name=cmc_mssm
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/cmc_mssm_%j.out
#SBATCH --error=logs/cmc_mssm_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/CMC_MSSM

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.cmc_mssm \
  -work-dir ${NF_DIR}/work/CMC_MSSM
