#!/bin/bash
#SBATCH --job-name=hbcc_1M
#SBATCH --partition=medium
#SBATCH --time=8:00:00       # 8 hours
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/hbcc_1M_%j.out
#SBATCH --error=logs/hbcc_1M_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/NIMH_HBCC_1M

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.nimh_hbcc_1M \
  -work-dir ${NF_DIR}/work/NIMH_HBCC_1M \
  -resume
