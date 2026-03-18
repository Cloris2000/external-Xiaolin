#!/bin/bash
#SBATCH --job-name=hbcc_h650
#SBATCH --partition=medium
#SBATCH --time=8:00:00       # 8 hours
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/hbcc_h650_%j.out
#SBATCH --error=logs/hbcc_h650_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/NIMH_HBCC_h650

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.nimh_hbcc_h650 \
  -work-dir ${NF_DIR}/work/NIMH_HBCC_h650 \
  -resume