#!/bin/bash
#SBATCH --job-name=amp_ad_mayo
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/AMP_AD_Mayo/amp_ad_mayo_%j.out
#SBATCH --error=logs/AMP_AD_Mayo/amp_ad_mayo_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/AMP_AD_Mayo

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.amp_ad_mayo \
  -work-dir ${NF_DIR}/work/AMP_AD_Mayo \
  -resume d1057b0a-fa42-492b-a8c9-96aea00afedd
