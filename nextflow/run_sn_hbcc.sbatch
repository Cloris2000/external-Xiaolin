#!/bin/bash
#SBATCH --job-name=sn_hbcc
#SBATCH --partition=mediumtmp
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/sn_hbcc_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/sn_hbcc_%j.err

set -euo pipefail

export JAVA_HOME=/nethome/kcni/xzhou/.anaconda3/envs/nextflow
export PATH="$JAVA_HOME/bin:$PATH"

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

mkdir -p logs/sn_NIMH_HBCC

nextflow run combined_pipeline_v2.nf \
    -c nextflow.config.combined.sn_hbcc \
    -w work/sn_NIMH_HBCC \
    -resume
