#!/bin/bash
#SBATCH --job-name=sn_hbcc
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/zhoux156/logs/sn_hbcc_%j.out
#SBATCH --error=/scratch/zhoux156/logs/sn_hbcc_%j.err
#SBATCH --account=rrg-shreejoy

set -euo pipefail

export PATH="$JAVA_HOME/bin:$PATH"

cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow

mkdir -p logs/sn_NIMH_HBCC

nextflow run combined_pipeline_v2.nf \
    -c nextflow.config.combined.sn_hbcc \
    -w work/sn_NIMH_HBCC \
    -resume
