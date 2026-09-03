#!/bin/bash
#SBATCH --job-name=meta_sens_hbcc_platforms_only
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/zhoux156/logs/meta_sensitivity/hbcc_platforms_only_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_sensitivity/hbcc_platforms_only_%j.err
#SBATCH --account=rrg-shreejoy

cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow
nextflow run /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/standalone_meta_analysis.nf \
    -c /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/sensitivity_configs/nextflow.config.sensitivity.hbcc_platforms_only \
    -work-dir /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/work/meta_sensitivity/hbcc_platforms_only \
    -resume
