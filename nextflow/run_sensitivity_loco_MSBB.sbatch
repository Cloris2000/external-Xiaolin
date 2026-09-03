#!/bin/bash
#SBATCH --job-name=meta_sens_loco_MSBB
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/meta_sensitivity/loco_MSBB_%j.out
#SBATCH --error=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/meta_sensitivity/loco_MSBB_%j.err
#SBATCH --account=rrg-shreejoy

cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow
nextflow run /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/standalone_meta_analysis.nf \
    -c /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/sensitivity_configs/nextflow.config.sensitivity.loco_MSBB \
    -work-dir /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/work/meta_sensitivity/loco_MSBB \
    -resume
