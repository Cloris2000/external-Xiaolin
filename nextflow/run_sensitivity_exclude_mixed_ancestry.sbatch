#!/bin/bash
#SBATCH --job-name=meta_sens_exclude_mixed_ancestry
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/meta_sensitivity/exclude_mixed_ancestry_%j.out
#SBATCH --error=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/meta_sensitivity/exclude_mixed_ancestry_%j.err
#SBATCH --account=rrg-shreejoy

cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow
nextflow run /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/standalone_meta_analysis.nf \
    -c /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/sensitivity_configs/nextflow.config.sensitivity.exclude_mixed_ancestry \
    -work-dir /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/work/meta_sensitivity/exclude_mixed_ancestry \
    -resume
