#!/bin/bash
#SBATCH --job-name=meta_anc_EUR
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/meta_sensitivity/ancestry_EUR_%j.out
#SBATCH --error=logs/meta_sensitivity/ancestry_EUR_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
mkdir -p "${NF_DIR}/logs/meta_sensitivity" "${NF_DIR}/run_dirs/ancestry_meta_EUR"
cd "${NF_DIR}/run_dirs/ancestry_meta_EUR"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/sensitivity_configs/nextflow.config.sensitivity.ancestry_EUR_meta" \
    -work-dir "${NF_DIR}/work/meta_sensitivity/ancestry_specific_EUR" \
    -resume
