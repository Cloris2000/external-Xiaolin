#!/bin/bash
#SBATCH --job-name=meta_strat_EUR
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_sensitivity/stratified_eur_homogeneous_%j.out
#SBATCH --error=logs/meta_sensitivity/stratified_eur_homogeneous_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/sensitivity_configs/nextflow.config.sensitivity.stratified_eur_homogeneous_15cohorts" \
    -work-dir "${NF_DIR}/work/meta_sensitivity/stratified_eur_homogeneous_15cohorts" \
    -resume
