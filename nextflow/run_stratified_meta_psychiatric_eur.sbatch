#!/bin/bash
#SBATCH --job-name=meta_psych_EUR
#SBATCH --partition=medium
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_sensitivity/stratified_psychiatric_eur_%j.out
#SBATCH --error=logs/meta_sensitivity/stratified_psychiatric_eur_%j.err
# Legacy single-job launcher (superseded by run_stratified_meta_psychiatric_eur_array.sbatch).
# Use: bash run_diagnosis_eur_analysis.sbatch

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
mkdir -p "${NF_DIR}/run_dirs/stratified_psychiatric_eur_15cohorts"
cd "${NF_DIR}/run_dirs/stratified_psychiatric_eur_15cohorts"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/sensitivity_configs/nextflow.config.sensitivity.stratified_psychiatric_eur_15cohorts" \
    -work-dir "${NF_DIR}/work/meta_sensitivity/stratified_psychiatric_eur_15cohorts" \
    -resume
