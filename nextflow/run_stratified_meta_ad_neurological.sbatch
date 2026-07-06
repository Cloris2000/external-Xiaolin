#!/bin/bash
#SBATCH --job-name=meta_strat_AD
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_sensitivity/stratified_ad_neurological_%j.out
#SBATCH --error=logs/meta_sensitivity/stratified_ad_neurological_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
mkdir -p "${NF_DIR}/work/meta_sensitivity/stratified_ad_neurological_15cohorts"
cd "${NF_DIR}"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/sensitivity_configs/nextflow.config.sensitivity.stratified_ad_neurological_15cohorts" \
    -work-dir "${NF_DIR}/work/meta_sensitivity/stratified_ad_neurological_15cohorts" \
    -resume
