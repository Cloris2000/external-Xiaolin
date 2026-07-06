#!/bin/bash
#SBATCH --job-name=meta_strat_ctrl
#SBATCH --partition=medium
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_sensitivity/stratified_neurologically_normal_%j.out
#SBATCH --error=logs/meta_sensitivity/stratified_neurologically_normal_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/sensitivity_configs/nextflow.config.sensitivity.stratified_neurologically_normal_15cohorts" \
    -work-dir "${NF_DIR}/work/meta_sensitivity/stratified_neurologically_normal_15cohorts" \
    -resume
