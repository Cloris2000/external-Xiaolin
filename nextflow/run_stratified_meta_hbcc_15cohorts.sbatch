#!/bin/bash
#SBATCH --job-name=meta_strat_HBCC
#SBATCH --partition=medium
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_sensitivity/stratified_hbcc_%j.out
#SBATCH --error=logs/meta_sensitivity/stratified_hbcc_%j.err

# AFR-enriched stratum: NIMH HBCC 3 platforms only.
# Run from a dedicated run_dir to avoid Nextflow history lock conflicts.

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
mkdir -p "${NF_DIR}/run_dirs/stratified_hbcc_15cohorts"
cd "${NF_DIR}/run_dirs/stratified_hbcc_15cohorts"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/sensitivity_configs/nextflow.config.sensitivity.stratified_hbcc_15cohorts" \
    -work-dir "${NF_DIR}/work/meta_sensitivity/stratified_hbcc_15cohorts" \
    -resume
