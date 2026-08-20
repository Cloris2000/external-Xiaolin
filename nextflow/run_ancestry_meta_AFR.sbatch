#!/bin/bash
#SBATCH --job-name=meta_anc_AFR
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/meta_sensitivity/ancestry_AFR_%j.out
#SBATCH --error=logs/meta_sensitivity/ancestry_AFR_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
# Launch from project root so -resume can reuse the original AFR session cache
cd "${NF_DIR}"

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/sensitivity_configs/nextflow.config.sensitivity.ancestry_AFR_meta" \
    -work-dir "${NF_DIR}/work/meta_sensitivity/ancestry_specific_AFR" \
    -resume
