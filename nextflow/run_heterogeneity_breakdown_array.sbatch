#!/bin/bash
#SBATCH --job-name=het_ct
#SBATCH --partition=compute
#SBATCH --array=0-18
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/zhoux156/logs/het_breakdown/het_breakdown_array_%A_%a.out
#SBATCH --error=/scratch/zhoux156/logs/het_breakdown/het_breakdown_array_%A_%a.err
#SBATCH --account=rrg-shreejoy

# One cell type per array task — heterogeneity breakdown worker.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
mkdir -p "${NF_DIR}/logs/het_breakdown"
cd "${NF_DIR}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

COHORTS=(
  AMP_AD_Mayo AMP_AD_Rush
  CMC_MSSM CMC_PENN CMC_PITT
  GTEx_v10 GVEX MSBB Mayo NABEC
  NIMH_HBCC_1M NIMH_HBCC_Omni5M NIMH_HBCC_h650
  ROSMAP ROSMAP_array
)

CELL_TYPE="${CELL_TYPES[$SLURM_ARRAY_TASK_ID]}"
OUT_DIR="${NF_DIR}/results/meta_analysis_15cohorts/heterogeneity_breakdown"

echo "Processing cell type: ${CELL_TYPE} (array task ${SLURM_ARRAY_TASK_ID})"

python3 "${NF_DIR}/scripts/heterogeneity_breakdown.py" \
    --cell-type "${CELL_TYPE}" \
    --meta-tbl-dir "${NF_DIR}/results/meta_analysis_15cohorts" \
    --cohort-sumstat-pattern "${NF_DIR}/results/{cohort}/regenie_step2/{cohort}_{cell_type}_step2.regenie.raw_p" \
    --cohort-metadata "${NF_DIR}/docs/meta_cohort_metadata.tsv" \
    --output-dir "${OUT_DIR}" \
    --cohorts "${COHORTS[@]}" \
    --sumstat-results-dir "${NF_DIR}/results"

echo "Done: ${CELL_TYPE}"
