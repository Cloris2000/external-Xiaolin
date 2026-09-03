#!/bin/bash
#SBATCH --job-name=anc_gwas
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-11
#SBATCH --output=/scratch/zhoux156/logs/ancestry_gwas_%A_%a.out
#SBATCH --error=/scratch/zhoux156/logs/ancestry_gwas_%A_%a.err
#SBATCH --account=rrg-shreejoy

# Ancestry-subset GWAS array.
# Tasks come from docs/ancestry_specific/ancestry_gwas_array_tasks.tsv (re_gwas only).
# Array indices 0..N-1 must match that TSV row order.
#
# Pilot single task:
#   sbatch --array=0 run_ancestry_gwas_array.sbatch
# Full:
#   sbatch run_ancestry_gwas_array.sbatch

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd "${NF_DIR}"
mkdir -p logs

TASKS="${NF_DIR}/docs/ancestry_specific/ancestry_gwas_array_tasks.tsv"
# Skip header; SLURM_ARRAY_TASK_ID is 0-based
LINE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i {print}' "${TASKS}")
if [[ -z "${LINE}" ]]; then
  echo "No task for array index ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

STUDY=$(echo "${LINE}" | cut -f1)
CONFIG=$(echo "${LINE}" | cut -f6)

echo "=== Ancestry GWAS task ${SLURM_ARRAY_TASK_ID}: ${STUDY} ==="
echo "Config: ${CONFIG}"
date

mkdir -p "run_dirs/${STUDY}" "work/${STUDY}" "logs/${STUDY}" "results/${STUDY}"
cd "run_dirs/${STUDY}"

nextflow run "${NF_DIR}/combined_pipeline_v2.nf" \
  -c "${NF_DIR}/${CONFIG}" \
  -work-dir "${NF_DIR}/work/${STUDY}" \
  -resume

echo "=== Done ${STUDY} ==="
date
