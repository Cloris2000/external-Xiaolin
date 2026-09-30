#!/bin/bash
#SBATCH --job-name=figure2_focus_rint_adj
#SBATCH --output=/scratch/zhoux156/logs/manuscript_figure/fig2focus_rint_adj_%j.log
#SBATCH --error=/scratch/zhoux156/logs/manuscript_figure/fig2focus_rint_adj_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute

# Figure 2 validation, SENSITIVITY MODE: compares the two sides on the scale
# actually used for the GWAS.
#   Bulk  - the stored GWAS phenotypes, results/<cohort>/phenotypes_RINT.txt,
#           matched per donor via each cohort's biospecimen/ID bridge. These
#           are NOT recomputed; they are what REGENIE saw.
#   snRNA - has no stored equivalent, so it is residualised on the same sex/age
#           design (pheno_cov_prep.R) and RINT-transformed within each
#           validation cohort x cell type.
# Writes *_rint_adj outputs only; the default Figure 2 files are not touched.
#
# Validated by a login-node dry run on 2026-09-17: 22030 of 22258 paired rows
# retained (losses are ROSMAP/Mathys donors with no WGS specimen in the stored
# phenotype file). Runtime was a few minutes, so this is light enough to have
# been checked interactively, but is submitted here to keep compute off the
# login node.

set -euo pipefail
RSCRIPT="/home/zhoux156/miniforge3/envs/test/bin/Rscript"
SCRIPT_DIR="/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"
SCRATCH_OUT="/scratch/zhoux156/manuscript_figure"
FINAL_OUT="${SCRIPT_DIR}"

mkdir -p "${SCRATCH_OUT}"

export FIG2_OUT_DIR="${SCRATCH_OUT}"
export FIG2_SENSITIVITY="rint_adj"
echo "Starting at $(date)"
echo "Rscript: ${RSCRIPT}"
echo "FIG2_OUT_DIR=${FIG2_OUT_DIR}"
echo "FIG2_SENSITIVITY=${FIG2_SENSITIVITY}"
"${RSCRIPT}" "${SCRIPT_DIR}/figure2_celltype_validation_focus.R"

echo "Copying outputs to ${FINAL_OUT}"
cp -f "${SCRATCH_OUT}/combined_bulk_snrna_paired_rint_adj.tsv" \
      "${SCRATCH_OUT}/figure2_celltype_accuracy_rint_adj.tsv" \
      "${SCRATCH_OUT}/figure2_celltype_validation_focus_rint_adj."{png,pdf,svg} \
      "${FINAL_OUT}/"
echo "Finished at $(date)"
