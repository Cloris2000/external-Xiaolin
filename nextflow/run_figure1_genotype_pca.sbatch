#!/bin/bash
#SBATCH --job-name=fig1_geno_pca
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/zhoux156/logs/figure1_genotype_pca/pca_%j.out
#SBATCH --error=/scratch/zhoux156/logs/figure1_genotype_pca/pca_%j.err
#
# Figure 1 genotype panels (A: genetic ancestry, B: cohort, D: composition)
# for the 15-cohort meta-analysis.
#
# Runs 1000 Genomes reference-projection PCA over the QC'd genotypes of all
# 15 cohorts, then draws the three panels as editable SVG/PDF/PNG.
#
# Read-only with respect to every pipeline input: it reads the per-cohort
# *.QC.final.{pgen,pvar,psam} and the 1000G reference panel, and writes only
# under ${RESULTS_ROOT}/figure1_genotype_pca.
#
# Runtime is dominated by the LD-pruning pass over the 3,202-sample 1000G pgen.
#
# Submit:
#   sbatch run_figure1_genotype_pca.sbatch
#
set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

mkdir -p "${LOG_ROOT}/figure1_genotype_pca"
cd "${NF_DIR}"

export PLINK_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export PLINK_MEM=600000          # MB, of the node's 767 GB
export OUT_DIR="${RESULTS_ROOT}/figure1_genotype_pca"

echo "=== Figure 1 genotype PCA ==="
echo "host      : $(hostname)"
echo "threads   : ${PLINK_THREADS}"
echo "output    : ${OUT_DIR}"
echo

bash "${NF_DIR}/scripts/figure1_genotype_pca.sh"

echo
echo "=== Drawing panels ==="
"${HOME}/miniforge3/envs/python_env/bin/python" \
    "${NF_DIR}/scripts/figure1_genotype_panels.py" \
    --samples "${OUT_DIR}/figure1_pca_samples.tsv" \
    --out-dir "${OUT_DIR}/figures"

echo
echo "Done.  Panels (SVG for Inkscape, plus PDF/PNG):"
echo "  ${OUT_DIR}/figures/figure1_panelA_ancestry.svg"
echo "  ${OUT_DIR}/figures/figure1_panelB_cohort.svg"
echo "  ${OUT_DIR}/figures/figure1_panelD_composition.svg"
echo
echo "Tables:"
echo "  ${OUT_DIR}/figure1_pca_samples.tsv   per-sample PCs + ancestry"
echo "  ${OUT_DIR}/figure1_pca_summary.tsv   per-cohort composition"
echo
echo "Scratch is purged periodically - sync anything you want to keep:"
echo "  bash scripts/sync_results_to_project.sh"
