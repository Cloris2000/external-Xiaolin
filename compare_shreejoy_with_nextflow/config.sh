#!/bin/bash
# Shared paths and locus windows for the pipeline comparison.
# Everything here is read-only against both result trees.

# Output root (this directory)
OUT_DIR="/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow"
DATA_DIR="${OUT_DIR}/data"
FIG_DIR="${OUT_DIR}/figures"

# Xiaolin's Nextflow pipeline (hg19, 15-cohort METAL meta-analysis)
MINE_ROOT="/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
MINE_META="${MINE_ROOT}/results/meta_analysis_15cohorts_hg19_v2"
MINE_HET="${MINE_META}/heterogeneity/meta_heterogeneity_summary.tsv"
# Stage 1: MGP bulk-vs-snRNA agreement per broad class
MINE_DECONV="${MINE_ROOT}/manuscript_figure/figure2_celltype_accuracy.tsv"
# Stage 4: coloc, one file per cell type
MINE_COLOC="/scratch/zhoux156/results/downstream_v2/coloc/coloc_results_full"

# Shreejoy's phase-2 CTP-GWAS (hg38, pooled mega-analysis)
HIS_ROOT="/scratch/shreejoy/ctpgwas/results"
HIS_STEP2="${HIS_ROOT}/step2"
HIS_LOCI="${HIS_ROOT}/loci/annotated_hits.csv"
HIS_TIERS="${HIS_ROOT}/tiers/tiers.csv"
HIS_LAMBDA="${HIS_ROOT}/lambda/lambda_by_trait.csv"
# Stage 1: CelMod bulk-vs-snRNA agreement per supertype
HIS_ARM_AGREE="${HIS_ROOT}/person_pheno/arm_agreement.csv"
# Stage 3: his own inverse-variance meta across ancestry strata, same people
HIS_META_HITS="${HIS_ROOT}/meta_ancestry/hits.csv"
HIS_META_LOCI="${HIS_ROOT}/meta_ancestry/loci.csv"
# Stage 4
HIS_COLOC="${HIS_ROOT}/coloc/coloc_results.tsv"
HIS_REPO="/project/rrg-shreejoy/zhoux156/shreejoy_pipeline"
HIS_TAXONOMY="${HIS_REPO}/celltype-composition/refs/taxonomy_DFC_2026.tsv"

# R interpreter used elsewhere in this project
RSCRIPT="/home/zhoux156/miniforge3/envs/test/bin/Rscript"

# Locus windows.
#
# TMEM106B  hg19 chr7:12,250,867-12,282,993   hg38 chr7:12,211,240-12,243,367
# GRN       hg19 chr17:42,422,491-42,430,470  hg38 chr17:44,345,086-44,353,106
#
# Windows are +/-350 kb around each gene, given separately per build because no
# liftover chain file is available on this system. 04_align_coords.py derives the
# hg19->hg38 offset empirically from variants that match on REF:ALT.
MY_WIN="7:11900000:12600000 17:42100000:42800000"
HIS_WIN="7:11860000:12560000 17:44020000:44720000"

mkdir -p "${DATA_DIR}" "${FIG_DIR}"
