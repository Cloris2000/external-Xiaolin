#!/usr/bin/env bash
# Transfer CMC_MSSM / CMC_PENN / CMC_PITT inputs CAMH -> Trillium.
#
# Prerequisites on CAMH:
#   bash stage_remaining_cohorts_on_camh.sh --cohort CMC
#
# Run on Trillium (tmux):
#   bash transfer_cmc.sh
#
# Dest: /project/rrg-shreejoy/CMC
# Size: ~42G genotypes + RNA/meta
# Note: PENN and PITT share the same imputed VCF directory (copied once).

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/CMC}"
STAGE_CMC="${STAGE_ROOT}/CMC"
# Original QC'd genotyping array (GRCh37 plink binary) — the TRUE raw input for
# all three CMC cohorts (MSSM+Penn+Pitt were genotyped/imputed as one batch).
# Synapse syn4600985/87/89. On netdata -> direct pull. NOTE: no remote wildcard
# (the CAMH rssh mover does not expand globs) — transfer each file explicitly.
CMC_RAW_ARRAY_BASE="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/CMC_genotypes/SNPs/Release1/QCd/CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_QCed"
MSSM_VCF="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/WGS/QC/CMC_MSSM_reimputed/imputed_chr_normalized"
PENN_PITT_VCF="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/CMC_imputation/imputed_chr_normalized"
# Raw (pre-normalization) imputed dose VCFs for PENN/PITT (~60G):
#   chr{N}.dose.vcf.gz + chr{N}.info.gz + chr_{N}.zip (raw imputation-server output)
PENN_PITT_VCF_RAW="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/CMC_imputation/imputed_chr"
# Per-cohort normalized dirs referenced by config `normalized_vcf_dir` (the active
# genotype input is vcf_pattern -> CMC_imputation, but copy these too for safety).
PENN_VCF_REIMP="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/WGS/QC/CMC_PENN_reimputed/imputed_chr_normalized"
PITT_VCF_REIMP="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/WGS/QC/CMC_PITT_reimputed/imputed_chr_normalized"
# CMC MSSM has no raw imputed-VCF dir on CAMH (only imputed_chr_normalized/ + QC
# pgen). Set MSSM_VCF_RAW to the correct path once located to enable its copy.
MSSM_VCF_RAW="${MSSM_VCF_RAW:-}"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p \
  "${DEST_ROOT}/RNA" \
  "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Genotype/CMC_raw_genotyping_array" \
  "${DEST_ROOT}/Genotype/CMC_MSSM_imputed_chr_normalized" \
  "${DEST_ROOT}/Genotype/CMC_PENN_PITT_imputed_chr_normalized"

trillium_rsync "cmc_rna" \
  "${STAGE_CMC}/RNA/" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "cmc_snp_meta" \
  "${STAGE_CMC}/Metadata/" \
  "${DEST_ROOT}/Metadata/"

# Original raw QC'd genotyping array (.bed/.bim/.fam + README) — shared MSSM+Penn+Pitt
# One explicit rsync per file (no remote glob, which rssh won't expand).
for ext in .bed .bim .fam _README.txt; do
  trillium_rsync "cmc_raw_array${ext//./_}" \
    "${CMC_RAW_ARRAY_BASE}${ext}" \
    "${DEST_ROOT}/Genotype/CMC_raw_genotyping_array/"
done

trillium_rsync "cmc_mssm_vcf" \
  "${MSSM_VCF}/" \
  "${DEST_ROOT}/Genotype/CMC_MSSM_imputed_chr_normalized/"

trillium_rsync "cmc_penn_pitt_vcf" \
  "${PENN_PITT_VCF}/" \
  "${DEST_ROOT}/Genotype/CMC_PENN_PITT_imputed_chr_normalized/"

# Raw (pre-normalization) imputed VCFs
trillium_rsync "cmc_penn_pitt_vcf_raw" \
  "${PENN_PITT_VCF_RAW}/" \
  "${DEST_ROOT}/Genotype/CMC_PENN_PITT_imputed_chr_raw/"

# Per-cohort reimputed normalized dirs (config normalized_vcf_dir references)
trillium_rsync "cmc_penn_reimputed_vcf" \
  "${PENN_VCF_REIMP}/" \
  "${DEST_ROOT}/Genotype/CMC_PENN_reimputed_imputed_chr_normalized/"

trillium_rsync "cmc_pitt_reimputed_vcf" \
  "${PITT_VCF_REIMP}/" \
  "${DEST_ROOT}/Genotype/CMC_PITT_reimputed_imputed_chr_normalized/"

if [[ -n "${MSSM_VCF_RAW}" ]]; then
  trillium_rsync "cmc_mssm_vcf_raw" \
    "${MSSM_VCF_RAW}/" \
    "${DEST_ROOT}/Genotype/CMC_MSSM_imputed_chr_raw/"
else
  echo "MSSM_VCF_RAW not set — skipping CMC MSSM raw (no raw dir located yet)"
fi

echo "CMC transfer finished -> ${DEST_ROOT}"
echo "  RNA/CMC_{MSSM,PENN,PITT}_{count_matrix,metadata}.csv"
echo "  Metadata/CMC_Human_SNP_metadata.csv"
echo "  Genotype/CMC_raw_genotyping_array/                     # ~140M original QC'd array (bed/bim/fam+README)"
echo "  Genotype/CMC_MSSM_imputed_chr_normalized/              # ~12G (active input)"
echo "  Genotype/CMC_PENN_PITT_imputed_chr_normalized/         # ~30G (active input, shared)"
echo "  Genotype/CMC_PENN_PITT_imputed_chr_raw/                # ~60G raw (shared)"
echo "  Genotype/CMC_PENN_reimputed_imputed_chr_normalized/    # ~3.8G (config ref)"
echo "  Genotype/CMC_PITT_reimputed_imputed_chr_normalized/    # ~6.5G (config ref)"
if [[ -n "${MSSM_VCF_RAW}" ]]; then
  echo "  Genotype/CMC_MSSM_imputed_chr_raw/                     # MSSM raw"
fi
exit 0
