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
MSSM_VCF="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/CMC_MSSM_reimputed/imputed_chr_normalized"
PENN_PITT_VCF="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/CMC_imputation/imputed_chr_normalized"

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
  "${DEST_ROOT}/Genotype/CMC_MSSM_imputed_chr_normalized" \
  "${DEST_ROOT}/Genotype/CMC_PENN_PITT_imputed_chr_normalized"

trillium_rsync "cmc_rna" \
  "${STAGE_CMC}/RNA/" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "cmc_snp_meta" \
  "${STAGE_CMC}/Metadata/" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "cmc_mssm_vcf" \
  "${MSSM_VCF}/" \
  "${DEST_ROOT}/Genotype/CMC_MSSM_imputed_chr_normalized/"

trillium_rsync "cmc_penn_pitt_vcf" \
  "${PENN_PITT_VCF}/" \
  "${DEST_ROOT}/Genotype/CMC_PENN_PITT_imputed_chr_normalized/"

echo "CMC transfer finished -> ${DEST_ROOT}"
echo "  RNA/CMC_{MSSM,PENN,PITT}_{count_matrix,metadata}.csv"
echo "  Metadata/CMC_Human_SNP_metadata.csv"
echo "  Genotype/CMC_MSSM_imputed_chr_normalized/           # ~12G"
echo "  Genotype/CMC_PENN_PITT_imputed_chr_normalized/      # ~30G (shared)"
