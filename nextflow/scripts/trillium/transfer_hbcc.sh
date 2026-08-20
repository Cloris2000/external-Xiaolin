#!/usr/bin/env bash
# Transfer NIMH_HBCC (3 platforms) inputs CAMH -> Trillium.
# RNA + imputed VCFs are on netdata. Skip re-copying SNP_array / BAMs already
# under /project/rrg-shreejoy/NIMH_HBCC/ if present.
#
# Run on Trillium (tmux):
#   bash transfer_hbcc.sh
#
# Dest: /project/rrg-shreejoy/NIMH_HBCC
# Size: ~13G imputed genotypes + RNA/meta

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/NIMH_HBCC}"
HBCC_RNA="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data_input/nimh_hbcc"
HBCC_VCF_BASE="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/CMC_HBCC_reimputed"
CMC_SNP_META="/external/rprshnas01/netdata_kcni/stlab/CMC_genotypes/SNPs/Release3/Metadata/CMC_Human_SNP_metadata.csv"

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
  "${DEST_ROOT}/Genotype/imputed_chr_normalized_1M" \
  "${DEST_ROOT}/Genotype/imputed_chr_normalized_h650" \
  "${DEST_ROOT}/Genotype/imputed_chr_normalized_Omni5M"

trillium_rsync "hbcc_rna" \
  "${HBCC_RNA}/" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "hbcc_snp_meta" \
  "${CMC_SNP_META}" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "hbcc_1M" \
  "${HBCC_VCF_BASE}/imputed_chr_normalized_1M/" \
  "${DEST_ROOT}/Genotype/imputed_chr_normalized_1M/"

trillium_rsync "hbcc_h650" \
  "${HBCC_VCF_BASE}/imputed_chr_normalized_h650/" \
  "${DEST_ROOT}/Genotype/imputed_chr_normalized_h650/"

trillium_rsync "hbcc_Omni5M" \
  "${HBCC_VCF_BASE}/imputed_chr_normalized_Omni5M/" \
  "${DEST_ROOT}/Genotype/imputed_chr_normalized_Omni5M/"

echo "NIMH_HBCC transfer finished -> ${DEST_ROOT}"
echo "  RNA/HBCC_count_matrix.csv, HBCC_metadata.csv, ..."
echo "  Metadata/CMC_Human_SNP_metadata.csv"
echo "  Genotype/imputed_chr_normalized_{1M,h650,Omni5M}/"
echo "Existing (not re-transferred): NIMH_HBCC_SNP_array, HBCC_bam_aligned if already present."
