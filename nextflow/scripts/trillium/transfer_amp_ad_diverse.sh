#!/usr/bin/env bash
# Transfer AMP-AD Diverse (Rush + Mayo) RNA/metadata + normalized genotypes.
#
# SKIP DivCo raw VCFs — already on Trillium at:
#   /project/rrg-shreejoy/AMP_AD_Diverse/Genotype/
#
# Run on Trillium (tmux):
#   bash transfer_amp_ad_diverse.sh
#
# Dest: /project/rrg-shreejoy/AMP_AD_Diverse
# Size: ~268G normalized genotypes + small RNA/meta

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/AMP_AD_Diverse}"
AMP_BASE="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/AMP_AD_Diverse"
DATA_INPUT="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/nextflow/data_input/amp_ad_diverse"
SKIP_NORMALIZED="${SKIP_NORMALIZED:-0}"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]

  SKIP_NORMALIZED=1  skip Genotype/normalized (~268G) if you will re-normalize
                     from existing DivCo raw VCFs on Trillium
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p \
  "${DEST_ROOT}/RNA" \
  "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Bulk_RNA_Seq_and_QC/Mayo_emory" \
  "${DEST_ROOT}/Genotype/normalized"

# Staged pipeline RNA/meta
trillium_rsync "amp_ad_rna_meta" \
  "${DATA_INPUT}/" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "amp_ad_biospec" \
  "${AMP_BASE}/Metadata/" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "amp_ad_mayo_qc" \
  "${AMP_BASE}/Bulk_RNA_Seq_and_QC/Mayo_emory/Combined_QC_Metrics.csv" \
  "${DEST_ROOT}/Bulk_RNA_Seq_and_QC/Mayo_emory/"

if [[ "${SKIP_NORMALIZED}" -eq 1 ]]; then
  echo "SKIP_NORMALIZED=1 — not transferring Genotype/normalized"
  echo "  Using existing DivCo raw under ${DEST_ROOT}/Genotype/ (re-normalize on Trillium if needed)"
else
  trillium_rsync "amp_ad_vcf_norm" \
    "${AMP_BASE}/Genotype/normalized/" \
    "${DEST_ROOT}/Genotype/normalized/"
fi

echo "AMP_AD_Diverse transfer finished -> ${DEST_ROOT}"
echo "  RNA/ (Rush + Mayo counts/metadata from data_input)"
echo "  Metadata/AMP-AD_DiverseCohorts_biospecimen_metadata.csv"
echo "  Bulk_RNA_Seq_and_QC/Mayo_emory/Combined_QC_Metrics.csv"
echo "  Genotype/normalized/   # pipeline-ready (unless skipped)"
echo "  Genotype/DivCo*.vcf.gz # already present — do not re-transfer"
