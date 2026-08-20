#!/usr/bin/env bash
# Transfer NABEC cohort inputs CAMH -> Trillium.
#
# Prerequisites on CAMH:
#   bash stage_remaining_cohorts_on_camh.sh --cohort NABEC
#
# Run on Trillium (tmux):
#   bash transfer_nabec.sh
#
# Dest: /project/rrg-shreejoy/NABEC
# Size: ~6.5G genotypes + RNA/meta

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/NABEC}"
STAGE_NABEC="${STAGE_ROOT}/NABEC"
VCF_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/NABEC_WGS/normalized_vcf"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p "${DEST_ROOT}/RNA" "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Genotype/normalized_vcf"

trillium_rsync "nabec_rna" \
  "${STAGE_NABEC}/RNA/" \
  "${DEST_ROOT}/RNA/"

# May be empty if biospecimen mapping was missing at staging time
trillium_rsync "nabec_meta" \
  "${STAGE_NABEC}/Metadata/" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "nabec_vcf" \
  "${VCF_DIR}/" \
  "${DEST_ROOT}/Genotype/normalized_vcf/"

echo "NABEC transfer finished -> ${DEST_ROOT}"
echo "  RNA/ (counts, metadata_combined, combined_metrics)"
echo "  Metadata/ (biospecimen mapping if staged)"
echo "  Genotype/normalized_vcf/   # ~6.5G"
echo "NOTE: If Metadata/NABEC_biospecimen_mapping.txt is missing, locate it on CAMH before running NABEC GWAS."
