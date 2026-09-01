#!/usr/bin/env bash
# Transfer Mayo cohort inputs CAMH -> Trillium.
# All Mayo inputs are already on netdata_kcni (no staging needed).
#
# Run on Trillium (tmux recommended):
#   bash transfer_mayo.sh
#   bash transfer_mayo.sh --dry-run
#
# Dest: /project/rrg-shreejoy/Mayo
# Size: ~217G genotypes + small RNA/meta

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/Mayo}"
META_PCA="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/metabrain_PCA/data"
WGS_BASE="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/WGS"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p "${DEST_ROOT}/RNA" "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Genotype/Mayo_joint_WGS_vcf" \
  "${DEST_ROOT}/Genotype/Mayo_joint_WGS_vcf_normalized"

trillium_rsync "mayo_counts" \
  "${META_PCA}/Mayo_raw_counts_Nov_18_ensembl.csv" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "mayo_meta" \
  "${META_PCA}/Mayo_meta_tissue_counts_Nov_18.csv" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "mayo_wgs_raw" \
  "${WGS_BASE}/Mayo_joint_WGS_vcf/" \
  "${DEST_ROOT}/Genotype/Mayo_joint_WGS_vcf/"

trillium_rsync "mayo_wgs_norm" \
  "${WGS_BASE}/Mayo_joint_WGS_vcf_normalized/" \
  "${DEST_ROOT}/Genotype/Mayo_joint_WGS_vcf_normalized/"

echo "Mayo transfer finished -> ${DEST_ROOT}"
echo "  RNA/Mayo_raw_counts_Nov_18_ensembl.csv"
echo "  Metadata/Mayo_meta_tissue_counts_Nov_18.csv"
echo "  Genotype/Mayo_joint_WGS_vcf/              # raw ~107G"
echo "  Genotype/Mayo_joint_WGS_vcf_normalized/   # normalized ~110G"
echo "Note: Mayo has no biospecimen file in the pipeline config."
