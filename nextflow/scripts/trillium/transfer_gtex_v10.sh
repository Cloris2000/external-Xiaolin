#!/usr/bin/env bash
# Transfer GTEx_v10 cohort inputs CAMH -> Trillium.
# All on netdata_kcni (no staging). Largest single genotype transfer (~766G).
#
# Run on Trillium (tmux):
#   bash transfer_gtex_v10.sh
#
# Dest: /project/rrg-shreejoy/GTEx_v10

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/GTEx_v10}"
RNA_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/data"
VCF_DIR="/external/rprshnas01/netdata_kcni/stlab/GTEx_v10/Genotype/split_by_chr"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p "${DEST_ROOT}/RNA" "${DEST_ROOT}/Metadata" "${DEST_ROOT}/Genotype/split_by_chr"

trillium_rsync "gtex_v10_counts" \
  "${RNA_DIR}/GTEx_v10_BA9_raw_count_matrix.csv" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "gtex_v10_meta" \
  "${RNA_DIR}/GTEx_v10_BA9_sample_metadata.csv" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "gtex_v10_vcf" \
  "${VCF_DIR}/" \
  "${DEST_ROOT}/Genotype/split_by_chr/"

echo "GTEx_v10 transfer finished -> ${DEST_ROOT}"
echo "  RNA/GTEx_v10_BA9_raw_count_matrix.csv"
echo "  Metadata/GTEx_v10_BA9_sample_metadata.csv"
echo "  Genotype/split_by_chr/   # ~766G WGS"
