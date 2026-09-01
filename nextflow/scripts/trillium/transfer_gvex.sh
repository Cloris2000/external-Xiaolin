#!/usr/bin/env bash
# Transfer GVEX (BrainGVEX) inputs CAMH -> Trillium.
#
# Prerequisites on CAMH:
#   bash stage_remaining_cohorts_on_camh.sh --cohort GVEX
#
# Run on Trillium (tmux):
#   bash transfer_gvex.sh
#
# Dest: /project/rrg-shreejoy/GVEX
# Size: ~12G normalized genotypes + RNA/meta (+ ~469G raw if included)
# Raw BrainGVEX dose VCFs (~469G) live under external_data and must be staged
# onto netdata first:  bash stage_remaining_cohorts_on_camh.sh --cohort GVEX_RAW
# Set SKIP_RAW=1 to transfer normalized-only.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/GVEX}"
STAGE_GVEX="${STAGE_ROOT}/GVEX"
VCF_DIR="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/WGS/QC/BrainGVEX_combined_vcf_normalized"
VCF_RAW_STAGED="${STAGE_GVEX}/Genotype_raw"   # staged from external_data (~469G)
SKIP_RAW="${SKIP_RAW:-0}"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p "${DEST_ROOT}/RNA" "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Genotype/BrainGVEX_combined_vcf_normalized"

trillium_rsync "gvex_rna" \
  "${STAGE_GVEX}/RNA/" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "gvex_meta" \
  "${STAGE_GVEX}/Metadata/" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "gvex_vcf" \
  "${VCF_DIR}/" \
  "${DEST_ROOT}/Genotype/BrainGVEX_combined_vcf_normalized/"

# Raw dose VCFs (~469G) — must be staged first (stage_remaining_cohorts_on_camh.sh --cohort GVEX_RAW)
if [[ "${SKIP_RAW}" -eq 1 ]]; then
  echo "SKIP_RAW=1 — skipping GVEX raw genotype transfer"
else
  trillium_rsync "gvex_vcf_raw" \
    "${VCF_RAW_STAGED}/" \
    "${DEST_ROOT}/Genotype/BrainGVEX_raw_dose_vcf/"
fi

echo "GVEX transfer finished -> ${DEST_ROOT}"
echo "  RNA/GVEX_count_matrix.csv"
echo "  Metadata/ (manifest, clinical, ID mapping)"
echo "  Genotype/BrainGVEX_combined_vcf_normalized/   # ~12G"
if [[ "${SKIP_RAW}" -ne 1 ]]; then
  echo "  Genotype/BrainGVEX_raw_dose_vcf/              # ~469G raw (staged)"
fi
exit 0
