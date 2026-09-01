#!/usr/bin/env bash
# Transfer shared pipeline references/tools CAMH -> Trillium.
#
# Prerequisites on CAMH:
#   bash stage_remaining_cohorts_on_camh.sh --cohort shared
#
# Run on Trillium:
#   bash transfer_shared_refs.sh
#
# Dest: /project/rrg-shreejoy/pipeline_refs

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/pipeline_refs}"
STAGE_SHARED="${STAGE_ROOT}/shared"
REF_DATA="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/nextflow/reference_data"
METAL="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/WGS/METAL/generic-metal/executables"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p "${DEST_ROOT}/markers" "${DEST_ROOT}/reference_data" \
  "${DEST_ROOT}/tools/metal" "${DEST_ROOT}/tools"

trillium_rsync "hgnc" \
  "${STAGE_SHARED}/hgnc_complete_set.txt" \
  "${DEST_ROOT}/markers/"

trillium_rsync "mgp_markers" \
  "${STAGE_SHARED}/new_MTGnCgG_lfct2.5_Publication.csv" \
  "${DEST_ROOT}/markers/"

trillium_rsync "liftover_chain" \
  "${REF_DATA}/" \
  "${DEST_ROOT}/reference_data/"

trillium_rsync "metal" \
  "${METAL}/" \
  "${DEST_ROOT}/tools/metal/"

trillium_rsync "plink_regenie" \
  "${STAGE_SHARED}/tools/" \
  "${DEST_ROOT}/tools/"

echo "Shared refs transfer finished -> ${DEST_ROOT}"
echo "  markers/hgnc_complete_set.txt"
echo "  markers/new_MTGnCgG_lfct2.5_Publication.csv"
echo "  reference_data/ (hg38ToHg19 chain)"
echo "  tools/plink2, tools/regenie, tools/metal/"
echo "Note: /project/rrg-shreejoy/Genomic_references already has genome fasta/GENCODE — not duplicated."
