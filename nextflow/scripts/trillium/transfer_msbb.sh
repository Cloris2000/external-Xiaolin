#!/usr/bin/env bash
# Transfer MSBB cohort inputs CAMH -> Trillium.
#
# Prerequisites on CAMH:
#   bash stage_remaining_cohorts_on_camh.sh --cohort MSBB
#
# Run on Trillium (tmux):
#   bash transfer_msbb.sh
#
# Dest: /project/rrg-shreejoy/MSBB
# Size: ~240G genotypes + RNA/meta

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/MSBB}"
META_PCA="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data"
WGS_BASE="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS"
STAGE_MSBB="${STAGE_ROOT}/MSBB"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p "${DEST_ROOT}/RNA" "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Genotype/MSBB_joint_WGS_vcf" \
  "${DEST_ROOT}/Genotype/MSBB_joint_WGS_vcf_normalized"

trillium_rsync "msbb_counts" \
  "${META_PCA}/MSBB_gene_all_counts_matrix_clean_gene_id_aggregated.csv" \
  "${DEST_ROOT}/RNA/"

# Raw / intermediate RNA matrices that generate the aggregated count matrix:
#   clean_gene_id.csv  -> raw per-sample counts (dup gene IDs)
#   _unique.csv        -> dedup step
#   _aggregated.csv    -> gene-level aggregate (pipeline input, above)
trillium_rsync "msbb_counts_raw" \
  "${META_PCA}/MSBB_gene_all_counts_matrix_clean_gene_id.csv" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "msbb_counts_unique" \
  "${META_PCA}/MSBB_gene_all_counts_matrix_clean_gene_id_unique.csv" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "msbb_meta" \
  "${META_PCA}/msbb_meta.csv" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "msbb_biospec" \
  "${STAGE_MSBB}/Metadata/" \
  "${DEST_ROOT}/Metadata/"

# Genotypes are large (~240G) and already transferred. Set SKIP_WGS=1 to
# skip re-scanning them (e.g. when only pulling the raw RNA matrices).
if [[ "${SKIP_WGS:-0}" -eq 1 ]]; then
  echo "SKIP_WGS=1 — skipping MSBB WGS genotype transfers"
else
  trillium_rsync "msbb_wgs_raw" \
    "${WGS_BASE}/MSBB_joint_WGS_vcf/" \
    "${DEST_ROOT}/Genotype/MSBB_joint_WGS_vcf/"

  trillium_rsync "msbb_wgs_norm" \
    "${WGS_BASE}/MSBB_joint_WGS_vcf_normalized/" \
    "${DEST_ROOT}/Genotype/MSBB_joint_WGS_vcf_normalized/"
fi

echo "MSBB transfer finished -> ${DEST_ROOT}"
echo "  RNA/ = aggregated + clean_gene_id (raw per-sample) + unique matrices"
echo "  Metadata/ (incl. biospecimen)"
if [[ "${SKIP_WGS:-0}" -ne 1 ]]; then
  echo "  Genotype/MSBB_joint_WGS_vcf/              # raw ~119G"
  echo "  Genotype/MSBB_joint_WGS_vcf_normalized/   # normalized ~121G"
fi
