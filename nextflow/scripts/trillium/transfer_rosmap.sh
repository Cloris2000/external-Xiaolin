#!/usr/bin/env bash
# Transfer ROSMAP + ROSMAP_array pipeline inputs from CAMH -> Trillium.
#
# Prerequisites on CAMH (staging of metadata + TOPmed raw + sample maps):
#   bash stage_rosmap_on_camh.sh
#
# Run on Trillium (tmux):
#   bash transfer_rosmap.sh
#   bash transfer_rosmap.sh --dry-run
#
# Transfers:
#   1. RNA pipeline inputs: ROSMAP_DLPFC_batch_all.csv (count matrix, ~448M) +
#      ROSMAP_combined_metrics.csv (metadata). Shared by ROSMAP & ROSMAP_array.
#   2. Metadata (staged): biospecimen, ROSmaster.rds, assay snpArray/wholeGenomeSeq
#   3. ROSMAP joint WGS raw VCFs           (~426G)
#   4. ROSMAP joint WGS normalized VCFs    (~437G)
#   5. ROSMAP_array raw TOPmed imputed VCFs (~370G, staged; SKIP_TOPMED=1 to skip)
#   6. ROSMAP_array normalized VCFs        (~439M)
#   7. ROSMAP_array sample ID maps
#
# Note: bulk RNA under /project/rrg-shreejoy/ROSMAP/ROSMAP_Raw_Counts_Bulk and
# ROSMAP_RNAseq_provenance are separate provenance files (not the pipeline
# count_matrix_file), so the pipeline RNA inputs above are still transferred.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/ROSMAP}"
META_PCA="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/metabrain_PCA/data"
WGS_BASE="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/WGS"
STAGE="${STAGE_ROOT}/ROSMAP"
STAGE_META="${STAGE}/Metadata"
STAGE_TOPMED="${STAGE}/Genotype/ROSMAP_array_TOPmed_imputed_vcf"
STAGE_SAMPLE_MAPS="${STAGE}/Genotype/ROSMAP_array_sample_maps"
SKIP_TOPMED="${SKIP_TOPMED:-0}"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]

  SKIP_TOPMED=1  skip the ~370G ROSMAP_array raw TOPmed imputed VCFs
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

# 0) RNA pipeline inputs (count matrix + metrics) — on netdata, direct
trillium_rsync "rosmap_counts" \
  "${META_PCA}/ROSMAP_DLPFC_batch_all.csv" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "rosmap_metrics" \
  "${META_PCA}/ROSMAP_combined_metrics.csv" \
  "${DEST_ROOT}/RNA/"

# 1) Metadata (staged on CAMH from external_data + nethome -> netdata)
trillium_rsync "metadata" \
  "${STAGE_META}/" \
  "${DEST_ROOT}/Metadata/"

# 2) ROSMAP WGS genotypes — raw (already on netdata)
trillium_rsync "wgs_raw" \
  "${WGS_BASE}/ROSMAP_joint_WGS_vcf/" \
  "${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf/"

# 3) ROSMAP WGS genotypes — normalized (already on netdata)
trillium_rsync "wgs_normalized" \
  "${WGS_BASE}/ROSMAP_joint_WGS_vcf_normalized/" \
  "${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf_normalized/"

# 4) ROSMAP_array genotypes — raw TOPmed (staged onto netdata on CAMH first)
if [[ "${SKIP_TOPMED}" -eq 1 ]]; then
  echo "SKIP_TOPMED=1 — skipping TOPmed raw transfer"
else
  trillium_rsync "array_topmed_raw" \
    "${STAGE_TOPMED}/" \
    "${DEST_ROOT}/Genotype/ROSMAP_array_TOPmed_imputed_vcf/"
fi

# 5) ROSMAP_array genotypes — normalized (already on netdata)
trillium_rsync "array_normalized" \
  "${WGS_BASE}/ROSMAP_array_vcf_normalized/" \
  "${DEST_ROOT}/Genotype/ROSMAP_array_vcf_normalized/"

# 6) Sample ID maps (staged)
trillium_rsync "array_sample_maps" \
  "${STAGE_SAMPLE_MAPS}/" \
  "${DEST_ROOT}/Genotype/ROSMAP_array_sample_maps/"

echo "All ROSMAP transfers finished -> ${DEST_ROOT}"
echo "  RNA/ROSMAP_DLPFC_batch_all.csv          # pipeline count matrix"
echo "  RNA/ROSMAP_combined_metrics.csv         # pipeline metadata"
echo "  Metadata/ (biospecimen, ROSmaster.rds, assay snpArray/wholeGenomeSeq)"
echo "  Genotype/ROSMAP_joint_WGS_vcf/                 # WGS raw"
echo "  Genotype/ROSMAP_joint_WGS_vcf_normalized/      # WGS normalized"
echo "  Genotype/ROSMAP_array_TOPmed_imputed_vcf/      # array raw (TOPmed)"
echo "  Genotype/ROSMAP_array_vcf_normalized/          # array normalized"
echo "  Genotype/ROSMAP_array_sample_maps/             # ID maps for prep"
