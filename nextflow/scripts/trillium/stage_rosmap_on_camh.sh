#!/usr/bin/env bash
# Run on CAMH (full shell, e.g. dev03) — NOT on Trillium.
#
# The Trillium rsync gateway (192.197.205.74) can see netdata_kcni but NOT
# external_data or (reliably) nethome. This script copies ROSMAP inputs that
# live outside netdata into:
#   .../nextflow/data_input/trillium_staging/ROSMAP/
#
# Then Trillium transfer_rosmap.sh pulls only from netdata-visible paths.

set -euo pipefail

STAGE="${STAGE:-/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data_input/trillium_staging/ROSMAP}"
EXT_META="/external/rprshnas01/external_data/rosmap/metadata"
TOPMED_VCF="/external/rprshnas01/external_data/rosmap/genotype/TOPmed_imputed/vcf"
BIOSPEC="/nethome/kcni/xzhou/GWAS_tut/AMP-AD/ROSMAP_biospecimen_metadata.csv"
ARRAY_RESULTS="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP_array"

STAGE_META=1
STAGE_TOPMED=1

usage() {
  cat <<EOF
Usage: $(basename "$0") [--metadata-only | --topmed-only]

  --metadata-only   Stage small metadata/sample maps only (~4MB)
  --topmed-only     Stage TOPmed imputed VCFs only (~433G)
  (default)         Stage both
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --metadata-only) STAGE_TOPMED=0; shift ;;
    --topmed-only) STAGE_META=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "${STAGE}/Metadata" \
         "${STAGE}/Genotype/ROSMAP_array_sample_maps" \
         "${STAGE}/Genotype/ROSMAP_array_TOPmed_imputed_vcf"

if [[ "${STAGE_META}" -eq 1 ]]; then
  echo "[$(date '+%F %T')] Staging metadata..."
  cp -v \
    "${EXT_META}/ROSmaster.rds" \
    "${EXT_META}/ROSMAP_assay_wholeGenomeSeq_metadata.csv" \
    "${EXT_META}/ROSMAP_assay_snpArray_metadata.csv" \
    "${BIOSPEC}" \
    "${STAGE}/Metadata/"

  cp -v \
    "${ARRAY_RESULTS}/samples_keep_original_ids.txt" \
    "${ARRAY_RESULTS}/samples_reheader_map.txt" \
    "${ARRAY_RESULTS}/samples_to_keep.txt" \
    "${STAGE}/Genotype/ROSMAP_array_sample_maps/"
  echo "[$(date '+%F %T')] Metadata staging done."
fi

if [[ "${STAGE_TOPMED}" -eq 1 ]]; then
  echo "[$(date '+%F %T')] Staging TOPmed VCFs (~433G) into netdata..."
  echo "  SRC : ${TOPMED_VCF}/"
  echo "  DEST: ${STAGE}/Genotype/ROSMAP_array_TOPmed_imputed_vcf/"
  rsync -avP "${TOPMED_VCF}/" "${STAGE}/Genotype/ROSMAP_array_TOPmed_imputed_vcf/"
  echo "[$(date '+%F %T')] TOPmed staging done."
fi

echo
echo "Staging complete under: ${STAGE}"
du -sh "${STAGE}" "${STAGE}"/* 2>/dev/null || true
