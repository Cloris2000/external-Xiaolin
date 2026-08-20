#!/usr/bin/env bash
# Run on CAMH (full shell, e.g. dev03) — NOT on Trillium.
#
# Stages files that live outside netdata_kcni (nethome / external_data / kcni)
# into netdata so the Trillium rsync gateway can see them.
#
#   bash stage_remaining_cohorts_on_camh.sh
#   bash stage_remaining_cohorts_on_camh.sh --cohort MSBB
#   bash stage_remaining_cohorts_on_camh.sh --cohort all

set -euo pipefail

STAGE_ROOT="${STAGE_ROOT:-/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data_input/trillium_staging}"
COHORT="all"

usage() {
  cat <<EOF
Usage: $(basename "$0") [--cohort NAME]

  NAME: all (default) | shared | MSBB | CMC | NABEC | GVEX
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort) COHORT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

stage_shared() {
  echo "[$(date '+%F %T')] Staging shared refs..."
  mkdir -p "${STAGE_ROOT}/shared/tools"
  cp -v /external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt \
        "${STAGE_ROOT}/shared/"
  cp -v /external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/new_MTGnCgG_lfct2.5_Publication.csv \
        "${STAGE_ROOT}/shared/"
  # Tools under /kcni may not be visible to the Trillium rsync gateway
  cp -v /external/rprshnas01/kcni/mwainberg/software/plink2 \
        /external/rprshnas01/kcni/mwainberg/software/regenie \
        "${STAGE_ROOT}/shared/tools/"
  # liftover chain + METAL already under netdata-visible paths
}

stage_msbb() {
  echo "[$(date '+%F %T')] Staging MSBB metadata..."
  mkdir -p "${STAGE_ROOT}/MSBB/Metadata"
  cp -v /nethome/kcni/xzhou/GWAS_tut/AMP-AD/MSBB_biospecimen_metadata.csv \
        "${STAGE_ROOT}/MSBB/Metadata/"
}

stage_cmc() {
  echo "[$(date '+%F %T')] Staging CMC RNA + SNP metadata..."
  mkdir -p "${STAGE_ROOT}/CMC/RNA" "${STAGE_ROOT}/CMC/Metadata"
  cp -v /nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_MSSM_count_matrix.csv \
        /nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_MSSM_metadata.csv \
        /nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_PENN_count_matrix.csv \
        /nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_PENN_metadata.csv \
        /nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_PITT_count_matrix.csv \
        /nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_PITT_metadata.csv \
        "${STAGE_ROOT}/CMC/RNA/"
  cp -v /external/rprshnas01/netdata_kcni/stlab/CMC_genotypes/SNPs/Release3/Metadata/CMC_Human_SNP_metadata.csv \
        "${STAGE_ROOT}/CMC/Metadata/"
}

stage_nabec() {
  echo "[$(date '+%F %T')] Staging NABEC RNA/metadata..."
  mkdir -p "${STAGE_ROOT}/NABEC/RNA" "${STAGE_ROOT}/NABEC/Metadata"
  cp -v /nethome/kcni/xzhou/GWAS_tut/NABEC/gene_count_matrix_geneid.csv \
        /nethome/kcni/xzhou/GWAS_tut/NABEC/NABEC_metadata_combined.csv \
        /nethome/kcni/xzhou/GWAS_tut/NABEC/combined_metrics.csv \
        "${STAGE_ROOT}/NABEC/RNA/"
  # Biospecimen mapping expected at nextflow/data/NABEC_biospecimen_mapping.txt
  # — may be missing; copy if present.
  local biospec="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data/NABEC_biospecimen_mapping.txt"
  if [[ -f "${biospec}" ]]; then
    cp -v "${biospec}" "${STAGE_ROOT}/NABEC/Metadata/"
  else
    echo "WARNING: ${biospec} not found — locate/regenerate before NABEC GWAS on Trillium"
  fi
}

stage_gvex() {
  echo "[$(date '+%F %T')] Staging GVEX RNA/metadata from external_data + kcni..."
  mkdir -p "${STAGE_ROOT}/GVEX/RNA" "${STAGE_ROOT}/GVEX/Metadata"
  cp -v /external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/GVEX_count_matrix.csv \
        "${STAGE_ROOT}/GVEX/RNA/"
  cp -v /external/rprshnas01/external_data/psychencode/PsychENCODE/BrainGVEX/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv \
        /external/rprshnas01/external_data/psychencode/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv \
        /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data_input/gvex/gvex_rna_wgs_id_mapping.tsv \
        "${STAGE_ROOT}/GVEX/Metadata/"
}

case "${COHORT}" in
  all)
    stage_shared; stage_msbb; stage_cmc; stage_nabec; stage_gvex
    ;;
  shared) stage_shared ;;
  MSBB|msbb) stage_msbb ;;
  CMC|cmc) stage_cmc ;;
  NABEC|nabec) stage_nabec ;;
  GVEX|gvex) stage_gvex ;;
  *) echo "Unknown cohort: ${COHORT}" >&2; usage; exit 1 ;;
esac

echo
echo "Staging done under ${STAGE_ROOT}"
du -sh "${STAGE_ROOT}"/* 2>/dev/null || true
