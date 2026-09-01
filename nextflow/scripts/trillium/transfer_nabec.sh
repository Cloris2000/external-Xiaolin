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
# Size: ~43G genotypes (normalized + raw PLINK/VCF/hg19 + standardized) + RNA/meta

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/NABEC}"
STAGE_NABEC="${STAGE_ROOT}/NABEC"
NABEC_WGS="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/NABEC_WGS"
VCF_DIR="${NABEC_WGS}/normalized_vcf"
# Raw genotypes (GRCh38, dbGaP freeze9):
#   PLINK2 pgen/psam/pvar = the true original dbGaP download (~20G)
#   gtonly VCF            = GT-only VCF version of the same (~5.8G)
#   hg19_vcf             = lifted-to-hg19 intermediate, input to genotyping-QC (~5.7G)
VCF_RAW_PLINK="${NABEC_WGS}/dbGaP_phs2636_NABEC_WGS_cohort_genotype_per_chromosome_PLINK"
VCF_RAW_DBGAP="${NABEC_WGS}/dbGaP_phs2636_NABEC_WGS_cohort_genotype_per_chromosome_VCF"
VCF_RAW_HG19="${NABEC_WGS}/hg19_vcf"
VCF_HG19_STD="${NABEC_WGS}/hg19_vcf_standardized"
# Sample-mapping / dbGaP file-description metadata (~224K)
SAMPLE_DIR="${NABEC_WGS}/sample"
# Pipeline metadata (on netdata, direct pull): biospecimen mapping used by config
# biospecimen_file, plus the RNA->WGS mapping source file (same content).
NF_META="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/nextflow/data/metadata"

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

# Pipeline metadata mappings (direct netdata pull, no staging needed)
trillium_rsync "nabec_biospecimen_map" \
  "${NF_META}/NABEC_biospecimen_mapping.txt" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "nabec_rna_to_wgs_map" \
  "${NF_META}/NABEC_RNA_to_WGS_mapping.txt" \
  "${DEST_ROOT}/Metadata/"

trillium_rsync "nabec_vcf" \
  "${VCF_DIR}/" \
  "${DEST_ROOT}/Genotype/normalized_vcf/"

# Raw genotypes
trillium_rsync "nabec_plink_raw" \
  "${VCF_RAW_PLINK}/" \
  "${DEST_ROOT}/Genotype/dbGaP_raw_plink/"

trillium_rsync "nabec_vcf_dbgap_raw" \
  "${VCF_RAW_DBGAP}/" \
  "${DEST_ROOT}/Genotype/dbGaP_raw_vcf/"

trillium_rsync "nabec_vcf_hg19" \
  "${VCF_RAW_HG19}/" \
  "${DEST_ROOT}/Genotype/hg19_vcf/"

trillium_rsync "nabec_vcf_hg19_std" \
  "${VCF_HG19_STD}/" \
  "${DEST_ROOT}/Genotype/hg19_vcf_standardized/"

# Sample-mapping / dbGaP file descriptions
trillium_rsync "nabec_sample_meta" \
  "${SAMPLE_DIR}/" \
  "${DEST_ROOT}/Genotype/sample/"

echo "NABEC transfer finished -> ${DEST_ROOT}"
echo "  RNA/ (counts, metadata_combined, combined_metrics)"
echo "  Metadata/NABEC_biospecimen_mapping.txt   # individualID -> specimenID (config biospecimen_file)"
echo "  Metadata/NABEC_RNA_to_WGS_mapping.txt    # RNA subject -> WGS sample (source mapping)"
echo "  Genotype/normalized_vcf/   # ~6.5G  (active GWAS input)"
echo "  Genotype/dbGaP_raw_plink/  # ~20G   raw (dbGaP freeze9 PLINK2 pgen/psam/pvar)"
echo "  Genotype/dbGaP_raw_vcf/    # ~5.8G  raw (dbGaP GT-only VCF, GRCh38)"
echo "  Genotype/hg19_vcf/              # ~5.7G  raw (hg19 liftover intermediate)"
echo "  Genotype/hg19_vcf_standardized/ # ~5.2G  standardized hg19 intermediate"
echo "  Genotype/sample/                # ~224K  sample-mapping / dbGaP file descriptions"
echo "NOTE: config reads biospecimen_file from <pipeline>/data/ — copy NABEC_biospecimen_mapping.txt there (or update the path) when running NABEC GWAS on Trillium."
