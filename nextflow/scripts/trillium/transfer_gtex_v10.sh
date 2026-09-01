#!/usr/bin/env bash
# Transfer GTEx_v10 cohort inputs CAMH -> Trillium.
# All on netdata_kcni (no staging). Largest cohort by far.
#
# Genotype note: GTEx v9 WGS VCFs are ALREADY normalized by the consortium
# (config: normalize_vcf=false, normalized_vcf_dir=""). The pipeline consumes
# split_by_chr/ directly, so those per-chr VCFs ARE the normalized input.
#
# Volumes (SKIP_RAW=1 skips the two huge raw archives, keeping split_by_chr):
#   split_by_chr/                     ~766G  (pipeline input / normalized)
#   extracted combined VCF            ~715G  (raw, untarred single VCF)
#   dbGaP .GRU.tar                    ~847G  (raw original archive)
#   genotype-qc + sample-info tars    ~927M
# Full run ~= 2.3 TB — ensure destination has space; run in tmux.
#
# Run on Trillium (tmux):
#   bash transfer_gtex_v10.sh              # everything (raw + normalized)
#   SKIP_RAW=1 bash transfer_gtex_v10.sh   # split_by_chr only
#
# Dest: /project/rrg-shreejoy/GTEx_v10

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=trillium_rsync_lib.sh
source "${SCRIPT_DIR}/trillium_rsync_lib.sh"

DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/GTEx_v10}"
RNA_DIR="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/data"
GENO_BASE="/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/GTEx_v10/Genotype"
VCF_DIR="${GENO_BASE}/split_by_chr"
EXTRACTED_DIR="${GENO_BASE}/extracted/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1"
RAW_TAR="${GENO_BASE}/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1.GRU.tar"
QC_TAR="${GENO_BASE}/phg001796.v1.GTEx_v9.genotype-qc.MULTI.tar"
SAMPLEINFO_TAR="${GENO_BASE}/phg001796.v1.GTEx_v9.sample-info.MULTI.tar"
SKIP_RAW="${SKIP_RAW:-0}"

trillium_usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]

  SKIP_RAW=1   skip the ~1.5T raw archives (dbGaP tar + extracted combined VCF),
               transferring only split_by_chr/ (the normalized pipeline input)
EOF
}

trillium_parse_args "$@"
trillium_init "${DEST_ROOT}"

mkdir -p "${DEST_ROOT}/RNA" "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Genotype/split_by_chr" \
  "${DEST_ROOT}/Genotype/raw_dbGaP" \
  "${DEST_ROOT}/Genotype/extracted_combined_vcf"

trillium_rsync "gtex_v10_counts" \
  "${RNA_DIR}/GTEx_v10_BA9_raw_count_matrix.csv" \
  "${DEST_ROOT}/RNA/"

trillium_rsync "gtex_v10_meta" \
  "${RNA_DIR}/GTEx_v10_BA9_sample_metadata.csv" \
  "${DEST_ROOT}/Metadata/"

# Processed / normalized (pipeline input) — per-chr VCFs
trillium_rsync "gtex_v10_vcf" \
  "${VCF_DIR}/" \
  "${DEST_ROOT}/Genotype/split_by_chr/"

# Small dbGaP metadata archives (always transferred)
trillium_rsync "gtex_v10_qc_tar" \
  "${QC_TAR}" \
  "${DEST_ROOT}/Genotype/raw_dbGaP/"

trillium_rsync "gtex_v10_sampleinfo_tar" \
  "${SAMPLEINFO_TAR}" \
  "${DEST_ROOT}/Genotype/raw_dbGaP/"

if [[ "${SKIP_RAW}" -eq 1 ]]; then
  echo "SKIP_RAW=1 — skipping dbGaP .GRU.tar (~847G) and extracted combined VCF (~715G)"
else
  # Raw original dbGaP archive
  trillium_rsync "gtex_v10_raw_tar" \
    "${RAW_TAR}" \
    "${DEST_ROOT}/Genotype/raw_dbGaP/"

  # Raw extracted single combined VCF (pre-split)
  trillium_rsync "gtex_v10_extracted_vcf" \
    "${EXTRACTED_DIR}/" \
    "${DEST_ROOT}/Genotype/extracted_combined_vcf/"
fi

echo "GTEx_v10 transfer finished -> ${DEST_ROOT}"
echo "  RNA/GTEx_v10_BA9_raw_count_matrix.csv"
echo "  Metadata/GTEx_v10_BA9_sample_metadata.csv"
echo "  Genotype/split_by_chr/            # ~766G  normalized pipeline input (per-chr VCF)"
echo "  Genotype/raw_dbGaP/               # dbGaP .GRU.tar (~847G) + qc/sample-info tars"
echo "  Genotype/extracted_combined_vcf/  # ~715G  raw untarred single VCF"
