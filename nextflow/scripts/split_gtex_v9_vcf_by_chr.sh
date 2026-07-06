#!/usr/bin/env bash
# split_gtex_v9_vcf_by_chr.sh
#
# Extracts the GTEx v9 WGS combined VCF from the dbGaP TAR archive and
# splits it into per-chromosome VCF files suitable for the Nextflow pipeline.
#
# The dbGaP TAR contains a single whole-genome VCF:
#   GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz
#
# Steps:
#   1. Rename .tmp -> .tar once download is verified complete
#   2. Extract the VCF from the TAR
#   3. Split into chr1..chr22 VCFs using bcftools
#   4. Index each output VCF
#
# Prerequisites:
#   - Complete download of phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1.GRU.tar
#   - bcftools available on PATH (or update BCFTOOLS below)
#   - ~1.2 TB free disk space (extracted VCF + per-chr VCFs)
#
# Usage:
#   bash split_gtex_v9_vcf_by_chr.sh [--check-only]
#   --check-only: only verify the download, do not extract or split

set -euo pipefail

# ============================================================================
# CONFIGURATION — edit paths as needed
# ============================================================================
GENOTYPE_DIR="/external/rprshnas01/netdata_kcni/stlab/GTEx_v10/Genotype"
TAR_TMP="${GENOTYPE_DIR}/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1.GRU.tar.tmp"
TAR_FINAL="${GENOTYPE_DIR}/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1.GRU.tar"

VCF_INNER_PATH="./phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz"
VCF_BASENAME="GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz"

EXTRACT_DIR="${GENOTYPE_DIR}/extracted/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1"
SPLIT_DIR="${GENOTYPE_DIR}/split_by_chr"

BCFTOOLS="${BCFTOOLS:-/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools}"

CHROMOSOMES=$(seq 1 22)
THREADS=8   # bcftools view threads per chromosome; adjust to available CPUs

# ============================================================================
# 0. Parse arguments
# ============================================================================
CHECK_ONLY=false
for arg in "$@"; do
  [[ "$arg" == "--check-only" ]] && CHECK_ONLY=true
done

# ============================================================================
# 1. Verify download completeness
# ============================================================================
echo "============================================================================="
echo "Step 1: Verifying download"
echo "============================================================================="

if [[ ! -f "$TAR_TMP" && ! -f "$TAR_FINAL" ]]; then
  echo "ERROR: TAR file not found at:"
  echo "  $TAR_TMP"
  echo "  $TAR_FINAL"
  echo ""
  echo "To resume download, run from $GENOTYPE_DIR:"
  echo "  module load sratoolkit  # or use conda env with sra-tools"
  echo "  prefetch --ngc prj_32731.ngc --cart cart_prj32731_202604141354.krt"
  exit 1
fi

TAR_FILE="${TAR_FINAL}"
if [[ ! -f "$TAR_FINAL" && -f "$TAR_TMP" ]]; then
  TAR_FILE="$TAR_TMP"
fi

echo "TAR file: $TAR_FILE"
TAR_SIZE=$(du -sh "$TAR_FILE" | cut -f1)
echo "Size: $TAR_SIZE"

# Test TAR integrity
echo "Testing TAR integrity (this reads the full archive - may take several minutes)..."
if tar -tf "$TAR_FILE" > /dev/null 2>&1; then
  echo "TAR integrity: OK"
  DOWNLOAD_OK=true
else
  echo ""
  echo "WARNING: TAR integrity check FAILED (likely incomplete download)."
  echo "  Error: $(tar -tf "$TAR_FILE" 2>&1 | tail -3)"
  echo ""
  echo "To resume/re-download, navigate to $GENOTYPE_DIR and run:"
  echo "  prefetch --ngc prj_32731.ngc --cart cart_prj32731_202604141354.krt"
  echo ""
  DOWNLOAD_OK=false
fi

if [[ "$CHECK_ONLY" == "true" ]]; then
  echo ""
  echo "Check-only mode. Exiting without extraction."
  exit 0
fi

if [[ "$DOWNLOAD_OK" != "true" ]]; then
  echo "ERROR: Download incomplete. Complete the download before extracting."
  exit 1
fi

# ============================================================================
# 2. Rename .tmp -> .tar (if needed) and extract VCF
# ============================================================================
echo ""
echo "============================================================================="
echo "Step 2: Extracting VCF from TAR"
echo "============================================================================="

if [[ ! -f "$TAR_FINAL" && -f "$TAR_TMP" ]]; then
  echo "Renaming $TAR_TMP -> $TAR_FINAL"
  mv "$TAR_TMP" "$TAR_FINAL"
  TAR_FILE="$TAR_FINAL"
fi

mkdir -p "$EXTRACT_DIR"

EXTRACTED_VCF="${EXTRACT_DIR}/${VCF_BASENAME}"
if [[ -f "$EXTRACTED_VCF" ]]; then
  echo "VCF already extracted: $EXTRACTED_VCF  (skipping extraction)"
else
  echo "Extracting: $VCF_INNER_PATH"
  echo "To: $EXTRACT_DIR"
  tar -xf "$TAR_FILE" -C "$EXTRACT_DIR" \
      --strip-components=1 \
      "$VCF_INNER_PATH"
  echo "Extraction done: $EXTRACTED_VCF"
fi

# Verify extraction
if [[ ! -f "$EXTRACTED_VCF" ]]; then
  echo "ERROR: Extracted VCF not found at $EXTRACTED_VCF"
  exit 1
fi

# Index if needed
VCF_INDEX="${EXTRACTED_VCF}.tbi"
if [[ ! -f "$VCF_INDEX" ]]; then
  echo "Indexing extracted VCF..."
  $BCFTOOLS index --tbi --threads "$THREADS" "$EXTRACTED_VCF"
  echo "Indexing done."
fi

# ============================================================================
# 3. Check available chromosomes
# ============================================================================
echo ""
echo "============================================================================="
echo "Step 3: Checking chromosomes in combined VCF"
echo "============================================================================="

echo "Chromosomes present in VCF header:"
$BCFTOOLS view --header-only "$EXTRACTED_VCF" | grep "^##contig" | head -30

# ============================================================================
# 4. Split by chromosome
# ============================================================================
echo ""
echo "============================================================================="
echo "Step 4: Splitting by chromosome"
echo "============================================================================="
mkdir -p "$SPLIT_DIR"

for CHR in $CHROMOSOMES; do
  # Support both "chr1" and "1" contig naming
  OUT_VCF="${SPLIT_DIR}/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv_chr${CHR}.vcf.gz"

  if [[ -f "$OUT_VCF" ]]; then
    echo "  chr${CHR}: already exists, skipping"
    continue
  fi

  echo "  chr${CHR}: splitting..."
  # Try chr-prefixed first, fall back to numeric
  if $BCFTOOLS view --header-only "$EXTRACTED_VCF" 2>/dev/null | grep -q "##contig=<ID=chr${CHR},"; then
    REGION="chr${CHR}"
  else
    REGION="${CHR}"
  fi

  $BCFTOOLS view \
    --regions "$REGION" \
    --output-type z \
    --threads "$THREADS" \
    --output-file "$OUT_VCF" \
    "$EXTRACTED_VCF"

  $BCFTOOLS index --tbi "$OUT_VCF"
  echo "  chr${CHR}: done -> $OUT_VCF"
done

echo ""
echo "============================================================================="
echo "All chromosomes split successfully."
echo ""
echo "Update vcf_pattern in nextflow.config.combined.gtex_v10 to:"
echo "  ${SPLIT_DIR}/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv_chr{chrom}.vcf.gz"
echo "============================================================================="
