#!/bin/bash
# ============================================================================
# HBCC VCF Pre-processing PER PLATFORM (Release 3): Normalize + Add 0_ prefix
#
# For per-platform GWAS + meta-analysis: process each array separately
# (no merging). Output files named CMC_HBCC.QC.{chrom}.normalized.vcf.gz
# so genotyping_qc_pipeline finds them via normalized_vcf_dir.
#
# Usage:
#   ./prep_hbcc_vcf_per_platform.sh [h650|1M|Omni5M]
#   ./prep_hbcc_vcf_per_platform.sh all    # run all three platforms
#
# Output dirs:
#   .../imputed_chr_normalized_h650/
#   .../imputed_chr_normalized_1M/
#   .../imputed_chr_normalized_Omni5M/
# ============================================================================

set -euo pipefail

BCFTOOLS="${BCFTOOLS:-/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools}"
HBCC_BASE="/external/rprshnas01/netdata_kcni/stlab/CMC_genotypes/SNPs/Release3/Imputed/HBCC"
BASE_OUT="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/CMC_HBCC_reimputed"
CHROMS=$(seq 1 22)
USE_SLURM=${USE_SLURM:-false}

process_platform() {
    local PLATFORM=$1
    local OUT_DIR="${BASE_OUT}/imputed_chr_normalized_${PLATFORM}"
    local WORK_DIR="${OUT_DIR}/tmp"
    mkdir -p "$WORK_DIR"

    case "$PLATFORM" in
        h650)
            VCF_PATTERN="${HBCC_BASE}/h650/CMC_HBCC_h650_ImputationHRC_chr%s.dose.vcf.gz"
            ;;
        1M)
            VCF_PATTERN="${HBCC_BASE}/1M/CMC_HBCC_1M_ImputationHRC_chr%s.dose.vcf.gz"
            ;;
        Omni5M)
            VCF_PATTERN="${HBCC_BASE}/Omni5M/CMC_HBCC_Omni5M_ImputationHRC_chr%s.dose.vcf.gz"
            ;;
        *)
            echo "Unknown platform: $PLATFORM. Use h650, 1M, or Omni5M." >&2
            exit 1
            ;;
    esac

    echo "=== Platform: ${PLATFORM} → ${OUT_DIR} ==="

    process_chrom() {
        local CHROM=$1
        local IN_VCF
        printf -v IN_VCF "$VCF_PATTERN" "$CHROM"
        local OUT_VCF="${OUT_DIR}/CMC_HBCC.QC.${CHROM}.normalized.vcf.gz"

        if [[ -f "${OUT_VCF}" && -f "${OUT_VCF}.tbi" ]]; then
            echo "  [SKIP] chr${CHROM}"
            return 0
        fi

        local NORM_VCF="${WORK_DIR}/chr${CHROM}.norm.vcf.gz"
        local SAMPLES_FILE="${WORK_DIR}/chr${CHROM}_samples.txt"
        local RENAMED_SAMPLES="${WORK_DIR}/chr${CHROM}_renamed.txt"

        $BCFTOOLS norm -m- "$IN_VCF" -Oz -o "$NORM_VCF" 2>/dev/null
        $BCFTOOLS index -t "$NORM_VCF"
        $BCFTOOLS query -l "$NORM_VCF" > "$SAMPLES_FILE"
        paste "$SAMPLES_FILE" <(awk '{print "0_"$0}' "$SAMPLES_FILE") > "$RENAMED_SAMPLES"
        $BCFTOOLS reheader -s "$RENAMED_SAMPLES" -o "$OUT_VCF" "$NORM_VCF"
        $BCFTOOLS index -t "$OUT_VCF"
        rm -f "$NORM_VCF" "$NORM_VCF".tbi "$SAMPLES_FILE" "$RENAMED_SAMPLES"
        echo "  ✓ chr${CHROM}"
    }

    for CHROM in $CHROMS; do
        process_chrom "$CHROM"
    done
    echo "  Done: ${OUT_DIR}"
}

# Main
if [[ "${1:-}" == "all" ]]; then
    for p in h650 1M Omni5M; do
        process_platform "$p"
    done
elif [[ -n "${1:-}" ]]; then
    process_platform "$1"
else
    echo "Usage: $0 <h650|1M|Omni5M|all>" >&2
    exit 1
fi
