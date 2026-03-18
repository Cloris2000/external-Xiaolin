#!/bin/bash
# ============================================================================
# HBCC VCF Pre-processing PER PLATFORM (Release 3): Normalize + Add 0_ prefix
#
# Processes each array (h650, 1M, Omni5M) separately — no merge.
# For per-platform GWAS + meta-analysis. Output files named
# CMC_HBCC.QC.{chrom}.normalized.vcf.gz so genotyping_qc_pipeline finds them
# via normalized_vcf_dir.
#
# Usage:
#   ./prep_hbcc_vcf.sh [h650|1M|Omni5M]      # one platform, all 22 chr (sequential)
#   ./prep_hbcc_vcf.sh <platform> <chrom>    # one chromosome only (for SLURM array)
#   ./prep_hbcc_vcf.sh all                   # all three platforms (default)
#   ./prep_hbcc_vcf.sh slurm                 # submit 3 SLURM array jobs (22 chr in parallel per platform)
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

get_vcf_pattern() {
    case "$1" in
        h650)  echo "${HBCC_BASE}/h650/CMC_HBCC_h650_ImputationHRC_chr%s.dose.vcf.gz" ;;
        1M)    echo "${HBCC_BASE}/1M/CMC_HBCC_1M_ImputationHRC_chr%s.dose.vcf.gz" ;;
        Omni5M) echo "${HBCC_BASE}/Omni5M/CMC_HBCC_Omni5M_ImputationHRC_chr%s.dose.vcf.gz" ;;
        *)     echo "Unknown platform: $1" >&2; return 1 ;;
    esac
}

# Process a single chromosome for one platform (used for parallel SLURM array)
process_one_chrom() {
    local PLATFORM=$1
    local CHROM=$2
    local OUT_DIR="${BASE_OUT}/imputed_chr_normalized_${PLATFORM}"
    local WORK_DIR="${OUT_DIR}/tmp"
    mkdir -p "$WORK_DIR"

    local VCF_PATTERN
    VCF_PATTERN=$(get_vcf_pattern "$PLATFORM") || exit 1
    local IN_VCF OUT_VCF
    printf -v IN_VCF "$VCF_PATTERN" "$CHROM"
    OUT_VCF="${OUT_DIR}/CMC_HBCC.QC.${CHROM}.normalized.vcf.gz"

    if [[ -f "${OUT_VCF}" && -f "${OUT_VCF}.tbi" ]]; then
        echo "[SKIP] ${PLATFORM} chr${CHROM}"
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
    echo "  ✓ ${PLATFORM} chr${CHROM}"
}

process_platform() {
    local PLATFORM=$1
    local OUT_DIR="${BASE_OUT}/imputed_chr_normalized_${PLATFORM}"
    mkdir -p "${OUT_DIR}/tmp"
    echo "=== Platform: ${PLATFORM} → ${OUT_DIR} ==="
    for CHROM in $CHROMS; do
        process_one_chrom "$PLATFORM" "$CHROM"
    done
    echo "  Done: ${OUT_DIR}"
}

# Main
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${BASE_OUT}/logs"
SLURM_TIME="${SLURM_TIME:-4:00:00}"
SLURM_MEM="${SLURM_MEM:-16G}"
SLURM_CPUS="${SLURM_CPUS:-2}"

if [[ "${1:-}" == "slurm" ]]; then
    mkdir -p "$LOG_DIR"
    # Submit 3 array jobs (one per platform), each with 22 tasks = chr 1–22 in parallel
    for p in h650 1M Omni5M; do
        sbatch --job-name="hbcc_${p}" \
               --array=1-22 \
               --time="${SLURM_TIME}" \
               --mem="${SLURM_MEM}" \
               --cpus-per-task="${SLURM_CPUS}" \
               --output="${LOG_DIR}/prep_${p}_%A_%a.log" \
               --error="${LOG_DIR}/prep_${p}_%A_%a.err" \
               --wrap="bash '${SCRIPT_DIR}/prep_hbcc_vcf.sh' ${p} \${SLURM_ARRAY_TASK_ID}"
        echo "  Submitted: hbcc_${p} (array 1-22)"
    done
    echo "Logs: ${LOG_DIR}/prep_<platform>_<jobid>_<arrayid>.log"
    exit 0
fi

# Two args: platform + single chromosome (for SLURM array tasks)
if [[ -n "${1:-}" && -n "${2:-}" ]]; then
    process_one_chrom "$1" "$2"
    exit 0
fi

if [[ "${1:-}" == "all" || -z "${1:-}" ]]; then
    for p in h650 1M Omni5M; do
        process_platform "$p"
    done
elif [[ -n "${1:-}" ]]; then
    process_platform "$1"
else
    echo "Usage: $0 [h650|1M|Omni5M] [chrom] | all | slurm" >&2
    exit 1
fi
