#!/bin/bash
# ============================================================================
# HBCC VCF Pre-processing (Release 3): Normalize + Add 0_ prefix
#
# Two modes:
#   Per-platform: processes each array (h650, 1M, Omni5M) separately.
#     For per-platform GWAS + meta-analysis.
#     Output: .../imputed_chr_normalized_{platform}/
#
#   Merged: bcftools merge on the 3 platforms' raw imputed VCFs per
#     chromosome (1M=242 + Omni5M=109 + h650=104 = 455 samples total),
#     then normalizes. All samples in one VCF. For single-cohort GWAS.
#     Output: .../imputed_chr_normalized/
#
# Usage:
#   ./prep_hbcc_vcf.sh [h650|1M|Omni5M]      # one platform, all 22 chr (sequential)
#   ./prep_hbcc_vcf.sh <platform> <chrom>    # one chromosome only (for SLURM array)
#   ./prep_hbcc_vcf.sh all                   # all three platforms separately (default)
#   ./prep_hbcc_vcf.sh slurm                 # submit 3 SLURM array jobs (per-platform)
#   ./prep_hbcc_vcf.sh merged                # merge all 3 platforms, all 22 chr
#   ./prep_hbcc_vcf.sh merged <chrom>        # merge all 3 platforms, one chr (for SLURM array)
#   ./prep_hbcc_vcf.sh merged slurm          # submit SLURM array job for merged mode
#
# Output dirs:
#   .../imputed_chr_normalized_h650/
#   .../imputed_chr_normalized_1M/
#   .../imputed_chr_normalized_Omni5M/
#   .../imputed_chr_normalized/              (merged mode)
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

# ── Merged mode: merge all 3 platforms per chromosome then normalize ─────────

process_merged_chrom() {
    local CHROM=$1
    local OUT_DIR="${BASE_OUT}/imputed_chr_normalized"
    local WORK_DIR="${OUT_DIR}/tmp"
    mkdir -p "$WORK_DIR"

    local OUT_VCF="${OUT_DIR}/CMC_HBCC.QC.${CHROM}.normalized.vcf.gz"

    if [[ -f "${OUT_VCF}" && -f "${OUT_VCF}.tbi" ]]; then
        echo "[SKIP] merged chr${CHROM}"
        return 0
    fi

    local VCF_1M VCF_OMNI VCF_H650
    printf -v VCF_1M    "${HBCC_BASE}/1M/CMC_HBCC_1M_ImputationHRC_chr%s.dose.vcf.gz"     "$CHROM"
    printf -v VCF_OMNI  "${HBCC_BASE}/Omni5M/CMC_HBCC_Omni5M_ImputationHRC_chr%s.dose.vcf.gz" "$CHROM"
    printf -v VCF_H650  "${HBCC_BASE}/h650/CMC_HBCC_h650_ImputationHRC_chr%s.dose.vcf.gz"  "$CHROM"

    local MERGED_VCF="${WORK_DIR}/chr${CHROM}.merged.vcf.gz"
    local NORM_VCF="${WORK_DIR}/chr${CHROM}.norm.vcf.gz"
    local SAMPLES_FILE="${WORK_DIR}/chr${CHROM}_samples.txt"
    local RENAMED_SAMPLES="${WORK_DIR}/chr${CHROM}_renamed.txt"

    local THREADS="${SLURM_CPUS_PER_TASK:-${SLURM_CPUS:-2}}"

    # Merge 3 platforms (non-overlapping samples, same variants from HRC imputation)
    $BCFTOOLS merge --missing-to-ref --threads "$THREADS" \
        -Oz -o "$MERGED_VCF" \
        "$VCF_1M" "$VCF_OMNI" "$VCF_H650"
    $BCFTOOLS index -t "$MERGED_VCF"

    # Normalize (split multi-allelic)
    $BCFTOOLS norm -m- --threads "$THREADS" "$MERGED_VCF" -Oz -o "$NORM_VCF" 2>/dev/null
    $BCFTOOLS index -t "$NORM_VCF"

    # Add 0_ prefix to sample IDs
    $BCFTOOLS query -l "$NORM_VCF" > "$SAMPLES_FILE"
    paste "$SAMPLES_FILE" <(awk '{print "0_"$0}' "$SAMPLES_FILE") > "$RENAMED_SAMPLES"
    $BCFTOOLS reheader -s "$RENAMED_SAMPLES" -o "$OUT_VCF" "$NORM_VCF"
    $BCFTOOLS index -t "$OUT_VCF"

    rm -f "$MERGED_VCF" "$MERGED_VCF".tbi "$NORM_VCF" "$NORM_VCF".tbi \
          "$SAMPLES_FILE" "$RENAMED_SAMPLES"
    echo "  ✓ merged chr${CHROM}"
}

process_merged_all() {
    local OUT_DIR="${BASE_OUT}/imputed_chr_normalized"
    mkdir -p "${OUT_DIR}/tmp"
    echo "=== Merged (all platforms) → ${OUT_DIR} ==="
    for CHROM in $CHROMS; do
        process_merged_chrom "$CHROM"
    done
    rmdir "${OUT_DIR}/tmp" 2>/dev/null || true
    echo "  Done: ${OUT_DIR}"
}

# ── Main ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${BASE_OUT}/logs"
SLURM_TIME="${SLURM_TIME:-4:00:00}"
SLURM_MEM="${SLURM_MEM:-16G}"
SLURM_CPUS="${SLURM_CPUS:-2}"

# merged slurm — submit one array job (22 chr in parallel)
if [[ "${1:-}" == "merged" && "${2:-}" == "slurm" ]]; then
    mkdir -p "$LOG_DIR"
    sbatch --job-name="hbcc_merged" \
           --array=1-22 \
           --time="${SLURM_TIME}" \
           --mem="${SLURM_MEM}" \
           --cpus-per-task="${SLURM_CPUS}" \
           --output="${LOG_DIR}/prep_merged_%A_%a.log" \
           --error="${LOG_DIR}/prep_merged_%A_%a.err" \
           --wrap="bash '${SCRIPT_DIR}/prep_hbcc_vcf.sh' merged \${SLURM_ARRAY_TASK_ID}"
    echo "Submitted: hbcc_merged (array 1-22)"
    echo "Logs: ${LOG_DIR}/prep_merged_<jobid>_<arrayid>.log"
    exit 0
fi

# merged <chrom> — single chromosome for SLURM array task
if [[ "${1:-}" == "merged" && -n "${2:-}" && "${2:-}" =~ ^[0-9]+$ ]]; then
    process_merged_chrom "$2"
    exit 0
fi

# merged — all chromosomes sequentially
if [[ "${1:-}" == "merged" ]]; then
    process_merged_all
    exit 0
fi

# per-platform slurm
if [[ "${1:-}" == "slurm" ]]; then
    mkdir -p "$LOG_DIR"
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
    echo "Usage: $0 [h650|1M|Omni5M] [chrom] | all | slurm | merged | merged <chrom> | merged slurm" >&2
    exit 1
fi
