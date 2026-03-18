#!/bin/bash
# ============================================================================
# GVEX combined VCF pre-processing: split multiallelics and bgzip/index output
#
# Input:
#   /external/rprshnas01/external_data/psychencode/PsychENCODE/genotypes_BrainGVEX/DNA/BrainGVEX_chr{chrom}.vcf
#
# Output:
#   /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/BrainGVEX_combined_vcf_normalized/GVEX.QC.{chrom}.normalized.vcf.gz
# ============================================================================

set -euo pipefail

BCFTOOLS="${BCFTOOLS:-/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools}"
GVEX_BASE="/external/rprshnas01/external_data/psychencode/PsychENCODE/genotypes_BrainGVEX/DNA"
OUT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/BrainGVEX_combined_vcf_normalized"
LOG_DIR="${OUT_DIR}/logs"
CHROMS=$(seq 1 22)
SLURM_TIME="${SLURM_TIME:-4:00:00}"
SLURM_MEM="${SLURM_MEM:-16G}"
SLURM_CPUS="${SLURM_CPUS:-2}"

process_one_chrom() {
    local CHROM=$1
    local IN_VCF="${GVEX_BASE}/BrainGVEX_chr${CHROM}.vcf"
    local OUT_VCF="${OUT_DIR}/GVEX.QC.${CHROM}.normalized.vcf.gz"

    mkdir -p "${OUT_DIR}" "${LOG_DIR}"

    if [[ -f "${OUT_VCF}" && -f "${OUT_VCF}.tbi" ]]; then
        echo "[SKIP] chr${CHROM}"
        return 0
    fi

    echo "[RUN] chr${CHROM}"
    "${BCFTOOLS}" norm -m- "${IN_VCF}" -Oz -o "${OUT_VCF}"
    "${BCFTOOLS}" index -t "${OUT_VCF}"
    echo "[DONE] chr${CHROM}"
}

if [[ "${1:-}" == "slurm" ]]; then
    mkdir -p "${OUT_DIR}" "${LOG_DIR}"
    sbatch \
        --job-name="prep_gvex" \
        --array=1-22 \
        --time="${SLURM_TIME}" \
        --mem="${SLURM_MEM}" \
        --cpus-per-task="${SLURM_CPUS}" \
        --output="${LOG_DIR}/prep_gvex_%A_%a.log" \
        --error="${LOG_DIR}/prep_gvex_%A_%a.err" \
        --wrap="bash \"$0\" \${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

if [[ -n "${1:-}" ]]; then
    process_one_chrom "$1"
    exit 0
fi

for CHROM in ${CHROMS}; do
    process_one_chrom "${CHROM}"
done
