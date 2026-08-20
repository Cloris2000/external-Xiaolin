#!/usr/bin/env bash
#SBATCH --job-name=psychad_radc_prep
#SBATCH --partition=mediumtmp
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/psychad_radc_prep_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/psychad_radc_prep_%j.err
#
# Resume PsychAD RADC VCF prep AFTER download is complete.
# Skips indexing if combined.vcf.gz.tbi already exists.
#
# Prerequisites:
#   data_input/sn_psychad_radc_hodge/genotype/combined.vcf.gz
#
# Usage:
#   # Option A — wait for interactive tabix to finish, then:
#   sbatch scripts/prep_psychad_radc_vcf.sbatch
#
#   # Option B — kill interactive tabix and redo indexing on SLURM:
#   # (Ctrl-C the interactive script, then)
#   sbatch scripts/prep_psychad_radc_vcf.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
WGS_DIR=${NF_DIR}/data_input/sn_psychad_radc_hodge/genotype
PROP_FILE=${NF_DIR}/data/snRNAseq_hodge_label_cell_prop/psychad_radc_cell_proportions_qc.csv
BCFTOOLS=/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools
TABIX=/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/tabix

VCF="${WGS_DIR}/combined.vcf.gz"
mkdir -p "${NF_DIR}/logs"

echo "============================================"
echo "PsychAD RADC VCF prep (post-download)"
echo "$(date)"
echo "Host: $(hostname)"
echo "WGS_DIR: ${WGS_DIR}"
echo "============================================"

[ -f "${VCF}" ] || { echo "ERROR: ${VCF} not found — download first"; exit 1; }
echo "VCF size: $(ls -lh ${VCF} | awk '{print $5}')"

# ── 1. Index (skip if already done) ──────────────────────────────────────────
echo ""
if [ -f "${VCF}.tbi" ]; then
    echo "[1/4] Index exists — skipping tabix"
    ls -lh "${VCF}.tbi"
else
    echo "[1/4] Indexing with tabix..."
    ${TABIX} -p vcf "${VCF}"
    ls -lh "${VCF}.tbi"
fi

# ── 2. Genome build check ────────────────────────────────────────────────────
echo ""
echo "[2/4] Genome build (from VCF header):"
${BCFTOOLS} view -h "${VCF}" \
    | grep -i 'reference\|assembly\|GRCh\|hg38\|hg19\|b37\|b38' \
    | head -8 \
    || echo "  (no explicit build tag — check contig lengths manually)"
echo ""
echo "  Contig lengths (chr1 / 1):"
${BCFTOOLS} view -h "${VCF}" | grep -E '^##contig=<ID=(chr)?1,' | head -2

# ── 3. Sample overlap ────────────────────────────────────────────────────────
echo ""
echo "[3/4] Checking sample overlap with proportion file..."
${BCFTOOLS} query -l "${VCF}" > "${WGS_DIR}/psychad_radc_samples.txt"
TOTAL_VCF=$(wc -l < "${WGS_DIR}/psychad_radc_samples.txt")
echo "  Total samples in VCF: ${TOTAL_VCF}"

awk -F',' 'NR>1{print $1}' "${PROP_FILE}" | sort > "${WGS_DIR}/prop_ids.txt"
sort "${WGS_DIR}/psychad_radc_samples.txt" > "${WGS_DIR}/vcf_ids_sorted.txt"
comm -12 "${WGS_DIR}/vcf_ids_sorted.txt" "${WGS_DIR}/prop_ids.txt" \
    > "${WGS_DIR}/psychad_radc_matched.txt"
N_MATCH=$(wc -l < "${WGS_DIR}/psychad_radc_matched.txt")
echo "  Matched RADC donors (proportion ∩ VCF): ${N_MATCH} / 101"

if [ "${N_MATCH}" -lt 10 ]; then
    echo "ERROR: fewer than 10 matched donors — check ID formats"
    echo "  First 5 VCF IDs: $(head -5 ${WGS_DIR}/psychad_radc_samples.txt | tr '\n' ' ')"
    echo "  First 5 prop IDs: $(head -5 ${WGS_DIR}/prop_ids.txt | tr '\n' ' ')"
    exit 1
fi

# ── 4. Subset ────────────────────────────────────────────────────────────────
echo ""
echo "[4/4] Subsetting VCF to ${N_MATCH} matched donors..."
${BCFTOOLS} view \
    -S "${WGS_DIR}/psychad_radc_matched.txt" \
    --threads 4 \
    -Oz \
    -o "${WGS_DIR}/psychad_radc_subset.vcf.gz" \
    "${VCF}"
${TABIX} -p vcf "${WGS_DIR}/psychad_radc_subset.vcf.gz"

echo ""
echo "============================================"
echo "Done: $(date)"
ls -lh "${WGS_DIR}"/psychad_radc_subset.vcf.gz \
        "${WGS_DIR}"/psychad_radc_matched.txt \
        "${WGS_DIR}"/psychad_radc_samples.txt \
        "${VCF}.tbi"
echo "============================================"
