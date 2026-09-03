#!/usr/bin/env bash
#
# Download PsychAD combined VCF (syn60552875) for RADC cohort,
# verify sample overlap with the 101-donor RADC proportion file,
# subset to the matched donors, and index the result.
#
# Run interactively (NOT via SLURM -- Synapse token input required):
#   bash scripts/download_psychad_radc_vcf.sh
#
# Synapse authentication -- do ONE of:
#   Option A (recommended, persists for session):
#     synapse login -u <YOUR_EMAIL> -p <YOUR_TOKEN>
#     bash scripts/download_psychad_radc_vcf.sh
#
#   Option B (inline, token stays in shell history -- use with care):
#     SYNAPSE_AUTH_TOKEN=<YOUR_TOKEN> bash scripts/download_psychad_radc_vcf.sh
#
# Outputs:
#   WGS_DIR/combined.vcf.gz                 -- original download
#   WGS_DIR/combined.vcf.gz.tbi
#   WGS_DIR/psychad_radc_samples.txt        -- all sample IDs in the VCF
#   WGS_DIR/psychad_radc_matched.txt        -- RADC donors matched to proportion file
#   WGS_DIR/psychad_radc_subset.vcf.gz      -- subset VCF (matched donors only)
#   WGS_DIR/psychad_radc_subset.vcf.gz.tbi
#
# IMPORTANT -- genome build:
#   The existing pipeline uses GRCh37/b37. The PsychAD combined VCF is
#   expected to be GRCh38/b38. Step [3] prints the build from the VCF header.
#   If b38, GWAS summary stats will need liftover to b37 before meta-analysis
#   with the bulk CTP cohorts. Use scripts/liftover_sumstats.py afterwards.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
WGS_DIR=${NF_DIR}/data_input/sn_psychad_radc_hodge/genotype
PROP_FILE=${NF_DIR}/data/snRNAseq_hodge_label_cell_prop/psychad_radc_cell_proportions_qc.csv
SYNAPSE=/home/zhoux156/miniforge3/bin/synapse
BCFTOOLS=/home/zhoux156/miniforge3/envs/bcftools_env/bin/bcftools
TABIX=/home/zhoux156/miniforge3/envs/bcftools_env/bin/tabix

mkdir -p "${WGS_DIR}"

echo "============================================"
echo "PsychAD RADC VCF download + prep"
echo "$(date)"
echo "WGS_DIR: ${WGS_DIR}"
echo "============================================"

# ── 1. Download combined.vcf.gz from Synapse ─────────────────────────────────
echo ""
echo "[1/5] Downloading combined.vcf.gz (syn60552875)..."
echo "      (this may take a while -- ~10-50 GB depending on VCF size)"
if [ -n "${SYNAPSE_AUTH_TOKEN:-}" ]; then
    ${SYNAPSE} -p "${SYNAPSE_AUTH_TOKEN}" get syn60552875 --downloadLocation "${WGS_DIR}"
else
    ${SYNAPSE} get syn60552875 --downloadLocation "${WGS_DIR}"
fi

VCF="${WGS_DIR}/combined.vcf.gz"
[ -f "${VCF}" ] || { echo "ERROR: combined.vcf.gz not found after download"; exit 1; }
echo "  Downloaded: $(ls -lh ${VCF} | awk '{print $5}')"

# ── 2. Index ──────────────────────────────────────────────────────────────────
echo ""
echo "[2/5] Indexing..."
${TABIX} -p vcf "${VCF}"

# ── 3. Genome build check ─────────────────────────────────────────────────────
echo ""
echo "[3/5] Genome build (from VCF header):"
${BCFTOOLS} view -h "${VCF}" \
    | grep -i 'reference\|assembly\|GRCh\|hg38\|hg19\|b37\|b38' \
    | head -5 \
    || echo "  (no explicit build tag -- check contig lengths manually)"
echo ""
echo "  *** If GRCh38/hg38: liftover will be needed before meta-analysis. ***"
echo "  *** Use scripts/liftover_sumstats.py after the sn GWAS completes.  ***"

# ── 4. Verify sample overlap ──────────────────────────────────────────────────
echo ""
echo "[4/5] Checking sample overlap with proportion file..."
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
    echo "ERROR: fewer than 10 matched donors -- check ID format in VCF vs proportion file"
    exit 1
fi

# ── 5. Subset VCF to matched donors ──────────────────────────────────────────
echo ""
echo "[5/5] Subsetting VCF to ${N_MATCH} matched donors..."
${BCFTOOLS} view \
    -S "${WGS_DIR}/psychad_radc_matched.txt" \
    -Oz \
    -o "${WGS_DIR}/psychad_radc_subset.vcf.gz" \
    "${VCF}"
${TABIX} -p vcf "${WGS_DIR}/psychad_radc_subset.vcf.gz"

echo ""
echo "============================================"
echo "Done: $(date)"
ls -lh "${WGS_DIR}"/psychad_radc_subset.vcf.gz \
        "${WGS_DIR}"/psychad_radc_matched.txt \
        "${WGS_DIR}"/psychad_radc_samples.txt
echo ""
echo "Next steps:"
echo "  1. Check build from step [3] above."
echo "  2. Share the build + matched count with the AI assistant."
echo "     It will then set up genotype QC config + sn GWAS Nextflow config."
echo "============================================"
