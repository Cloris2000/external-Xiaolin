#!/usr/bin/env bash
#SBATCH --job-name=psychad_radc_download
#SBATCH --partition=long
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/psychad_radc_download_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/psychad_radc_download_%j.err
#
# Download PsychAD combined VCF (syn60552875) for RADC cohort,
# verify sample overlap with the 101-donor RADC proportion file,
# subset to the matched donors, and index the result.
#
# Outputs:
#   WGS_DIR/combined.vcf.gz                  — original download
#   WGS_DIR/combined.vcf.gz.tbi
#   WGS_DIR/psychad_radc_samples.txt         — all R* IDs in the VCF
#   WGS_DIR/psychad_radc_89_samples.txt      — RADC donors matched to proportion file
#   WGS_DIR/psychad_radc_subset.vcf.gz       — subset VCF (89 donors)
#   WGS_DIR/psychad_radc_subset.vcf.gz.tbi
#
# IMPORTANT — genome build check:
#   The existing pipeline uses GRCh37/b37. The PsychAD combined VCF is
#   expected to be GRCh38/b38. After download this script prints the build
#   from the VCF header. If b38, the genotype QC pipeline will need to be
#   run with a b38 reference and the GWAS summary stats will require
#   liftover to b37 before meta-analysis (scripts/liftover_sumstats.py).
#
# Usage:
#   sbatch scripts/download_psychad_radc_vcf.sbatch
#   # or interactively:
#   bash scripts/download_psychad_radc_vcf.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
WGS_DIR=${NF_DIR}/data_input/sn_psychad_radc_hodge/genotype
PROP_FILE=${NF_DIR}/data/snRNAseq_hodge_label_cell_prop/psychad_radc_cell_proportions_qc.csv
SYNAPSE=/nethome/kcni/xzhou/.anaconda3/bin/synapse
BCFTOOLS=/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools
TABIX=/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/tabix

mkdir -p "${WGS_DIR}" "${NF_DIR}/logs"

echo "============================================"
echo "PsychAD RADC VCF download + prep"
echo "$(date)"
echo "WGS_DIR: ${WGS_DIR}"
echo "============================================"

# ── 1. Download combined.vcf.gz from Synapse ─────────────────────────────────
echo ""
echo "[1/6] Downloading combined.vcf.gz (syn60552875)..."
cd "${WGS_DIR}"
${SYNAPSE} get syn60552875 --downloadLocation "${WGS_DIR}"

VCF="${WGS_DIR}/combined.vcf.gz"
[ -f "${VCF}" ] || { echo "ERROR: combined.vcf.gz not found after download"; exit 1; }

# ── 2. Index ──────────────────────────────────────────────────────────────────
echo ""
echo "[2/6] Indexing VCF..."
${TABIX} -p vcf "${VCF}"

# ── 3. Genome build check ─────────────────────────────────────────────────────
echo ""
echo "[3/6] Genome build check:"
${BCFTOOLS} view -h "${VCF}" | grep -i 'reference\|assembly\|contig.*length\|GRCh\|hg38\|hg19\|b37\|b38' \
  | head -5 || echo "  (no explicit build tag in header — check contig lengths manually)"
echo ""
echo "  *** If GRCh38/hg38: GWAS summary stats will need liftover to b37    ***"
echo "  *** before meta-analysis with the bulk CTP cohorts (all b37).        ***"
echo "  *** Use: scripts/liftover_sumstats.py after the sn GWAS run.         ***"

# ── 4. Verify sample overlap ──────────────────────────────────────────────────
echo ""
echo "[4/6] Verifying sample overlap with RADC proportion file..."
${BCFTOOLS} query -l "${VCF}" > "${WGS_DIR}/psychad_radc_samples.txt"
TOTAL_VCF=$(wc -l < "${WGS_DIR}/psychad_radc_samples.txt")
R_STAR=$(grep -c '^R' "${WGS_DIR}/psychad_radc_samples.txt" || true)
echo "  Total samples in VCF: ${TOTAL_VCF}"
echo "  R* samples in VCF:    ${R_STAR}"

# Get the 101 RADC proportion IDs
awk -F',' 'NR>1{print $1}' "${PROP_FILE}" | sort > "${WGS_DIR}/prop_ids.txt"
sort "${WGS_DIR}/psychad_radc_samples.txt" > "${WGS_DIR}/vcf_ids_sorted.txt"

comm -12 "${WGS_DIR}/vcf_ids_sorted.txt" "${WGS_DIR}/prop_ids.txt" \
  > "${WGS_DIR}/psychad_radc_89_samples.txt"
N_MATCH=$(wc -l < "${WGS_DIR}/psychad_radc_89_samples.txt")
echo "  Matched RADC donors (proportion ∩ VCF): ${N_MATCH} / 101"

if [ "${N_MATCH}" -lt 10 ]; then
  echo "ERROR: fewer than 10 matched donors — check ID format in VCF vs proportion file"
  exit 1
fi

# ── 5. Subset VCF to matched RADC donors ─────────────────────────────────────
echo ""
echo "[5/6] Subsetting VCF to ${N_MATCH} matched RADC donors..."
${BCFTOOLS} view \
  -S "${WGS_DIR}/psychad_radc_89_samples.txt" \
  -Oz \
  -o "${WGS_DIR}/psychad_radc_subset.vcf.gz" \
  "${VCF}"

${TABIX} -p vcf "${WGS_DIR}/psychad_radc_subset.vcf.gz"

# ── 6. Summary ────────────────────────────────────────────────────────────────
echo ""
echo "[6/6] Done. Outputs:"
ls -lh "${WGS_DIR}"/psychad_radc_subset.vcf.gz \
        "${WGS_DIR}"/psychad_radc_89_samples.txt \
        "${WGS_DIR}"/psychad_radc_samples.txt
echo ""
echo "Next steps:"
echo "  1. Check build from step [3] above."
echo "     - If GRCh38: run genotype QC with a b38 reference genome."
echo "     - Add liftover step after sn GWAS via scripts/liftover_sumstats.py."
echo "  2. Run genotype QC pipeline on psychad_radc_subset.vcf.gz:"
echo "     nextflow run genotyping_qc_pipeline.nf \\"
echo "       -c nextflow.config.genotyping_qc.psychad_radc"
echo "  3. Run sn GWAS:"
echo "     nextflow run combined_pipeline_v2.nf \\"
echo "       -c nextflow.config.combined.sn_psychad_radc_hodge"
echo ""
echo "============================================"
echo "Finished: $(date)"
echo "============================================"
