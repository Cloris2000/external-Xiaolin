#!/usr/bin/env bash
#SBATCH --job-name=psychad_radc_split
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-22
#SBATCH --output=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/psychad_radc_split_%A_%a.out
#SBATCH --error=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs/psychad_radc_split_%A_%a.err
#SBATCH --account=rrg-shreejoy
#
# Split PsychAD combined.vcf.gz (GRCh38, chr*-prefixed contigs) into
# per-chromosome VCFs with numeric contig names (1..22) for the
# genotyping QC pipeline.
#
# Usage:
#   sbatch scripts/split_psychad_radc_vcf_by_chr.sbatch

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
GENO=${NF_DIR}/data_input/sn_psychad_radc_hodge/genotype
OUT_DIR=${GENO}/by_chr
VCF=${GENO}/combined.vcf.gz
RENAME=${GENO}/chr_rename_hg38_to_numeric.txt
BCFTOOLS=/home/zhoux156/miniforge3/envs/bcftools_env/bin/bcftools
TABIX=/home/zhoux156/miniforge3/envs/bcftools_env/bin/tabix

CHROM=${SLURM_ARRAY_TASK_ID}
OUT=${OUT_DIR}/PsychAD_RADC.chr${CHROM}.vcf.gz

mkdir -p "${OUT_DIR}" "${NF_DIR}/logs"

echo "[$(date)] chr${CHROM}: extracting + renaming contigs"
echo "  Host: $(hostname)"

[ -f "${VCF}.tbi" ] || { echo "ERROR: ${VCF}.tbi missing — index combined.vcf.gz first"; exit 1; }
[ -f "${RENAME}" ] || { echo "ERROR: rename map missing: ${RENAME}"; exit 1; }

# Extract chrN region, rename chrN -> N for pipeline ID convention
${BCFTOOLS} view -r "chr${CHROM}" "${VCF}" \
  | ${BCFTOOLS} annotate --rename-chrs "${RENAME}" -Oz -o "${OUT}" \
  --threads 4

${TABIX} -p vcf "${OUT}"

echo "[$(date)] chr${CHROM}: done → ${OUT} ($(ls -lh ${OUT} | awk '{print $5}'))"
${BCFTOOLS} query -l "${OUT}" | wc -l | awk '{print "  samples:", $1}'
${BCFTOOLS} view -h "${OUT}" | grep -E '^##contig=<ID=' | head -2
