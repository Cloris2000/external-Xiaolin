#!/bin/bash
#SBATCH --job-name=hbcc_merge_vcf
#SBATCH --partition=medium
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/hbcc_merge_vcf_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/hbcc_merge_vcf_%j.err

set -euo pipefail

source ~/.anaconda3/etc/profile.d/conda.sh && conda activate bcftools_env

DIR_1M="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/CMC_HBCC_reimputed/imputed_chr_normalized_1M"
DIR_OMNI="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/CMC_HBCC_reimputed/imputed_chr_normalized_Omni5M"
DIR_H650="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/CMC_HBCC_reimputed/imputed_chr_normalized_h650"
OUT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/CMC_HBCC_reimputed/imputed_chr_normalized"

echo "Merging HBCC VCFs (1M + Omni5M + h650) into: $OUT_DIR"
echo "Started: $(date)"

for chrom in $(seq 1 22); do
    vcf_1M="${DIR_1M}/CMC_HBCC.QC.${chrom}.normalized.vcf.gz"
    vcf_omni="${DIR_OMNI}/CMC_HBCC.QC.${chrom}.normalized.vcf.gz"
    vcf_h650="${DIR_H650}/CMC_HBCC.QC.${chrom}.normalized.vcf.gz"
    out="${OUT_DIR}/CMC_HBCC.QC.${chrom}.normalized.vcf.gz"

    echo "  chr${chrom}: merging..."
    bcftools merge \
        --missing-to-ref \
        --threads 4 \
        -Oz -o "$out" \
        "$vcf_1M" "$vcf_omni" "$vcf_h650"
    bcftools index --tbi "$out"
    echo "  chr${chrom}: done ($(bcftools query -l $out | wc -l) samples)"
done

echo "Merge complete: $(date)"
echo "Total files in output dir: $(ls $OUT_DIR/*.vcf.gz | wc -l)"
