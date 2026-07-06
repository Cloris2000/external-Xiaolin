#!/bin/bash
#SBATCH --job-name=gtex_vcf_split
#SBATCH --array=1-22
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=mediumtmp
#SBATCH --time=4:00:00
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/gtex_vcf_split_chr%a_%j.log
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/gtex_vcf_split_chr%a_%j.err

CHR=${SLURM_ARRAY_TASK_ID}

EXTRACTED_VCF="/external/rprshnas01/netdata_kcni/stlab/GTEx_v10/Genotype/extracted/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz"
SPLIT_DIR="/external/rprshnas01/netdata_kcni/stlab/GTEx_v10/Genotype/split_by_chr"
BCFTOOLS="/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools"
THREADS=8

OUT_VCF="${SPLIT_DIR}/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv_chr${CHR}.vcf.gz"

mkdir -p "$SPLIT_DIR"

if [[ -f "$OUT_VCF" ]]; then
    echo "chr${CHR}: already exists, skipping"
    exit 0
fi

echo "chr${CHR}: starting split at $(date)"
echo "  Input:  $EXTRACTED_VCF"
echo "  Output: $OUT_VCF"

$BCFTOOLS view \
    --regions "chr${CHR}" \
    --output-type z \
    --threads "$THREADS" \
    --output-file "$OUT_VCF" \
    "$EXTRACTED_VCF"

echo "chr${CHR}: indexing..."
$BCFTOOLS index --tbi "$OUT_VCF"

echo "chr${CHR}: done at $(date)"
echo "  Size: $(ls -lh $OUT_VCF | awk '{print $5}')"
