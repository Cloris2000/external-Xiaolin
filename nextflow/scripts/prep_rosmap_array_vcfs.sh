#!/bin/bash
# =============================================================================
# prep_rosmap_array_vcfs.sh
#
# Splits the merged_overlap_rs BED (all 2067 samples, R2>0.8, MAF>0.0025,
# overlapping SNPs from both Illumina and Affymetrix platforms) into
# per-chromosome VCFs restricted to the 170 ROSMAP_array samples.
#
# Two-step process per chromosome:
#   1. plink2 --keep (original BED IDs: MAP#####/11AD#####) -> raw VCF
#   2. bcftools reheader to rename samples to Study+projid FIDs
#      (e.g. 11AD39717_11AD39717 -> MAP45800230)
#      so the pipeline's --double-id step produces FID=IID=Study+projid
#
# Output files follow the pipeline's expected naming convention:
#   {normalized_vcf_dir}/{study}.QC.{chrom}.normalized.vcf.gz
# =============================================================================

set -euo pipefail

PLINK=/external/rprshnas01/kcni/mwainberg/software/plink2
BCFTOOLS=/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools

BED=/external/rprshnas01/external_data/rosmap/genotype/TOPmed_imputed/vcf/merged/merged_overlap_rs
RESULTS_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP_array
OUTDIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_array_vcf_normalized
STUDY=ROSMAP_array

# Keep file uses ORIGINAL BED IDs (topmed_specimenID) so all 170 match
KEEP="${RESULTS_DIR}/samples_keep_original_ids.txt"
# Reheader map: old "SPEC_SPEC" -> new "Study+projid" (for bcftools reheader)
REHEADER="${RESULTS_DIR}/samples_reheader_map.txt"

mkdir -p "$OUTDIR"

echo "================================================"
echo "Splitting merged_overlap_rs -> per-chr VCFs"
echo "Samples: $(wc -l < "$KEEP")"
echo "Output:  $OUTDIR"
echo "================================================"

for CHR in $(seq 1 22); do
    OUT_FINAL="${OUTDIR}/${STUDY}.QC.${CHR}.normalized.vcf.gz"

    if [[ -f "$OUT_FINAL" ]] && [[ -f "${OUT_FINAL}.tbi" ]]; then
        N_SAM=$($BCFTOOLS query -l "$OUT_FINAL" | wc -l)
        echo "Chr ${CHR}: already exists (${N_SAM} samples), skipping"
        continue
    fi

    TMP_VCF="${OUTDIR}/${STUDY}.QC.${CHR}.tmp.vcf.gz"

    echo "Chr ${CHR}: exporting from BED..."
    $PLINK \
        --bfile "$BED" \
        --chr "$CHR" \
        --keep "$KEEP" \
        --export vcf bgz \
        --out "${OUTDIR}/${STUDY}.QC.${CHR}.tmp" \
        --no-psam-pheno \
        --threads 4

    echo "Chr ${CHR}: reheadering sample IDs to Study+projid FIDs..."
    $BCFTOOLS reheader \
        --samples "$REHEADER" \
        --output "$OUT_FINAL" \
        "$TMP_VCF"

    $BCFTOOLS index -t "$OUT_FINAL"

    # Clean up tmp file
    rm -f "$TMP_VCF" "${OUTDIR}/${STUDY}.QC.${CHR}.tmp.vcf.gz.tbi" \
          "${OUTDIR}/${STUDY}.QC.${CHR}.tmp.log" \
          "${OUTDIR}/${STUDY}.QC.${CHR}.tmp.pvar" \
          "${OUTDIR}/${STUDY}.QC.${CHR}.tmp.psam"

    N_SAMPLES=$($BCFTOOLS query -l "$OUT_FINAL" | wc -l)
    FIRST_SAMPLE=$($BCFTOOLS query -l "$OUT_FINAL" | head -1)
    echo "  Chr ${CHR}: ${N_SAMPLES} samples (e.g. ${FIRST_SAMPLE})"
done

echo ""
echo "================================================"
echo "Done. Summary:"
echo "================================================"
N_DONE=$(ls "${OUTDIR}/${STUDY}.QC."*.normalized.vcf.gz 2>/dev/null | wc -l)
echo "${N_DONE} chromosome VCFs written"
echo ""
echo "Sample IDs in chr1 VCF (first 5):"
$BCFTOOLS query -l "${OUTDIR}/${STUDY}.QC.1.normalized.vcf.gz" | head -5
