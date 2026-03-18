#!/bin/bash
#
# Debug Meta-Analysis Differences Script
# Compares original and Nextflow pipeline results
#
# Usage: bash debug_meta_analysis_differences.sh

set -e

echo "============================================================================="
echo "Meta-Analysis Debugging Script"
echo "============================================================================="
echo ""

# Define paths
ORIGINAL_PHENO="/nethome/kcni/xzhou/GWAS_tut/AMP-AD2/1438_AMP-AD_MGP_estimations_techAdj_sexAge_RINT.txt"
NEXTFLOW_PHENO="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/phenotypes_RINT.txt"
ORIGINAL_COV="/nethome/kcni/xzhou/GWAS_tut/ROSMAP/rosmap_wgs_cov.txt"
NEXTFLOW_COV="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/covariates.txt"
ORIGINAL_GWAS="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2/ROSMAP_WGS_step2_update_SST.regenie.raw_p"
NEXTFLOW_GWAS="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/step2/ROSMAP_SST_step2.regenie.raw_p"
ORIGINAL_META="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/METAL/Joint_AMP_AD_meta_with_CMC_NABEC_GTEx/SST_meta_analysis1.tbl"
NEXTFLOW_META="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis/SST_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP.tbl"

# Create output directory
OUTPUT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/debugging_outputs"
mkdir -p ${OUTPUT_DIR}

echo "1. Comparing Sample Counts"
echo "============================================================================="
echo ""
echo "ROSMAP Phenotype Files:"
echo "  Original: $(wc -l < ${ORIGINAL_PHENO}) lines ($(expr $(wc -l < ${ORIGINAL_PHENO}) - 1) samples)"
echo "  Nextflow: $(wc -l < ${NEXTFLOW_PHENO}) lines ($(expr $(wc -l < ${NEXTFLOW_PHENO}) - 1) samples)"
echo ""
echo "ROSMAP Covariate Files:"
echo "  Original: $(wc -l < ${ORIGINAL_COV}) lines ($(expr $(wc -l < ${ORIGINAL_COV}) - 1) samples)"
echo "  Nextflow: $(wc -l < ${NEXTFLOW_COV}) lines ($(expr $(wc -l < ${NEXTFLOW_COV}) - 1) samples)"
echo ""

echo "2. Extracting Sample IDs"
echo "============================================================================="
echo ""
# Extract sample IDs from phenotype files
tail -n +2 ${ORIGINAL_PHENO} | cut -f2 | sort > ${OUTPUT_DIR}/original_pheno_samples.txt
tail -n +2 ${NEXTFLOW_PHENO} | cut -f2 | sort > ${OUTPUT_DIR}/nextflow_pheno_samples.txt

# Extract sample IDs from covariate files
tail -n +2 ${ORIGINAL_COV} | cut -f2 | sort > ${OUTPUT_DIR}/original_cov_samples.txt
tail -n +2 ${NEXTFLOW_COV} | cut -f2 | sort > ${OUTPUT_DIR}/nextflow_cov_samples.txt

echo "Sample IDs saved to:"
echo "  ${OUTPUT_DIR}/original_pheno_samples.txt"
echo "  ${OUTPUT_DIR}/nextflow_pheno_samples.txt"
echo "  ${OUTPUT_DIR}/original_cov_samples.txt"
echo "  ${OUTPUT_DIR}/nextflow_cov_samples.txt"
echo ""

echo "3. Finding Missing Samples"
echo "============================================================================="
echo ""
# Samples in original but not in nextflow (phenotypes)
comm -23 ${OUTPUT_DIR}/original_pheno_samples.txt ${OUTPUT_DIR}/nextflow_pheno_samples.txt > ${OUTPUT_DIR}/samples_missing_in_nextflow_pheno.txt
echo "Phenotype samples in original but NOT in Nextflow: $(wc -l < ${OUTPUT_DIR}/samples_missing_in_nextflow_pheno.txt)"
echo "  Saved to: ${OUTPUT_DIR}/samples_missing_in_nextflow_pheno.txt"
echo "  First 10 missing samples:"
head -10 ${OUTPUT_DIR}/samples_missing_in_nextflow_pheno.txt | sed 's/^/    /'
echo ""

# Samples in nextflow but not in original (phenotypes)
comm -13 ${OUTPUT_DIR}/original_pheno_samples.txt ${OUTPUT_DIR}/nextflow_pheno_samples.txt > ${OUTPUT_DIR}/samples_extra_in_nextflow_pheno.txt
echo "Phenotype samples in Nextflow but NOT in original: $(wc -l < ${OUTPUT_DIR}/samples_extra_in_nextflow_pheno.txt)"
echo "  Saved to: ${OUTPUT_DIR}/samples_extra_in_nextflow_pheno.txt"
if [ $(wc -l < ${OUTPUT_DIR}/samples_extra_in_nextflow_pheno.txt) -gt 0 ]; then
    echo "  First 10 extra samples:"
    head -10 ${OUTPUT_DIR}/samples_extra_in_nextflow_pheno.txt | sed 's/^/    /'
fi
echo ""

# Common samples
comm -12 ${OUTPUT_DIR}/original_pheno_samples.txt ${OUTPUT_DIR}/nextflow_pheno_samples.txt > ${OUTPUT_DIR}/common_pheno_samples.txt
echo "Common phenotype samples: $(wc -l < ${OUTPUT_DIR}/common_pheno_samples.txt)"
echo "  Saved to: ${OUTPUT_DIR}/common_pheno_samples.txt"
echo ""

echo "4. Comparing Phenotype Values for Common Samples"
echo "============================================================================="
echo ""
# Get first 5 common samples
COMMON_SAMPLES=$(head -5 ${OUTPUT_DIR}/common_pheno_samples.txt)

echo "Comparing SST phenotype values for first 5 common samples:"
echo ""
printf "%-15s %-20s %-20s %-20s\n" "Sample_ID" "Original_SST" "Nextflow_SST" "Difference"
printf "%-15s %-20s %-20s %-20s\n" "---------------" "--------------------" "--------------------" "--------------------"

for sample in ${COMMON_SAMPLES}; do
    # Get SST column (column 18) from original
    original_sst=$(grep -w "^[^\t]*\t${sample}\t" ${ORIGINAL_PHENO} | cut -f18)
    # Get SST column from nextflow (need to find which column SST is in)
    nextflow_header=$(head -1 ${NEXTFLOW_PHENO})
    sst_col=$(echo "${nextflow_header}" | tr '\t' '\n' | grep -n "^SST$" | cut -d: -f1)
    nextflow_sst=$(grep -w "^[^\t]*\t${sample}\t" ${NEXTFLOW_PHENO} | cut -f${sst_col})
    
    # Calculate difference if both values exist
    if [ ! -z "${original_sst}" ] && [ ! -z "${nextflow_sst}" ]; then
        diff=$(echo "${original_sst} - ${nextflow_sst}" | bc -l)
        printf "%-15s %-20s %-20s %-20s\n" "${sample}" "${original_sst}" "${nextflow_sst}" "${diff}"
    else
        printf "%-15s %-20s %-20s %-20s\n" "${sample}" "${original_sst:-NA}" "${nextflow_sst:-NA}" "NA"
    fi
done
echo ""

echo "5. Comparing GWAS Results for Specific SNPs"
echo "============================================================================="
echo ""
if [ -f "${NEXTFLOW_GWAS}" ]; then
    echo "Comparing chr11:12541586:A:G (from meta-analysis comparison):"
    echo ""
    echo "Original GWAS:"
    grep "chr11:12541586:A:G" ${ORIGINAL_GWAS} | head -1
    echo ""
    echo "Nextflow GWAS:"
    grep "chr11:12541586:A:G" ${NEXTFLOW_GWAS} | head -1
    echo ""
    
    echo "Comparing chr4:70553471:A:G:"
    echo ""
    echo "Original GWAS:"
    grep "chr4:70553471:A:G" ${ORIGINAL_GWAS} | head -1
    echo ""
    echo "Nextflow GWAS:"
    grep "chr4:70553471:A:G" ${NEXTFLOW_GWAS} | head -1
    echo ""
else
    echo "WARNING: Nextflow GWAS file not found at ${NEXTFLOW_GWAS}"
    echo "Searching for alternative locations..."
    find /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow -name "ROSMAP_SST_step2.regenie.raw_p" -type f 2>/dev/null | head -5
    echo ""
fi

echo "6. Comparing Meta-Analysis Results"
echo "============================================================================="
echo ""
echo "Top 10 SNPs by p-value in Original Meta-Analysis:"
tail -n +2 ${ORIGINAL_META} | sort -k10,10g | head -10 | cut -f1,10 | column -t
echo ""
echo "Top 10 SNPs by p-value in Nextflow Meta-Analysis:"
tail -n +2 ${NEXTFLOW_META} | sort -k10,10g | head -10 | cut -f1,10 | column -t
echo ""

echo "7. Comparing Specific SNPs in Meta-Analysis"
echo "============================================================================="
echo ""
echo "chr11:12541586:A:G:"
echo "  Original:"
grep "chr11:12541586:A:G" ${ORIGINAL_META} | head -1
echo "  Nextflow:"
grep "chr11:12541586:A:G" ${NEXTFLOW_META} | head -1
echo ""

echo "8. Summary Statistics"
echo "============================================================================="
echo ""
echo "Total number of SNPs:"
echo "  Original meta-analysis: $(tail -n +2 ${ORIGINAL_META} | wc -l)"
echo "  Nextflow meta-analysis: $(tail -n +2 ${NEXTFLOW_META} | wc -l)"
echo ""

echo "Genome-wide significant SNPs (p < 5e-8):"
orig_sig=$(tail -n +2 ${ORIGINAL_META} | awk '$10 < 5e-8' | wc -l)
next_sig=$(tail -n +2 ${NEXTFLOW_META} | awk '$10 < 5e-8' | wc -l)
echo "  Original meta-analysis: ${orig_sig}"
echo "  Nextflow meta-analysis: ${next_sig}"
echo ""

echo "Suggestive SNPs (p < 1e-5):"
orig_sug=$(tail -n +2 ${ORIGINAL_META} | awk '$10 < 1e-5' | wc -l)
next_sug=$(tail -n +2 ${NEXTFLOW_META} | awk '$10 < 1e-5' | wc -l)
echo "  Original meta-analysis: ${orig_sug}"
echo "  Nextflow meta-analysis: ${next_sug}"
echo ""

echo "============================================================================="
echo "Debugging complete! All outputs saved to: ${OUTPUT_DIR}"
echo "============================================================================="
echo ""
echo "Next steps:"
echo "1. Review the missing samples file to understand why they were filtered"
echo "2. Check if phenotype values differ for common samples"
echo "3. Compare the QC parameters between original and Nextflow pipelines"
echo "4. Review the detailed report: META_ANALYSIS_DEBUGGING_REPORT.md"
echo ""


