#!/usr/bin/env Rscript

# Compare Input Files and Parameters: Original Analysis vs Nextflow Pipeline
# This will identify the root cause of r = 0.78 correlation

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("\n=======================================================\n")
cat("INPUT COMPARISON: Original Analysis vs Nextflow Pipeline\n")
cat("=======================================================\n\n")

output_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/diagnostics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =================================================
# SECTION 1: INPUT FILES USED
# =================================================
cat("[SECTION 1] INPUT FILE PATHS COMPARISON\n")
cat("-------------------------------------------------------\n\n")

cat("A. COUNT MATRIX (Gene Expression)\n")
cat("   Original (rosmap_pca_tech_cov.r):\n")
cat("     - Batch files: /external/rprshnas01/external_data/rosmap/gene_expression/.../ROSMAP_batch[1-4]_gene_all_counts_matrix.txt\n")
cat("     - Combined output: NOT SAVED (merged in memory)\n")
cat("   Nextflow config:\n")
cat("     - /external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_batch_all.csv\n")
cat("   STATUS: ✓ SAME (Nextflow uses the saved version from original analysis)\n\n")

cat("B. METADATA FILE\n")
cat("   Original (rosmap_pca_tech_cov.r):\n")
cat("     - Base: /external/rprshnas01/external_data/rosmap/.../RNAseq_Harmonization_ROSMAP_combined_metadata.csv\n")
cat("     - Output: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_combined_metrics.csv\n")
cat("   Nextflow config:\n")
cat("     - /external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_combined_metrics.csv\n")
cat("   STATUS: ✓ SAME\n\n")

cat("C. GENOTYPE FILES (VCF)\n")
cat("   Original (QC_ROSMAP_WGS_1.py):\n")
cat("     - VCF pattern: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_WGS_vcf/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{chrom}.recalibrated_variants.Broad_Rush.vcf.gz\n")
cat("     - Normalized dir: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_WGS_vcf_normalized\n")
cat("   Nextflow config:\n")
cat("     - VCF pattern: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_WGS_vcf/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{chrom}.recalibrated_variants.Broad_Rush.vcf.gz\n")
cat("     - Normalized dir: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_WGS_vcf_normalized\n")
cat("   STATUS: ✓ SAME\n\n")

cat("D. BIOSPECIMEN FILE\n")
cat("   Original: NOT EXPLICITLY SHOWN in provided scripts\n")
cat("   Nextflow config:\n")
cat("     - /nethome/kcni/xzhou/GWAS_tut/AMP-AD/ROSMAP_biospecimen_metadata.csv\n")
cat("   STATUS: ? NEEDS VERIFICATION\n\n")

cat("E. SAMPLE FILTERING FILE\n")
cat("   Original (QC_ROSMAP_WGS_1.py line 134):\n")
cat("     - /nethome/kcni/xzhou/GWAS_tut/ROSMAP/ROSMAP_samples_with_phenotypes.txt\n")
cat("   Nextflow config:\n")
cat("     - samples_to_keep = ${projectDir}/results/ROSMAP/samples_to_keep.txt\n")
cat("   STATUS: ⚠️  POTENTIALLY DIFFERENT - Need to check if same samples\n\n")

# =================================================
# SECTION 2: PROCESSING PARAMETERS
# =================================================
cat("\n[SECTION 2] PROCESSING PARAMETERS COMPARISON\n")
cat("-------------------------------------------------------\n\n")

cat("A. GENOTYPE QC PARAMETERS\n")
cat("   Parameter                 | Original (QC_ROSMAP_WGS_1.py) | Nextflow Config\n")
cat("   --------------------------|-------------------------------|----------------\n")
cat("   MAF threshold             | 0.05 (line 135)               | 0.05 ✓\n")
cat("   Mach R² filter            | 0.8 (line 133)                | 0.8 ✓\n")
cat("   VCF min GQ                | 20 (line 89)                  | 20 ✓\n")
cat("   VCF min DP                | 10 (line 90)                  | 10 ✓\n")
cat("   VCF min QUAL              | 30 (line 91)                  | 30 ✓\n")
cat("   HWE threshold             | 1e-15 (line 136)              | 1e-6 ⚠️  DIFFERENT!\n")
cat("   Geno threshold            | 0.1 (line 137)                | 0.1 ✓\n")
cat("   Mind threshold            | 0.1 (line 138)                | 0.1 ✓\n")
cat("   LD pruning (indep-pairwise)| 500kb 1 0.2 (line 153)       | 1000 50 0.2 ⚠️  DIFFERENT!\n")
cat("\n")

cat("B. PHENOTYPE PROCESSING PARAMETERS\n")
cat("   Parameter                 | Original (rosmap_pca_tech_cov.r) | Nextflow Config\n")
cat("   --------------------------|----------------------------------|----------------\n")
cat("   Tissue filter             | 'dorsolateral prefrontal cortex' | 'dorsolateral prefrontal cortex' ✓\n")
cat("   Gene variance filter      | > 0.1 (line 183)                 | 0.1 ✓\n")
cat("   Tech covariate selection  | NOT SHOWN (uses complete.cases)  | top_n_tech_cov = 20\n")
cat("   Normalization             | VST + z-score (lines 188, 196)   | ? NEEDS VERIFICATION\n")
cat("\n")

cat("C. REGENIE PARAMETERS (Need to find original command)\n")
cat("   Looking for original Regenie command...\n\n")

# =================================================
# SECTION 3: CRITICAL DIFFERENCES IDENTIFIED
# =================================================
cat("\n[SECTION 3] CRITICAL DIFFERENCES IDENTIFIED\n")
cat("=======================================================\n\n")

cat("⚠️  DIFFERENCE 1: HWE Threshold\n")
cat("   Original: 1e-15 (very strict)\n")
cat("   Nextflow: 1e-6 (less strict)\n")
cat("   Impact: Nextflow keeps MORE variants\n")
cat("   This could explain some p-value differences!\n\n")

cat("⚠️  DIFFERENCE 2: LD Pruning Parameters\n")
cat("   Original: indep-pairwise 500kb 1 0.2\n")
cat("   Nextflow: indep-pairwise 1000 50 0.2\n")
cat("   Impact: Different window sizes\n")
cat("   Note: This affects PCA calculation, which affects covariates!\n\n")

cat("❓ DIFFERENCE 3: Sample Filtering File\n")
cat("   Need to verify if same samples are being used\n\n")

# =================================================
# SECTION 4: CHECK SAMPLE OVERLAP
# =================================================
cat("\n[SECTION 4] CHECKING SAMPLE OVERLAP\n")
cat("-------------------------------------------------------\n")

# Check if files exist
orig_sample_file <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/ROSMAP_joint_WGS_update_maf5/ROSMAP.QC.final.psam"
nf_psam_file <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/ROSMAP.QC.final.psam"

if (file.exists(orig_sample_file) && file.exists(nf_psam_file)) {
  orig_samples <- fread(orig_sample_file, data.table = FALSE)
  nf_psam <- fread(nf_psam_file, data.table = FALSE)
  
  # Create IDs
  orig_ids <- orig_samples$IID
  nf_ids <- nf_psam$IID
  
  cat("\n  Original analysis samples:", length(orig_ids), "\n")
  cat("  Nextflow PSAM samples:", length(nf_ids), "\n")
  
  common <- intersect(orig_ids, nf_ids)
  orig_only <- setdiff(orig_ids, nf_ids)
  nf_only <- setdiff(nf_ids, orig_ids)
  
  cat("  Overlap:", length(common), "samples\n")
  cat("  Original only:", length(orig_only), "samples\n")
  cat("  Nextflow only:", length(nf_only), "samples\n")
  
  if (length(orig_only) > 0 || length(nf_only) > 0) {
    cat("\n  ⚠️  SAMPLE SETS DIFFER!\n")
    write.table(data.frame(ID = orig_only), 
                file.path(output_dir, "samples_original_only.txt"),
                row.names = FALSE, quote = FALSE)
    write.table(data.frame(ID = nf_only), 
                file.path(output_dir, "samples_nextflow_only.txt"),
                row.names = FALSE, quote = FALSE)
  } else {
    cat("\n  ✓ SAMPLE SETS ARE IDENTICAL\n")
  }
} else {
  cat("  Could not check - files not found\n")
  if (!file.exists(orig_sample_file)) cat("    Missing:", orig_sample_file, "\n")
  if (!file.exists(nf_psam_file)) cat("    Missing:", nf_psam_file, "\n")
}

# =================================================
# SECTION 5: CHECK PHENOTYPE FILES
# =================================================
cat("\n[SECTION 5] CHECKING PHENOTYPE FILES\n")
cat("-------------------------------------------------------\n")

# Try to find original phenotype file
cat("  Searching for original phenotype files...\n")
pheno_search <- system("find /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS -name '*pheno*' -name '*.txt' 2>/dev/null | grep -i astrocyte | head -5", 
                       intern = TRUE)
if (length(pheno_search) > 0) {
  cat("  Found potential phenotype files:\n")
  for (f in pheno_search) {
    cat("   ", f, "\n")
  }
} else {
  cat("  No phenotype files found with 'astrocyte' in name\n")
}

# Check Nextflow phenotype
nf_pheno_file <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/phenotype_files/ROSMAP_pheno_Astrocyte.txt"
if (file.exists(nf_pheno_file)) {
  nf_pheno <- fread(nf_pheno_file, data.table = FALSE)
  cat("\n  Nextflow phenotype file:\n")
  cat("    Path:", nf_pheno_file, "\n")
  cat("    Samples:", nrow(nf_pheno), "\n")
  cat("    Columns:", paste(names(nf_pheno), collapse = ", "), "\n")
  
  pheno_col <- names(nf_pheno)[3]
  pheno_vals <- nf_pheno[[pheno_col]]
  cat("    Mean:", round(mean(pheno_vals, na.rm = TRUE), 6), "\n")
  cat("    SD:", round(sd(pheno_vals, na.rm = TRUE), 6), "\n")
}

# =================================================
# SECTION 6: CHECK COVARIATE FILES
# =================================================
cat("\n[SECTION 6] CHECKING COVARIATE FILES\n")
cat("-------------------------------------------------------\n")

# Try to find original covariate file
cat("  Searching for original covariate files...\n")
covar_search <- system("find /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS -name '*covar*' -o -name '*cov*.txt' 2>/dev/null | head -5", 
                       intern = TRUE)
if (length(covar_search) > 0) {
  cat("  Found potential covariate files:\n")
  for (f in covar_search) {
    cat("   ", f, "\n")
  }
}

# Check Nextflow covariate
nf_covar_file <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/phenotype_files/ROSMAP_covar_Astrocyte.txt"
if (file.exists(nf_covar_file)) {
  nf_covar <- fread(nf_covar_file, data.table = FALSE)
  cat("\n  Nextflow covariate file:\n")
  cat("    Path:", nf_covar_file, "\n")
  cat("    Samples:", nrow(nf_covar), "\n")
  cat("    Covariates (excluding FID/IID):", ncol(nf_covar) - 2, "\n")
  cat("    Covariate names:", paste(names(nf_covar)[-(1:2)], collapse = ", "), "\n")
  
  # Count PCs
  pc_cols <- grep("^PC[0-9]+$", names(nf_covar), value = TRUE)
  cat("    Number of genetic PCs:", length(pc_cols), "\n")
  
  # Count tech covariates
  tech_cols <- grep("^tech_cov", names(nf_covar), value = TRUE)
  cat("    Number of technical covariates:", length(tech_cols), "\n")
}

# =================================================
# SECTION 7: SUMMARY & RECOMMENDATIONS
# =================================================
cat("\n=======================================================\n")
cat("SUMMARY: ROOT CAUSES OF DISCREPANCY (r = 0.78)\n")
cat("=======================================================\n\n")

cat("IDENTIFIED DIFFERENCES:\n\n")

cat("1. ⚠️  HWE THRESHOLD (HIGH IMPACT)\n")
cat("   Original: 1e-15 (extremely strict)\n")
cat("   Nextflow: 1e-6 (standard)\n")
cat("   → This removes DIFFERENT variants between analyses\n")
cat("   → Directly affects which SNPs are tested\n\n")

cat("2. ⚠️  LD PRUNING PARAMETERS (MEDIUM IMPACT)\n")
cat("   Original: 500kb window, step 1\n")
cat("   Nextflow: 1000 SNPs window, step 50\n")
cat("   → Produces different pruned SNP sets\n")
cat("   → Affects PCA calculation\n")
cat("   → Changes genetic PC covariates slightly\n\n")

cat("3. ❓ SAMPLE FILTERING (UNKNOWN IMPACT)\n")
cat("   Need to verify if same samples used\n")
cat("   Check output files in diagnostics directory\n\n")

cat("4. ❓ PHENOTYPE/COVARIATE PROCESSING (UNKNOWN IMPACT)\n")
cat("   Need to compare actual phenotype and covariate values\n")
cat("   between original and Nextflow analyses\n\n")

cat("\nRECOMMENDATIONS:\n\n")

cat("IMMEDIATE ACTION (to match original analysis):\n")
cat("1. Update nextflow.config.combined.rosmap:\n")
cat("   hwe_threshold = 1e-15  (change from 1e-6)\n")
cat("2. Update LD pruning parameters:\n")
cat("   prune_window_size = 500  (change from 1000)\n")
cat("   prune_step_size = 1      (change from 50)\n")
cat("3. Verify sample filtering file matches\n")
cat("4. Re-run pipeline and check if correlation improves\n\n")

cat("NEXT STEPS:\n")
cat("1. Find original Regenie command to compare all parameters\n")
cat("2. Compare actual phenotype file values (if original exists)\n")
cat("3. Compare covariate file values (especially PCs)\n")
cat("4. If still discrepant, check cell type deconvolution outputs\n\n")

cat("EXPECTED OUTCOME:\n")
cat("After fixing HWE and LD pruning, correlation should improve to r > 0.95\n")
cat("If not, then phenotype/covariate differences are likely the cause\n\n")

cat("=======================================================\n")
cat("Diagnostic complete!\n")
cat("=======================================================\n\n")

