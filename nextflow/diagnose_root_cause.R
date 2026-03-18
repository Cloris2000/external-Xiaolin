#!/usr/bin/env Rscript

# Root Cause Analysis: Compare INPUTS to Regenie
# This will tell us WHY the results differ, not just THAT they differ

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("\n=================================================\n")
cat("ROOT CAUSE DIAGNOSTIC: REGENIE INPUT COMPARISON\n")
cat("=================================================\n\n")

cell_type <- "Astrocyte"
output_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/diagnostics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =================================================
# STEP 1: Locate and Compare Input Files
# =================================================
cat("[STEP 1] Identifying Input Files\n")
cat("-------------------------------------------------\n")

# Previous analysis paths (need to find these)
cat("Previous analysis inputs:\n")
cat("  Looking for phenotype file...\n")
cat("  Looking for covariate file...\n")
cat("  Looking for genotype files...\n\n")

# Nextflow pipeline paths
nf_pheno_file <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/phenotype_files/ROSMAP_pheno_Astrocyte.txt"
nf_covar_file <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/phenotype_files/ROSMAP_covar_Astrocyte.txt"
nf_psam_file <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/genotyping_qc/ROSMAP.psam"

cat("Nextflow pipeline inputs:\n")
cat("  Phenotype:", nf_pheno_file, "\n")
cat("  Covariate:", nf_covar_file, "\n")
cat("  PSAM:", nf_psam_file, "\n\n")

# Check if files exist
if (!file.exists(nf_pheno_file)) {
  cat("  ✗ Phenotype file not found!\n")
} else {
  cat("  ✓ Phenotype file exists\n")
}

if (!file.exists(nf_covar_file)) {
  cat("  ✗ Covariate file not found!\n")
} else {
  cat("  ✓ Covariate file exists\n")
}

if (!file.exists(nf_psam_file)) {
  cat("  ✗ PSAM file not found!\n")
} else {
  cat("  ✓ PSAM file exists\n")
}

# =================================================
# STEP 2: Compare Phenotype Files
# =================================================
cat("\n[STEP 2] Analyzing Nextflow Phenotype File\n")
cat("-------------------------------------------------\n")

if (file.exists(nf_pheno_file)) {
  nf_pheno <- fread(nf_pheno_file, data.table = FALSE)
  cat("  Columns:", paste(names(nf_pheno), collapse = ", "), "\n")
  cat("  Number of samples:", nrow(nf_pheno), "\n")
  cat("  Phenotype column:", names(nf_pheno)[3], "\n")
  
  # Summary statistics
  pheno_col <- names(nf_pheno)[3]
  pheno_vals <- nf_pheno[[pheno_col]]
  cat("\n  Phenotype statistics:\n")
  cat("    Mean:", round(mean(pheno_vals, na.rm = TRUE), 4), "\n")
  cat("    SD:", round(sd(pheno_vals, na.rm = TRUE), 4), "\n")
  cat("    Min:", round(min(pheno_vals, na.rm = TRUE), 4), "\n")
  cat("    Max:", round(max(pheno_vals, na.rm = TRUE), 4), "\n")
  cat("    NA values:", sum(is.na(pheno_vals)), "\n")
  
  # Check if RINT normalized (should be approximately N(0,1))
  cat("\n  Checking normalization:\n")
  cat("    Expected for RINT: Mean ≈ 0, SD ≈ 1\n")
  if (abs(mean(pheno_vals, na.rm = TRUE)) < 0.01 && abs(sd(pheno_vals, na.rm = TRUE) - 1) < 0.1) {
    cat("    ✓ Appears to be RINT normalized\n")
  } else {
    cat("    ⚠️  May not be properly normalized\n")
  }
  
  # Save sample IDs
  write.table(nf_pheno[, 1:2], 
              file.path(output_dir, "nextflow_phenotype_samples.txt"),
              row.names = FALSE, quote = FALSE)
}

# =================================================
# STEP 3: Compare Covariate Files
# =================================================
cat("\n[STEP 3] Analyzing Nextflow Covariate File\n")
cat("-------------------------------------------------\n")

if (file.exists(nf_covar_file)) {
  nf_covar <- fread(nf_covar_file, data.table = FALSE)
  cat("  Columns:", paste(names(nf_covar), collapse = ", "), "\n")
  cat("  Number of samples:", nrow(nf_covar), "\n")
  cat("  Number of covariates:", ncol(nf_covar) - 2, "\n")
  
  # List covariate names
  covar_names <- names(nf_covar)[-(1:2)]
  cat("\n  Covariate list:\n")
  for (i in seq_along(covar_names)) {
    cat("    ", i, ". ", covar_names[i], "\n", sep = "")
  }
  
  # Check for PCs
  pc_cols <- grep("^PC[0-9]+$", covar_names, value = TRUE)
  if (length(pc_cols) > 0) {
    cat("\n  Genetic PCs included:", length(pc_cols), "(", paste(pc_cols, collapse = ", "), ")\n")
  }
  
  # Check for technical covariates
  tech_cols <- grep("^tech_cov", covar_names, value = TRUE)
  if (length(tech_cols) > 0) {
    cat("  Technical covariates:", length(tech_cols), "\n")
  }
  
  # Save sample IDs
  write.table(nf_covar[, 1:2], 
              file.path(output_dir, "nextflow_covariate_samples.txt"),
              row.names = FALSE, quote = FALSE)
}

# =================================================
# STEP 4: Compare Sample Sets
# =================================================
cat("\n[STEP 4] Sample Set Comparison\n")
cat("-------------------------------------------------\n")

if (file.exists(nf_pheno_file) && file.exists(nf_covar_file) && file.exists(nf_psam_file)) {
  nf_psam <- fread(nf_psam_file, data.table = FALSE)
  
  cat("  Genotype file (PSAM):", nrow(nf_psam), "samples\n")
  cat("  Phenotype file:", nrow(nf_pheno), "samples\n")
  cat("  Covariate file:", nrow(nf_covar), "samples\n")
  
  # Check overlap
  psam_ids <- paste(nf_psam$FID, nf_psam$IID, sep = "_")
  pheno_ids <- paste(nf_pheno[[1]], nf_pheno[[2]], sep = "_")
  covar_ids <- paste(nf_covar[[1]], nf_covar[[2]], sep = "_")
  
  common_all <- Reduce(intersect, list(psam_ids, pheno_ids, covar_ids))
  cat("\n  Common samples across all files:", length(common_all), "\n")
  
  # Check mismatches
  psam_only <- setdiff(psam_ids, union(pheno_ids, covar_ids))
  pheno_only <- setdiff(pheno_ids, union(psam_ids, covar_ids))
  covar_only <- setdiff(covar_ids, union(psam_ids, pheno_ids))
  
  if (length(psam_only) > 0) {
    cat("  Samples in PSAM only:", length(psam_only), "\n")
  }
  if (length(pheno_only) > 0) {
    cat("  Samples in phenotype only:", length(pheno_only), "\n")
  }
  if (length(covar_only) > 0) {
    cat("  Samples in covariate only:", length(covar_only), "\n")
  }
  
  # This tells us the effective N for Regenie
  cat("\n  ✓ Effective sample size for Regenie:", length(common_all), "\n")
}

# =================================================
# STEP 5: Find Previous Analysis Files
# =================================================
cat("\n[STEP 5] Searching for Previous Analysis Input Files\n")
cat("-------------------------------------------------\n")

# Search common locations
search_dirs <- c(
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS",
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA",
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/GWAS",
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/phenotype"
)

cat("Searching for phenotype files containing 'Astrocyte' or 'pheno'...\n")
for (dir in search_dirs) {
  if (dir.exists(dir)) {
    files <- list.files(dir, pattern = "pheno.*txt|.*Astrocyte.*txt", 
                        recursive = TRUE, full.names = TRUE)
    if (length(files) > 0) {
      cat("\n  Found in", dir, ":\n")
      for (f in head(files, 5)) {
        cat("    ", f, "\n")
      }
    }
  }
}

cat("\nSearching for covariate files...\n")
for (dir in search_dirs) {
  if (dir.exists(dir)) {
    files <- list.files(dir, pattern = "covar.*txt|cov.*txt", 
                        recursive = TRUE, full.names = TRUE)
    if (length(files) > 0) {
      cat("\n  Found in", dir, ":\n")
      for (f in head(files, 5)) {
        cat("    ", f, "\n")
      }
    }
  }
}

# =================================================
# STEP 6: Check Regenie Command/Log Files
# =================================================
cat("\n[STEP 6] Searching for Regenie Command History\n")
cat("-------------------------------------------------\n")

# Look for .command.sh files in Nextflow work directory
nf_work_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/work_combined_ROSMAP"
if (dir.exists(nf_work_dir)) {
  cat("Searching Nextflow work directory for Regenie commands...\n")
  system(paste0("find ", nf_work_dir, 
                " -name '.command.sh' -path '*/REGENIE_STEP2*' 2>/dev/null | head -5"))
}

# Look for log files in previous analysis
prev_log_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2"
if (dir.exists(prev_log_dir)) {
  cat("\nSearching previous analysis directory for logs...\n")
  system(paste0("ls -lh ", prev_log_dir, "/*.log 2>/dev/null | head -5"))
  system(paste0("ls -lh ", prev_log_dir, "/*.sh 2>/dev/null | head -5"))
}

# =================================================
# STEP 7: Recommendations
# =================================================
cat("\n=================================================\n")
cat("RECOMMENDATIONS FOR NEXT STEPS\n")
cat("=================================================\n\n")

cat("To identify the root cause, you need to:\n\n")

cat("1. FIND THE PREVIOUS ANALYSIS INPUT FILES:\n")
cat("   Look for the exact phenotype and covariate files used in:\n")
cat("   /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2/\n\n")

cat("2. CHECK THE REGENIE COMMAND USED PREVIOUSLY:\n")
cat("   Look for .sh or .log files containing the Regenie command\n")
cat("   Key parameters to check:\n")
cat("     - Which covariates were included (--covarColList)\n")
cat("     - Phenotype column name\n")
cat("     - Sample filtering options\n")
cat("     - Step 1 prediction file used\n\n")

cat("3. COMPARE PHENOTYPE VALUES:\n")
cat("   Once you find the previous phenotype file:\n")
cat("   - Check if normalization method is the same (RINT)\n")
cat("   - Compare mean and SD\n")
cat("   - Check sample overlap\n\n")

cat("4. COMPARE COVARIATE VALUES:\n")
cat("   Once you find the previous covariate file:\n")
cat("   - Count number of PCs included\n")
cat("   - Check if technical covariates match\n")
cat("   - Verify PC values are identical\n\n")

cat("5. COMPARE GENOTYPE FILES:\n")
cat("   - Check if same QC filters were applied\n")
cat("   - Verify MAF thresholds\n")
cat("   - Check imputation quality filters\n\n")

cat("LIKELY CAUSES (based on r = 0.78):\n")
cat("  • Different number or values of covariates (most common)\n")
cat("  • Different phenotype normalization\n")
cat("  • Different sample sets\n")
cat("  • Different genotype QC filters\n\n")

cat("NEXT ACTION:\n")
cat("  Please provide the location of your previous analysis\n")
cat("  phenotype and covariate files, and I can do a direct\n")
cat("  comparison.\n\n")

cat("Output files saved to:", output_dir, "\n")
cat("  - nextflow_phenotype_samples.txt\n")
cat("  - nextflow_covariate_samples.txt\n\n")

cat("=================================================\n\n")

