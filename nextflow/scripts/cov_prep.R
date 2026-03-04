#!/usr/bin/env Rscript
# Covariate Preparation for GWAS
# Merges PCA with clinical covariates
# This step comes AFTER genotype QC and PCA calculation

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--pca_file"), type="character", default=NULL,
              help="Path to PCA file (from genotype QC)"),
  make_option(c("--clinical_cov_file"), type="character", default=NULL,
              help="Path to clinical covariates file (from pheno_prep)"),
  make_option(c("--samples_file"), type="character", default=NULL,
              help="Path to samples with phenotypes file (for validation)"),
  make_option(c("--output_file"), type="character", default="covariates.txt",
              help="Output filename for final covariate file"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("\n=============================================================================\n")
cat("COVARIATE PREPARATION (Step 4 of Pipeline)\n")
cat("=============================================================================\n\n")
cat("This step merges:\n")
cat("  1. PCA (from genotype QC)\n")
cat("  2. Clinical covariates (sex, age from pheno_prep)\n\n")
cat("NOTE: This step requires PCA to be calculated AFTER genotype QC!\n")
cat("=============================================================================\n\n")

# Normalize duplicated IDs (X_X -> X) for matching (e.g. NABEC pipeline vs individual scripts)
norm_id_vec <- function(v) {
  s <- as.character(v)
  sub("^(.+)_\\1$", "\\1", s)
}

# Load files
cat("Loading PCA file:", opt$pca_file, "\n")
pca <- read.csv(opt$pca_file)
cat("  PCA samples:", nrow(pca), "\n")
cat("  PCA columns:", paste(colnames(pca), collapse=", "), "\n\n")

# PCA may use FID or #FID; normalize duplicated IDs for merge
fid_col_pca <- if ("FID" %in% colnames(pca)) "FID" else if ("X.FID" %in% colnames(pca)) "X.FID" else names(pca)[1]
if (fid_col_pca %in% colnames(pca)) {
  pca[[fid_col_pca]] <- norm_id_vec(pca[[fid_col_pca]])
}
if ("IID" %in% colnames(pca)) pca$IID <- norm_id_vec(pca$IID)

cat("Loading clinical covariates:", opt$clinical_cov_file, "\n")
clinical_cov <- read.table(opt$clinical_cov_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
cat("  Clinical covariate samples:", nrow(clinical_cov), "\n")
cat("  Clinical covariate columns:", paste(colnames(clinical_cov), collapse=", "), "\n\n")

# Normalize genotyping_FID for merge (handles X_X format from pipeline)
if ("genotyping_FID" %in% colnames(clinical_cov)) {
  clinical_cov$genotyping_FID <- norm_id_vec(clinical_cov$genotyping_FID)
}

cat("Loading samples file:", opt$samples_file, "\n")
samples_with_pheno <- read.table(opt$samples_file, header=FALSE, skip=1, sep="\t", stringsAsFactors=FALSE)
colnames(samples_with_pheno) <- c("FID", "IID")
samples_with_pheno$FID <- norm_id_vec(samples_with_pheno$FID)
samples_with_pheno$IID <- norm_id_vec(samples_with_pheno$IID)
cat("  Expected samples:", nrow(samples_with_pheno), "\n\n")

# Check if clinical_cov has genotyping_FID mapping
if ("genotyping_FID" %in% colnames(clinical_cov)) {
  cat("Found genotyping_FID mapping column\n")
  cat("  This indicates some samples had FID remapping (clinical ID -> genotyping ID)\n\n")
  
  # PCA uses genotyping FIDs, clinical_cov has clinical FIDs
  # We need to map: PCA (genotyping FID) -> clinical_cov (using genotyping_FID column)
  
  # Merge PCA with clinical covariates using genotyping_FID
  pca_cov_merged <- merge(pca, 
                          clinical_cov[, c("genotyping_FID", "msex", "age_death", "age_death_sex",
                                          "age_death2", "age_death2_sex")],
                          by.x = fid_col_pca,
                          by.y = "genotyping_FID",
                          all.x = FALSE, all.y = FALSE)
  # Ensure FID column exists for downstream (rename if needed)
  if (fid_col_pca != "FID" && fid_col_pca %in% colnames(pca_cov_merged)) {
    colnames(pca_cov_merged)[colnames(pca_cov_merged) == fid_col_pca] <- "FID"
  }
  
  cat("Merged PCA + clinical covariates:\n")
  cat("  Total samples:", nrow(pca_cov_merged), "\n")
  
} else {
  cat("No genotyping_FID mapping - all samples matched directly\n")
  cat("  Merging PCA with clinical covariates by FID\n\n")
  
  # Simple merge by FID
  pca_cov_merged <- merge(pca,
                          clinical_cov[, c("FID", "msex", "age_death", "age_death_sex",
                                          "age_death2", "age_death2_sex")],
                          by = "FID",
                          all.x = FALSE, all.y = FALSE)
  
  cat("Merged PCA + clinical covariates:\n")
  cat("  Total samples:", nrow(pca_cov_merged), "\n")
}

# Check for duplicate column names from merge
if (any(grepl("\\.x$|\\.y$", colnames(pca_cov_merged)))) {
  cat("\nWarning: Found duplicate columns from merge. Cleaning up...\n")
  
  # Keep .x columns (from PCA) and rename them
  dup_cols <- unique(gsub("\\.x$|\\.y$", "", grep("\\.x$|\\.y$", colnames(pca_cov_merged), value=TRUE)))
  for (col in dup_cols) {
    if (paste0(col, ".x") %in% colnames(pca_cov_merged)) {
      pca_cov_merged[[col]] <- pca_cov_merged[[paste0(col, ".x")]]
    }
  }
  
  # Remove .x and .y columns
  pca_cov_merged <- pca_cov_merged[, !grepl("\\.x$|\\.y$", colnames(pca_cov_merged))]
}

# Reorder columns: FID, IID, PC1-PC10, msex, age_death, age_death_sex, age_death2, age_death2_sex
pc_cols <- grep("^PC[0-9]+$", colnames(pca_cov_merged), value = TRUE)
other_cols <- c("msex", "age_death", "age_death_sex", "age_death2", "age_death2_sex")
other_cols <- other_cols[other_cols %in% colnames(pca_cov_merged)]

pca_cov_merged <- pca_cov_merged[, c("FID", "IID", pc_cols, other_cols)]

cat("  Final columns:", paste(colnames(pca_cov_merged), collapse=", "), "\n\n")

# Validation
cat("Validation:\n")
cat("  Samples in covariates:", nrow(pca_cov_merged), "\n")
cat("  Samples expected (with phenotypes):", nrow(samples_with_pheno), "\n")

overlap <- sum(pca_cov_merged$FID %in% samples_with_pheno$FID)
cat("  Overlap:", overlap, "\n")

if (overlap < nrow(pca_cov_merged) * 0.95) {
  cat("\n⚠ Warning: Less than 95% of covariate samples have phenotypes!\n")
  cat("  This may indicate a mismatch between phenotype and genotype data.\n")
} else {
  cat("\n✓ Covariate samples match phenotype samples\n")
}

# Save covariate file
# Always write to current directory - Nextflow will publish to final location
output_path <- opt$output_file

cat("\nSaving covariate file:", output_path, "\n")
write.table(pca_cov_merged,
            output_path,
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\n=============================================================================\n")
cat("COVARIATE PREPARATION COMPLETED!\n")
cat("=============================================================================\n")
cat("Output file:", output_path, "\n")
cat("Total samples:", nrow(pca_cov_merged), "\n")
cat("Covariates:", paste(setdiff(colnames(pca_cov_merged), c("FID", "IID")), collapse=", "), "\n\n")
cat("Next step: Run GWAS with phenotypes and covariates\n")
cat("=============================================================================\n\n")

