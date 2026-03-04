#!/usr/bin/env Rscript
# Phenotype Preparation for GWAS
# Creates phenotypes, sample list, and clinical covariates (NO PCA NEEDED)
# This step comes FIRST in the pipeline, before genotype QC

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyverse)
  library(broom)
  library(RNOmni)
})

# Helper function to find column names with flexible matching
find_column <- function(df, possible_names, required = FALSE) {
  for (name in possible_names) {
    if (name %in% colnames(df)) {
      return(name)
    }
  }
  if (required) {
    stop(paste("Required column not found. Tried:", paste(possible_names, collapse = ", ")))
  }
  return(NULL)
}

# Parse command line arguments
option_list <- list(
  make_option(c("--cell_proportions"), type="character", default=NULL,
              help="Path to cell proportions CSV file"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Path to cleaned metadata CSV file"),
  make_option(c("--biospecimen_file"), type="character", default=NULL,
              help="Path to biospecimen metadata file (optional)"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--study"), type="character", default=NULL,
              help="Study name (e.g., ROSMAP, CMC, etc.)"),
  make_option(c("--phenotype_output"), type="character", default="phenotypes_RINT.txt",
              help="Output filename for phenotype file"),
  make_option(c("--samples_output"), type="character", default="samples_with_phenotypes.txt",
              help="Output filename for samples file"),
  make_option(c("--clinical_cov_output"), type="character", default="clinical_covariates.txt",
              help="Output filename for clinical covariates (sex, age)"),
  make_option(c("--col_msex"), type="character", default=NULL,
              help="Column name for sex (will auto-detect if not specified)"),
  make_option(c("--col_age"), type="character", default=NULL,
              help="Column name for age_death (will auto-detect if not specified)"),
  make_option(c("--col_individualID"), type="character", default=NULL,
              help="Column name for individualID (will auto-detect if not specified)"),
  make_option(c("--col_study"), type="character", default=NULL,
              help="Column name for Study (will auto-detect if not specified)"),
  make_option(c("--col_projid"), type="character", default=NULL,
              help="Column name for projid (will auto-detect if not specified)"),
  make_option(c("--col_specimenID"), type="character", default=NULL,
              help="Column name for specimenID (will auto-detect if not specified)"),
  make_option(c("--fid_method"), type="character", default="study_projid",
              help="Method for creating FID: study_projid, individual_id, specimen_id, etc."),
  make_option(c("--fid_format"), type="character", default=NULL,
              help="Format string for custom FID creation"),
  make_option(c("--biospec_col_individual"), type="character", default="individualID",
              help="Column name in biospecimen file for individual ID"),
  make_option(c("--biospec_col_specimen"), type="character", default="specimenID",
              help="Column name in biospecimen file for specimen ID")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("\n=============================================================================\n")
cat("PHENOTYPE PREPARATION (Step 1 of Pipeline)\n")
cat("=============================================================================\n\n")
cat("This step creates:\n")
cat("  1. Phenotypes (RINT-transformed cell proportions)\n")
cat("  2. Sample list (samples with phenotype data)\n")
cat("  3. Clinical covariates (sex, age, age interactions)\n\n")
cat("NOTE: This step does NOT require genotype data or PCA!\n")
cat("=============================================================================\n\n")

cat("Loading cell proportions and metadata...\n")
# Load cell proportions
ROSMAP_estimations_scaled <- read.csv(opt$cell_proportions)

# Load metadata - handle both RData and CSV formats
if (endsWith(opt$metadata, ".RData") || endsWith(opt$metadata, ".rds")) {
  ROSMAP_meta_temp_cleaned <- readRDS(opt$metadata)
  if (inherits(ROSMAP_meta_temp_cleaned, "tbl_df") || inherits(ROSMAP_meta_temp_cleaned, "tbl")) {
    ROSMAP_meta_temp_cleaned <- as.data.frame(ROSMAP_meta_temp_cleaned)
  }
} else {
  ROSMAP_meta_temp_cleaned <- read.csv(opt$metadata)
}

# Remove excluded samples
if ("exclude" %in% colnames(ROSMAP_meta_temp_cleaned)) {
  ROSMAP_meta_temp_cleaned$exclude[is.na(ROSMAP_meta_temp_cleaned$exclude)] <- FALSE
  if (!is.logical(ROSMAP_meta_temp_cleaned$exclude)) {
    ROSMAP_meta_temp_cleaned$exclude <- as.logical(ROSMAP_meta_temp_cleaned$exclude)
  }
  ROSMAP_meta_temp_cleaned <- ROSMAP_meta_temp_cleaned[!ROSMAP_meta_temp_cleaned$exclude, ]
}

cat("Number of samples in metadata:", nrow(ROSMAP_meta_temp_cleaned), "\n")
cat("DEBUG: Column names in metadata:", paste(head(colnames(ROSMAP_meta_temp_cleaned), 10), collapse=", "), "\n")

# Find column names using flexible matching
col_msex <- if (!is.null(opt$col_msex)) opt$col_msex else find_column(ROSMAP_meta_temp_cleaned, c("msex", "sex", "gender", "Sex", "Gender", "SEX.y", "SEX.x", "reportedSex", "reportedGender", "biological_sex", "phenotypic_sex", "GENDER", "SEX"), required = TRUE)
col_age <- if (!is.null(opt$col_age)) opt$col_age else find_column(ROSMAP_meta_temp_cleaned, c("age_death", "age_at_death", "age", "Age", "AGE.y", "AGE.x", "ageDeath", "AGE"), required = TRUE)
col_individualID <- if (!is.null(opt$col_individualID)) opt$col_individualID else find_column(ROSMAP_meta_temp_cleaned, c("individualID", "individual_id", "subjectID", "subject_id", "Individual"), required = TRUE)
col_study <- if (!is.null(opt$col_study)) opt$col_study else find_column(ROSMAP_meta_temp_cleaned, c("Study", "study", "cohort", "Cohort"))
col_projid <- if (!is.null(opt$col_projid)) opt$col_projid else find_column(ROSMAP_meta_temp_cleaned, c("projid", "proj_id", "project_id"))
col_specimenID <- if (!is.null(opt$col_specimenID)) opt$col_specimenID else find_column(ROSMAP_meta_temp_cleaned, c("specimenID", "specimen_id", "sampleID", "sample_id", "SAMPID"))
col_synapseID <- find_column(ROSMAP_meta_temp_cleaned, c("synapseID", "synapse_id", "SynapseID", "SYNAPSE_ID"), required = FALSE)

cat("Using columns:\n")
cat("  Sex:", col_msex, "\n")
cat("  Age:", col_age, "\n")
cat("  IndividualID:", col_individualID, "\n")
if (!is.null(col_study)) cat("  Study:", col_study, "\n")
if (!is.null(col_projid)) cat("  ProjID:", col_projid, "\n")
if (!is.null(col_specimenID)) cat("  SpecimenID:", col_specimenID, "\n")
if (!is.null(col_synapseID)) cat("  SynapseID:", col_synapseID, "\n")

# Prepare metadata subset with required columns
ROSMAP_meta_sub <- ROSMAP_meta_temp_cleaned[, c(col_msex, col_age, col_individualID)]
colnames(ROSMAP_meta_sub) <- c("msex", "age_death", "individualID")

# Convert sex to numeric (0/1) for mathematical operations
# This is CRITICAL for creating age*sex interaction terms
# NABEC/phs001300 use 1=male, 2=female; map both 0/1 and 1/2 codings
ROSMAP_meta_sub$msex <- tolower(trimws(as.character(ROSMAP_meta_sub$msex)))
ROSMAP_meta_sub$msex <- ifelse(ROSMAP_meta_sub$msex %in% c("male", "m", "1"), 1,
                               ifelse(ROSMAP_meta_sub$msex %in% c("female", "f", "0", "2"), 0, NA))
cat("  - Converted msex to numeric (0=female, 1=male)\n")
cat("  - Samples with valid sex:", sum(!is.na(ROSMAP_meta_sub$msex)), "\n")

# Add Study and projid if available
if (!is.null(col_study) && col_study %in% colnames(ROSMAP_meta_temp_cleaned)) {
  ROSMAP_meta_sub$Study <- ROSMAP_meta_temp_cleaned[[col_study]]
}
if (!is.null(col_projid) && col_projid %in% colnames(ROSMAP_meta_temp_cleaned)) {
  ROSMAP_meta_sub$projid <- ROSMAP_meta_temp_cleaned[[col_projid]]
}

# Create FID based on fid_method
if (opt$fid_method == "study_projid") {
  if (!"Study" %in% colnames(ROSMAP_meta_sub) || !"projid" %in% colnames(ROSMAP_meta_sub)) {
    stop("fid_method='study_projid' requires Study and projid columns")
  }
  ROSMAP_meta_sub$FID <- paste0(ROSMAP_meta_sub$Study, ROSMAP_meta_sub$projid)
} else if (opt$fid_method == "individual_id") {
  ROSMAP_meta_sub$FID <- ROSMAP_meta_sub$individualID
} else {
  stop(paste("Unknown fid_method:", opt$fid_method))
}

# Clean age_death
# Handle age ranges (e.g., "50-59" → 54.5)
ROSMAP_meta_sub$age_death <- gsub("90\\+", "90", ROSMAP_meta_sub$age_death)
if (any(grepl("-", ROSMAP_meta_sub$age_death, fixed = TRUE))) {
  cat("  - Detected age ranges, converting to midpoints...\n")
  ROSMAP_meta_sub$age_death <- sapply(ROSMAP_meta_sub$age_death, function(x) {
    if (grepl("-", x, fixed = TRUE)) {
      parts <- as.numeric(strsplit(x, "-")[[1]])
      if (length(parts) == 2) return(mean(parts))
    }
    return(as.numeric(x))
  })
} else {
  ROSMAP_meta_sub$age_death <- as.numeric(gsub("[^0-9.]", NA, ROSMAP_meta_sub$age_death))
}

# Create age interaction terms
ROSMAP_meta_sub$age_death2 <- ROSMAP_meta_sub$age_death^2
ROSMAP_meta_sub$age_death_sex <- ROSMAP_meta_sub$age_death * ROSMAP_meta_sub$msex
ROSMAP_meta_sub$age_death2_sex <- ROSMAP_meta_sub$age_death2 * ROSMAP_meta_sub$msex

# Remove rows with any NA values
ROSMAP_meta_sub_clean <- ROSMAP_meta_sub[complete.cases(ROSMAP_meta_sub), ]

# Keep only the first occurrence of each FID
ROSMAP_meta_sub_unique <- ROSMAP_meta_sub_clean[!duplicated(ROSMAP_meta_sub_clean$FID), ]
rownames(ROSMAP_meta_sub_unique) <- ROSMAP_meta_sub_unique$FID

cat("Samples with complete clinical data:", nrow(ROSMAP_meta_sub_unique), "\n")

# Remove columns with zero variance
num_unique_values <- sapply(ROSMAP_meta_temp_cleaned, function(x) length(unique(x)))
cols_to_keep <- num_unique_values > 1
ROSMAP_meta_temp_cleaned_rmVar <- ROSMAP_meta_temp_cleaned[, cols_to_keep]

# Deduplicate using synapseID (priority), specimenID, or individualID
if (!is.null(col_synapseID) && col_synapseID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  # Use synapseID for cohorts with Synapse data (ROSMAP, Mayo, MSBB)
  ROSMAP_meta_temp_cleaned_rmVar <- ROSMAP_meta_temp_cleaned_rmVar[!duplicated(ROSMAP_meta_temp_cleaned_rmVar[[col_synapseID]]), ]
  rownames(ROSMAP_meta_temp_cleaned_rmVar) <- ROSMAP_meta_temp_cleaned_rmVar[[col_synapseID]]
} else if (!is.null(col_specimenID) && col_specimenID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  ROSMAP_meta_temp_cleaned_rmVar <- ROSMAP_meta_temp_cleaned_rmVar[!duplicated(ROSMAP_meta_temp_cleaned_rmVar[[col_specimenID]]), ]
  rownames(ROSMAP_meta_temp_cleaned_rmVar) <- ROSMAP_meta_temp_cleaned_rmVar[[col_specimenID]]
} else {
  # Use individualID as fallback
  ROSMAP_meta_temp_cleaned_rmVar <- ROSMAP_meta_temp_cleaned_rmVar[!duplicated(ROSMAP_meta_temp_cleaned_rmVar[[col_individualID]]), ]
  rownames(ROSMAP_meta_temp_cleaned_rmVar) <- ROSMAP_meta_temp_cleaned_rmVar[[col_individualID]]
}

# Prepare estimations
# Build column list based on what's available (priority: synapseID > specimenID > individualID)
matcher_cols <- c()
if (!is.null(col_synapseID) && col_synapseID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  matcher_cols <- c(matcher_cols, col_synapseID)
} else if (!is.null(col_specimenID) && col_specimenID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  matcher_cols <- c(matcher_cols, col_specimenID)
}
matcher_cols <- c(matcher_cols, col_individualID)

ROSMAP_meta_matcher <- ROSMAP_meta_temp_cleaned_rmVar[, matcher_cols, drop = FALSE]
if (!is.null(col_study)) ROSMAP_meta_matcher$Study <- ROSMAP_meta_temp_cleaned_rmVar[[col_study]]
if (!is.null(col_projid)) ROSMAP_meta_matcher$projid <- ROSMAP_meta_temp_cleaned_rmVar[[col_projid]]

# Rename columns appropriately based on what ID type was found
if (!is.null(col_synapseID) && col_synapseID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  # For Synapse cohorts: synapseID is the merge key
  if (length(matcher_cols) == 2) {
    colnames(ROSMAP_meta_matcher)[1:2] <- c("specimenID", "individualID")  # synapseID → specimenID for merging
  } else {
    colnames(ROSMAP_meta_matcher)[1] <- "specimenID"
    ROSMAP_meta_matcher$individualID <- ROSMAP_meta_matcher$specimenID
  }
} else if (!is.null(col_specimenID) && col_specimenID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  colnames(ROSMAP_meta_matcher)[1:2] <- c("specimenID", "individualID")
} else {
  colnames(ROSMAP_meta_matcher)[1] <- "individualID"
  ROSMAP_meta_matcher$specimenID <- ROSMAP_meta_matcher$individualID  # Use individualID as specimenID
}

# Find matching column in cell proportions
cell_id_col <- if (!is.null(col_specimenID) && col_specimenID %in% colnames(ROSMAP_estimations_scaled)) {
  col_specimenID
} else if ("specimenID" %in% colnames(ROSMAP_estimations_scaled)) {
  "specimenID"
} else if (col_individualID %in% colnames(ROSMAP_estimations_scaled)) {
  col_individualID
} else {
  colnames(ROSMAP_estimations_scaled)[1]
}

cat("DEBUG: Merging cell proportions with metadata...\n")
cat("  - Metadata specimenID sample (first 3):", head(ROSMAP_meta_matcher$specimenID, 3), "\n")
cat("  - Cell proportion ID column:", cell_id_col, "\n")
cat("  - Cell proportion ID sample (first 3):", head(ROSMAP_estimations_scaled[[cell_id_col]], 3), "\n")

ROSMAP_estimations_combined <- merge(ROSMAP_meta_matcher, ROSMAP_estimations_scaled,
                                     by.x = "specimenID", by.y = cell_id_col)

# Create FID for estimations
if (opt$fid_method == "study_projid" && "Study" %in% colnames(ROSMAP_estimations_combined)) {
  ROSMAP_estimations_combined$FID <- paste0(ROSMAP_estimations_combined$Study, ROSMAP_estimations_combined$projid)
} else {
  ROSMAP_estimations_combined$FID <- ROSMAP_estimations_combined$individualID
}
ROSMAP_estimations_combined$IID <- ROSMAP_estimations_combined$FID

# Select cell type columns
exclude_cols <- c("specimenID", "individualID", "Study", "projid", "FID", "IID")
cell_cols <- setdiff(colnames(ROSMAP_estimations_combined), exclude_cols)
ROSMAP_estimations_named <- ROSMAP_estimations_combined[, c("FID", "IID", cell_cols)]

# Merge with metadata for RINT transformation
combined_df <- merge(ROSMAP_meta_sub_unique, ROSMAP_estimations_named, by.x = "FID", by.y = "FID")
combined_df <- combined_df[!duplicated(combined_df$FID), ]

cat("Samples with both phenotype and clinical data:", nrow(combined_df), "\n")

# Perform linear regression to regress out sex and age, then apply RINT
cell_types <- cell_cols
pheno_lms <- lapply(cell_types, function(cell_type) {
  lm_formula <- paste0("scale(", cell_type, ") ~ msex + scale(age_death) + scale(age_death_sex) + scale(age_death2) + scale(age_death2_sex)")
  results <- lm(lm_formula, data = combined_df) %>%
    residuals()
  results <- results[!is.na(results)] %>%
    RNOmni::RankNorm() %>%
    as.data.frame()
  results$cell_type <- cell_type
  results$FID <- combined_df$FID[!is.na(results$.)]
  return(results)
}) %>% bind_rows()

# Format phenotypes
names(pheno_lms)[1] <- "transformed_residuals"
pivot_df <- pivot_wider(data = pheno_lms,
                        names_from = "cell_type",
                        values_from = "transformed_residuals")
pivot_df$IID <- pivot_df$FID
pivot_df <- pivot_df[, c(ncol(pivot_df), 1:(ncol(pivot_df) - 1))]
colnames(pivot_df)[1] <- "FID"
colnames(pivot_df)[2] <- "IID"

cat("Total samples with phenotypes:", nrow(pivot_df), "\n")

# Use ALL samples with bulk RNA-seq (no WGS/genotype filtering at this stage)
# Genotype QC downstream will filter to samples present in the VCF; REGENIE uses intersection
samples_final <- pivot_df
cat("\nIncluding all samples with bulk RNA-seq (no genotype overlap filter):", nrow(samples_final), "\n")

# Save phenotype file
cat("\nSaving phenotype file:", opt$phenotype_output, "\n")
write.table(samples_final,
            opt$phenotype_output,
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save samples file (for genotype QC filtering)
cat("Saving samples file:", opt$samples_output, "\n")
samples_with_phenotypes <- samples_final[, c("FID", "IID"), drop = FALSE]
samples_with_phenotypes <- rbind(c("FID", "IID"), samples_with_phenotypes)
write.table(samples_with_phenotypes,
            opt$samples_output,
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Save clinical covariates (for later merging with PCA)
cat("Saving clinical covariates:", opt$clinical_cov_output, "\n")

# Clinical covariates for all samples (genotyping_FID = FID; no ID mapping)
clinical_cov_df <- combined_df[combined_df$FID %in% samples_final$FID,
  c("FID", "IID", "msex", "age_death", "age_death_sex", "age_death2", "age_death2_sex"), drop = FALSE]
clinical_cov_df$genotyping_FID <- clinical_cov_df$FID

write.table(clinical_cov_df,
            opt$clinical_cov_output,
            row.names = FALSE, quote = FALSE, sep = "\t")

# Copy to output_dir if different
if (opt$output_dir != "." && opt$output_dir != getwd()) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(opt$phenotype_output, file.path(opt$output_dir, opt$phenotype_output), overwrite = TRUE)
  file.copy(opt$samples_output, file.path(opt$output_dir, opt$samples_output), overwrite = TRUE)
  file.copy(opt$clinical_cov_output, file.path(opt$output_dir, opt$clinical_cov_output), overwrite = TRUE)
}

cat("\n=============================================================================\n")
cat("PHENOTYPE PREPARATION COMPLETED!\n")
cat("=============================================================================\n")
cat("Output files:\n")
cat("  1. Phenotypes:", opt$phenotype_output, "(", nrow(samples_final), "samples )\n")
cat("  2. Sample list:", opt$samples_output, "\n")
cat("  3. Clinical covariates:", opt$clinical_cov_output, "\n\n")
cat("Next step: Run genotype QC using the sample list to filter to phenotyped samples\n")
cat("=============================================================================\n\n")

