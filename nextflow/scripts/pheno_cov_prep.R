#!/usr/bin/env Rscript
# Phenotype and Covariate File Preparation for GWAS
# Prepares final phenotype and covariate files with RINT transformation

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyverse)
  library(broom)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(ggsignif)
  library(modelr)
  library(cowplot)
  library(gridExtra)
  library(RColorBrewer)
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
  make_option(c("--wgs_psam_file"), type="character", default=NULL,
              help="Path to WGS PSAM file"),
  make_option(c("--pca_file"), type="character", default=NULL,
              help="Path to PCA file"),
  make_option(c("--biospecimen_file"), type="character", default=NULL,
              help="Path to biospecimen metadata file"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--study"), type="character", default=NULL,
              help="Study name (e.g., ROSMAP, CMC, etc.) - used for cohort-specific ID mapping logic"),
  make_option(c("--phenotype_output"), type="character", default="phenotypes_RINT.txt",
              help="Output filename for phenotype file"),
  make_option(c("--covariate_output"), type="character", default="covariates.txt",
              help="Output filename for covariate file"),
  make_option(c("--samples_output"), type="character", default="samples_with_phenotypes.txt",
              help="Output filename for samples file"),
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
              help="Method for creating FID: study_projid, individual_id, specimen_id, submitted_subject_id, genotyping_sample_id, subject_id_from_sampid"),
  make_option(c("--fid_format"), type="character", default=NULL,
              help="Format string for custom FID creation (e.g., '{study}_{projid}'). Only used if fid_method='custom'"),
  make_option(c("--biospec_col_individual"), type="character", default="individualID",
              help="Column name in biospecimen file for individual/subject ID (default: individualID)"),
  make_option(c("--biospec_col_specimen"), type="character", default="specimenID",
              help="Column name in biospecimen file for specimen/genotyping ID (default: specimenID)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Loading cell proportions and metadata...\n")
# Load cell proportions
ROSMAP_estimations_scaled <- read.csv(opt$cell_proportions)

# Load metadata - handle both RData and CSV formats
if (endsWith(opt$metadata, ".RData") || endsWith(opt$metadata, ".rds")) {
  ROSMAP_meta_temp_cleaned <- readRDS(opt$metadata)
  # Convert to data.frame if it's a tibble
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

cat("Preparing metadata for covariate file...\n")

# Detect column names
# For GTEx: after joining with subject metadata, sex/age columns are SEX.y and AGE.y
col_msex <- if (!is.null(opt$col_msex)) {
  opt$col_msex
} else {
  find_column(ROSMAP_meta_temp_cleaned, c("msex", "Msex", "sex", "Sex", "SEX.y", "SEX", "gender", "Gender"), required = TRUE)
}

col_age <- if (!is.null(opt$col_age)) {
  opt$col_age
} else {
  find_column(ROSMAP_meta_temp_cleaned, c("age_death", "ageDeath", "age_at_death", "Age_death", "AgeDeath", "AGE.y", "AGE", "age"), required = TRUE)
}

col_individualID <- if (!is.null(opt$col_individualID)) {
  opt$col_individualID
} else {
  find_column(ROSMAP_meta_temp_cleaned, c("individualID", "IndividualID", "individual_id", "subjectID", "SubjectID"))
}

cat("DEBUG: opt$col_study =", if(is.null(opt$col_study)) "NULL" else opt$col_study, "\n")
col_study <- if (!is.null(opt$col_study)) {
  cat("DEBUG: Using opt$col_study:", opt$col_study, "\n")
  opt$col_study
} else {
  result <- find_column(ROSMAP_meta_temp_cleaned, c("Study", "study", "cohort", "Cohort"))
  cat("DEBUG: find_column returned:", if(is.null(result)) "NULL" else result, "\n")
  result
}

col_projid <- if (!is.null(opt$col_projid)) {
  opt$col_projid
} else {
  find_column(ROSMAP_meta_temp_cleaned, c("projid", "Projid", "proj_id", "subject_id", "SubjectID"))
}

cat("Detected column names:\n")
cat("  sex:", col_msex, "\n")
cat("  age:", col_age, "\n")
cat("  individualID:", ifelse(is.null(col_individualID), "NOT FOUND", col_individualID), "\n")
cat("  study (col_study):", ifelse(is.null(col_study), "NOT FOUND", col_study), "\n")
cat("  projid:", ifelse(is.null(col_projid), "NOT FOUND", col_projid), "\n")
cat("DEBUG: Columns in ROSMAP_meta_temp_cleaned:", paste(colnames(ROSMAP_meta_temp_cleaned), collapse=", "), "\n")
cat("DEBUG: Is 'Study' in columns?", "Study" %in% colnames(ROSMAP_meta_temp_cleaned), "\n")

# Extract metadata for covariates
required_cols <- c(col_msex, col_age)
new_colnames <- c("msex", "age_death")

if (!is.null(col_individualID)) {
  required_cols <- c(required_cols, col_individualID)
  new_colnames <- c(new_colnames, "individualID")
}
if (!is.null(col_study)) {
  required_cols <- c(required_cols, col_study)
  new_colnames <- c(new_colnames, "Study")
}
if (!is.null(col_projid)) {
  required_cols <- c(required_cols, col_projid)
  new_colnames <- c(new_colnames, "projid")
}

ROSMAP_meta_sub <- ROSMAP_meta_temp_cleaned[, required_cols, drop = FALSE]
colnames(ROSMAP_meta_sub) <- new_colnames
cat("  Renamed columns:", paste(paste(required_cols, "->", new_colnames), collapse=", "), "\n")

# Create FID for metadata using parameterized method (same as later in script)
# This FID is for matching clinical metadata with cell proportions (both use same method)
# Move create_fid function definition earlier so we can use it here
# For now, use simple logic that matches the parameterized method
cat("Creating FID for metadata using method:", opt$fid_method, "\n")
if (opt$fid_method == "study_projid") {
  if (!("Study" %in% colnames(ROSMAP_meta_sub) && "projid" %in% colnames(ROSMAP_meta_sub))) {
    stop("ERROR: fid_method='study_projid' requires Study and projid columns in metadata")
  }
  ROSMAP_meta_sub$FID <- paste0(ROSMAP_meta_sub$Study, ROSMAP_meta_sub$projid)
} else if (opt$fid_method == "individual_id") {
  if (!("individualID" %in% colnames(ROSMAP_meta_sub))) {
    stop("ERROR: fid_method='individual_id' requires individualID column in metadata")
  }
  ROSMAP_meta_sub$FID <- ROSMAP_meta_sub$individualID
} else if (opt$fid_method == "specimen_id") {
  # For specimen_id, we need to find the specimenID column in the original metadata
  # But ROSMAP_meta_sub might not have it, so we'll handle this in the later merge
  # For now, try to use individualID if available
  if ("individualID" %in% colnames(ROSMAP_meta_sub)) {
    ROSMAP_meta_sub$FID <- ROSMAP_meta_sub$individualID
  } else {
    stop("ERROR: fid_method='specimen_id' requires individualID or specimenID column")
  }
} else if (opt$fid_method == "subject_id_from_sampid") {
  # GTEx: Extract subject ID from SAMPID (first two parts separated by "-")
  # Need to find SAMPID in the original metadata before subsetting
  sampid_col <- NULL
  for (col in c("SAMPID", "sampid", "specimenID", "SpecimenID", "specimen_id", 
                "synapseID", "SynapseID", "sampleID", "SampleID")) {
    if (col %in% colnames(ROSMAP_meta_temp_cleaned)) {
      sampid_col <- col
      break
    }
  }
  if (is.null(sampid_col)) {
    stop("ERROR: fid_method='subject_id_from_sampid' requires SAMPID or specimenID column in metadata")
  }
  # Extract subject ID: split by "-" and take first two parts
  fid_from_meta <- sapply(strsplit(as.character(ROSMAP_meta_temp_cleaned[[sampid_col]]), "-"), function(x) {
    if (length(x) >= 2) {
      paste(x[1], x[2], sep = "-")
    } else {
      NA_character_
    }
  })
  # Match the FIDs to ROSMAP_meta_sub by row order (they should be in the same order)
  if (nrow(ROSMAP_meta_sub) == length(fid_from_meta)) {
    ROSMAP_meta_sub$FID <- fid_from_meta
  } else {
    stop("ERROR: Row count mismatch when creating FID from SAMPID. Metadata has ", 
         nrow(ROSMAP_meta_sub), " rows but SAMPID has ", length(fid_from_meta), " values")
  }
} else {
  stop("ERROR: fid_method='", opt$fid_method, "' not yet supported for early FID creation. Supported: study_projid, individual_id, subject_id_from_sampid")
}
cat("  Created", length(unique(ROSMAP_meta_sub$FID)), "unique FIDs for metadata\n")

# Normalize age
ROSMAP_meta_sub$age_death <- gsub("90\\+", "90", ROSMAP_meta_sub$age_death)
ROSMAP_meta_sub$age_death <- as.numeric(gsub("[^0-9.]", "", ROSMAP_meta_sub$age_death))

# Convert sex to numeric (0/1)
# Handle various formats: M/F, Male/Female, 0/1, GTEx format (1=male, 2=female), etc.
if (is.character(ROSMAP_meta_sub$msex) || is.factor(ROSMAP_meta_sub$msex)) {
  msex_values <- as.character(ROSMAP_meta_sub$msex)
  # Convert to uppercase for matching
  msex_upper <- toupper(msex_values)
  # Map to 0/1: 
  # - GTEx format: "2" = female (0), "1" = male (1)
  # - Standard: M/Male -> 1, F/Female -> 0
  ROSMAP_meta_sub$msex <- ifelse(msex_upper %in% c("M", "MALE", "1", "TRUE"), 1,
                                 ifelse(msex_upper %in% c("F", "FEMALE", "0", "FALSE", "2"), 0, NA))
  cat("  Converted sex column to numeric (0=Female, 1=Male)\n")
  cat("  Sex distribution:", table(ROSMAP_meta_sub$msex, useNA="ifany"), "\n")
} else {
  # Already numeric, but ensure it's 0/1
  # Handle GTEx format where 2 = female, 1 = male
  ROSMAP_meta_sub$msex <- as.numeric(ROSMAP_meta_sub$msex)
  if (any(ROSMAP_meta_sub$msex == 2, na.rm=TRUE)) {
    cat("  Detected GTEx sex format (1=male, 2=female), converting 2->0...\n")
    ROSMAP_meta_sub$msex[ROSMAP_meta_sub$msex == 2] <- 0
  }
  if (any(ROSMAP_meta_sub$msex > 1 | ROSMAP_meta_sub$msex < 0, na.rm=TRUE)) {
    cat("  WARNING: Sex values outside 0/1 range, normalizing...\n")
    ROSMAP_meta_sub$msex[ROSMAP_meta_sub$msex > 1] <- 1
    ROSMAP_meta_sub$msex[ROSMAP_meta_sub$msex < 0] <- 0
  }
}

# Create age interaction terms
ROSMAP_meta_sub$age_death2 <- ROSMAP_meta_sub$age_death^2
ROSMAP_meta_sub$age_death_sex <- ROSMAP_meta_sub$age_death * ROSMAP_meta_sub$msex
ROSMAP_meta_sub$age_death2_sex <- ROSMAP_meta_sub$age_death2 * ROSMAP_meta_sub$msex

# Remove rows with NA
ROSMAP_meta_sub_clean <- ROSMAP_meta_sub[complete.cases(ROSMAP_meta_sub), ]

# Keep unique FID
ROSMAP_meta_sub_unique <- ROSMAP_meta_sub_clean[!duplicated(ROSMAP_meta_sub_clean$FID), ]
rownames(ROSMAP_meta_sub_unique) <- ROSMAP_meta_sub_unique$FID

cat("Preparing metadata for cell proportion matching...\n")
# Prepare metadata for matching
num_unique_values <- sapply(ROSMAP_meta_temp_cleaned, function(x) length(unique(x)))
cols_to_keep <- num_unique_values > 1
ROSMAP_meta_temp_cleaned_rmVar <- ROSMAP_meta_temp_cleaned[, cols_to_keep, drop = FALSE]

# Detect specimenID column
col_specimenID <- if (!is.null(opt$col_specimenID)) {
  opt$col_specimenID
} else {
  find_column(ROSMAP_meta_temp_cleaned_rmVar, c("specimenID", "SpecimenID", "specimen_id", "synapseID", "SynapseID", 
                                                 "sampleID", "SampleID", "id", "ID"), required = TRUE)
}

# Remove duplicates
ROSMAP_meta_temp_cleaned_rmVar <- ROSMAP_meta_temp_cleaned_rmVar[!duplicated(ROSMAP_meta_temp_cleaned_rmVar[[col_specimenID]]), ]
rownames(ROSMAP_meta_temp_cleaned_rmVar) <- ROSMAP_meta_temp_cleaned_rmVar[[col_specimenID]]

# Match cell proportions with metadata
# Following reference script: select synapseID (or specimenID), individualID, Study, projid
# First, ensure we have the right column names in the metadata
# The reference script uses "synapseID" in metadata, but it might be "specimenID" in our case
synapse_col_in_meta <- ifelse("synapseID" %in% colnames(ROSMAP_meta_temp_cleaned_rmVar),
                              "synapseID",
                              ifelse(col_specimenID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar),
                                     col_specimenID, NULL))

if (is.null(synapse_col_in_meta)) {
  stop("ERROR: Cannot find synapseID or specimenID column in metadata")
}

# Build matcher columns - following reference script exactly
matcher_cols <- c(synapse_col_in_meta)
if (!is.null(col_individualID) && col_individualID %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  matcher_cols <- c(matcher_cols, col_individualID)
}
if (!is.null(col_study) && col_study %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  matcher_cols <- c(matcher_cols, col_study)
}
if (!is.null(col_projid) && col_projid %in% colnames(ROSMAP_meta_temp_cleaned_rmVar)) {
  matcher_cols <- c(matcher_cols, col_projid)
}

ROSMAP_meta_matcher <- ROSMAP_meta_temp_cleaned_rmVar[, matcher_cols, drop = FALSE]
# Rename to match reference script: synapseID, individualID, Study, projid
colnames(ROSMAP_meta_matcher)[colnames(ROSMAP_meta_matcher) == synapse_col_in_meta] <- "synapseID"
if (!is.null(col_individualID) && col_individualID %in% colnames(ROSMAP_meta_matcher) && col_individualID != "individualID") {
  colnames(ROSMAP_meta_matcher)[colnames(ROSMAP_meta_matcher) == col_individualID] <- "individualID"
}
if (!is.null(col_study) && col_study %in% colnames(ROSMAP_meta_matcher) && col_study != "Study") {
  colnames(ROSMAP_meta_matcher)[colnames(ROSMAP_meta_matcher) == col_study] <- "Study"
}
if (!is.null(col_projid) && col_projid %in% colnames(ROSMAP_meta_matcher) && col_projid != "projid") {
  colnames(ROSMAP_meta_matcher)[colnames(ROSMAP_meta_matcher) == col_projid] <- "projid"
}

cat("  Metadata matcher rows:", nrow(ROSMAP_meta_matcher), "\n")
cat("  Metadata matcher columns:", paste(colnames(ROSMAP_meta_matcher), collapse=", "), "\n")
cat("  First few synapseIDs in metadata:", paste(head(ROSMAP_meta_matcher$synapseID, 5), collapse=", "), "\n")

# Determine specimenID column name in cell proportions
specimen_col <- ifelse("specimenID" %in% colnames(ROSMAP_estimations_scaled), 
                       "specimenID",
                       ifelse("synapseID" %in% colnames(ROSMAP_estimations_scaled),
                              "synapseID", colnames(ROSMAP_estimations_scaled)[1]))

cat("  Cell proportions specimen column:", specimen_col, "\n")
cat("  Cell proportions rows:", nrow(ROSMAP_estimations_scaled), "\n")
cat("  Cell proportions columns:", paste(head(colnames(ROSMAP_estimations_scaled), 5), collapse=", "), "\n")
if (specimen_col %in% colnames(ROSMAP_estimations_scaled)) {
  cat("  First few specimenIDs in cell props:", paste(head(ROSMAP_estimations_scaled[[specimen_col]], 5), collapse=", "), "\n")
  cat("  Unique specimenIDs in cell props:", length(unique(ROSMAP_estimations_scaled[[specimen_col]])), "\n")
}

# Check overlap
if (specimen_col %in% colnames(ROSMAP_estimations_scaled)) {
  matching_ids <- intersect(ROSMAP_meta_matcher$synapseID, ROSMAP_estimations_scaled[[specimen_col]])
  cat("  Matching specimenIDs:", length(matching_ids), "out of", 
      length(unique(ROSMAP_meta_matcher$synapseID)), "metadata and",
      length(unique(ROSMAP_estimations_scaled[[specimen_col]])), "cell props\n")
  if (length(matching_ids) == 0) {
    cat("  WARNING: No matching specimenIDs! First 5 metadata IDs:", 
        paste(head(ROSMAP_meta_matcher$synapseID, 5), collapse=", "), "\n")
    cat("  First 5 cell prop IDs:", paste(head(ROSMAP_estimations_scaled[[specimen_col]], 5), collapse=", "), "\n")
  }
}

ROSMAP_estimations_combined <- merge(ROSMAP_meta_matcher, 
                                     ROSMAP_estimations_scaled, 
                                     by.x = "synapseID", 
                                     by.y = specimen_col,
                                     all.x = FALSE, all.y = TRUE)

cat("  After merge with metadata - rows:", nrow(ROSMAP_estimations_combined), "\n")
cat("  After merge - columns:", paste(head(colnames(ROSMAP_estimations_combined), 5), collapse=", "), "\n")

# Check for missing Study/projid
if ("Study" %in% colnames(ROSMAP_estimations_combined)) {
  cat("  Study column - unique values:", length(unique(ROSMAP_estimations_combined$Study)), "\n")
  cat("  Study column - NA count:", sum(is.na(ROSMAP_estimations_combined$Study)), "\n")
} else {
  cat("  WARNING: Study column not found after merge!\n")
}

if ("projid" %in% colnames(ROSMAP_estimations_combined)) {
  cat("  projid column - unique values:", length(unique(ROSMAP_estimations_combined$projid)), "\n")
  cat("  projid column - NA count:", sum(is.na(ROSMAP_estimations_combined$projid)), "\n")
} else {
  cat("  WARNING: projid column not found after merge!\n")
}

# Create FID and IID - parameterized for different cohorts
cat("Creating FID using method:", opt$fid_method, "\n")

# Helper function to create FID based on method
create_fid <- function(df, method, format_str = NULL) {
  if (method == "study_projid") {
    # ROSMAP: Study + projid
    if (!("Study" %in% colnames(df) && "projid" %in% colnames(df))) {
      stop("ERROR: fid_method='study_projid' requires Study and projid columns")
    }
    valid_rows <- !is.na(df$Study) & !is.na(df$projid)
    if (sum(valid_rows) == 0) {
      stop("ERROR: No valid rows with both Study and projid. Check specimenID matching.")
    }
    fid <- NA_character_
    fid[valid_rows] <- paste0(df$Study[valid_rows], df$projid[valid_rows])
    return(list(fid = fid, valid_rows = valid_rows))
    
  } else if (method == "individual_id") {
    # Mayo, MSBB, CommonMind: Use individualID directly
    if (!("individualID" %in% colnames(df))) {
      stop("ERROR: fid_method='individual_id' requires individualID column")
    }
    valid_rows <- !is.na(df$individualID)
    if (sum(valid_rows) == 0) {
      stop("ERROR: No valid rows with individualID")
    }
    fid <- df$individualID
    return(list(fid = fid, valid_rows = valid_rows))
    
  } else if (method == "specimen_id") {
    # MSBB (after ID matching): Use specimenID
    # Find specimenID column
    specimen_col_fid <- NULL
    for (col in c("specimenID", "SpecimenID", "specimen_id", "synapseID", "SynapseID", 
                  "sampleID", "SampleID", "id", "ID")) {
      if (col %in% colnames(df)) {
        specimen_col_fid <- col
        break
      }
    }
    if (is.null(specimen_col_fid)) {
      stop("ERROR: fid_method='specimen_id' requires specimenID column")
    }
    valid_rows <- !is.na(df[[specimen_col_fid]])
    if (sum(valid_rows) == 0) {
      stop("ERROR: No valid rows with specimenID")
    }
    fid <- df[[specimen_col_fid]]
    return(list(fid = fid, valid_rows = valid_rows))
    
  } else if (method == "submitted_subject_id") {
    # NABEC: Use submitted_subject_id
    if (!("submitted_subject_id" %in% colnames(df))) {
      stop("ERROR: fid_method='submitted_subject_id' requires submitted_subject_id column")
    }
    valid_rows <- !is.na(df$submitted_subject_id)
    if (sum(valid_rows) == 0) {
      stop("ERROR: No valid rows with submitted_subject_id")
    }
    fid <- df$submitted_subject_id
    return(list(fid = fid, valid_rows = valid_rows))
    
  } else if (method == "genotyping_sample_id") {
    # CommonMind: Use Genotyping_Sample_ID (first column or specific column)
    genotyping_col <- NULL
    for (col in c("Genotyping_Sample_ID", "genotyping_sample_id", "GenotypingSampleID")) {
      if (col %in% colnames(df)) {
        genotyping_col <- col
        break
      }
    }
    if (is.null(genotyping_col)) {
      # Try first column if it looks like an ID
      if (ncol(df) > 0) {
        genotyping_col <- colnames(df)[1]
        cat("  WARNING: Using first column as Genotyping_Sample_ID:", genotyping_col, "\n")
      } else {
        stop("ERROR: fid_method='genotyping_sample_id' requires Genotyping_Sample_ID column")
      }
    }
    valid_rows <- !is.na(df[[genotyping_col]])
    if (sum(valid_rows) == 0) {
      stop("ERROR: No valid rows with Genotyping_Sample_ID")
    }
    fid <- df[[genotyping_col]]
    return(list(fid = fid, valid_rows = valid_rows))
    
  } else if (method == "subject_id_from_sampid") {
    # GTEx: Extract subject ID from SAMPID (first two parts separated by "-")
    sampid_col <- NULL
    for (col in c("SAMPID", "sampid", "specimenID", "SpecimenID", "specimen_id", 
                  "synapseID", "SynapseID", "sampleID", "SampleID")) {
      if (col %in% colnames(df)) {
        sampid_col <- col
        break
      }
    }
    if (is.null(sampid_col)) {
      stop("ERROR: fid_method='subject_id_from_sampid' requires SAMPID or specimenID column")
    }
    valid_rows <- !is.na(df[[sampid_col]])
    if (sum(valid_rows) == 0) {
      stop("ERROR: No valid rows with SAMPID")
    }
    # Extract subject ID: split by "-" and take first two parts
    fid <- sapply(strsplit(as.character(df[[sampid_col]]), "-"), function(x) {
      if (length(x) >= 2) {
        paste(x[1], x[2], sep = "-")
      } else {
        NA_character_
      }
    })
    valid_rows <- valid_rows & !is.na(fid)
    return(list(fid = fid, valid_rows = valid_rows))
    
  } else if (method == "custom") {
    # Custom format (future implementation)
    if (is.null(format_str)) {
      stop("ERROR: fid_method='custom' requires fid_format parameter")
    }
    stop("ERROR: Custom FID format not yet implemented. Please use one of the predefined methods.")
    
  } else {
    stop("ERROR: Unknown fid_method: ", method, 
         ". Valid options: study_projid, individual_id, specimen_id, submitted_subject_id, genotyping_sample_id, subject_id_from_sampid")
  }
}

# Create FID using the specified method
fid_result <- create_fid(ROSMAP_estimations_combined, opt$fid_method, opt$fid_format)

# Apply FID and filter valid rows
ROSMAP_estimations_combined$FID <- fid_result$fid
ROSMAP_estimations_combined <- ROSMAP_estimations_combined[fid_result$valid_rows, , drop = FALSE]

cat("  Created FID using method:", opt$fid_method, "\n")
cat("  Valid rows:", sum(fid_result$valid_rows), "out of", length(fid_result$fid), "\n")
cat("  Rows after filtering:", nrow(ROSMAP_estimations_combined), "\n")

ROSMAP_estimations_combined$IID <- ROSMAP_estimations_combined$FID
cat("  Unique FIDs created:", length(unique(ROSMAP_estimations_combined$FID)), "\n")
if (length(unique(ROSMAP_estimations_combined$FID)) == 1) {
  cat("  WARNING: Only 1 unique FID! This suggests all samples have the same ID, or merge failed.\n")
  cat("  First few FIDs:", paste(head(ROSMAP_estimations_combined$FID, 5), collapse=", "), "\n")
}

# Select cell type columns (exclude metadata columns)
cell_type_cols <- setdiff(colnames(ROSMAP_estimations_combined), 
                          c("synapseID", "individualID", "Study", "projid", "FID", "IID", specimen_col))
ROSMAP_estimations_named <- ROSMAP_estimations_combined[, c("FID", "IID", cell_type_cols), drop = FALSE]

cat("Performing RINT transformation...\n")
# Debug: Check dimensions before merge
cat("  Metadata rows:", nrow(ROSMAP_meta_sub_unique), "\n")
cat("  Cell proportions rows:", nrow(ROSMAP_estimations_named), "\n")
cat("  Metadata FIDs:", length(unique(ROSMAP_meta_sub_unique$FID)), "\n")
cat("  Cell prop FIDs:", length(unique(ROSMAP_estimations_named$FID)), "\n")

# Merge with metadata for regression
combined_df <- merge(ROSMAP_meta_sub_unique, ROSMAP_estimations_named, by = "FID", all = FALSE)
cat("  After merge rows:", nrow(combined_df), "\n")

if (nrow(combined_df) == 0) {
  stop("ERROR: No matching samples found between metadata and cell proportions. Check FID matching.")
}

# Remove duplicates
combined_df <- combined_df[!duplicated(combined_df$FID), ]
cat("  After removing duplicates:", nrow(combined_df), "\n")

# Check for required columns
required_vars <- c("msex", "age_death", "age_death2", "age_death_sex", "age_death2_sex")
missing_vars <- setdiff(required_vars, colnames(combined_df))
if (length(missing_vars) > 0) {
  stop(paste("ERROR: Missing required variables:", paste(missing_vars, collapse = ", ")))
}

# Check for valid cases
valid_cases <- complete.cases(combined_df[, required_vars])
cat("  Valid cases (complete):", sum(valid_cases), "out of", length(valid_cases), "\n")

if (sum(valid_cases) == 0) {
  stop("ERROR: No valid cases after checking for complete data. All samples have missing values in required variables.")
}

# Get cell types
cell_types <- setdiff(colnames(ROSMAP_estimations_named), c("FID", "IID"))
cat("  Cell types to process:", length(cell_types), "\n")

# Perform linear regression and RINT transformation
pheno_lms <- lapply(cell_types, function(cell_type) {
  if (cell_type %in% colnames(combined_df)) {
    # Check for valid data for this cell type
    valid_data <- complete.cases(combined_df[, c(required_vars, cell_type)])
    if (sum(valid_data) == 0) {
      warning(paste("No valid cases for cell type:", cell_type))
      return(NULL)
    }
    
    # Create formula
    formula_str <- paste0("scale(", cell_type, ") ~ msex + scale(age_death) + scale(age_death_sex) + scale(age_death2) + scale(age_death2_sex)")
    
    # Fit model with only valid cases
    model_data <- combined_df[valid_data, ]
    if (nrow(model_data) == 0) {
      warning(paste("No valid data for regression for cell type:", cell_type))
      return(NULL)
    }
    
    # Fit model
    model <- lm(as.formula(formula_str), data = model_data)
    
    # Get residuals
    residuals <- residuals(model)
    residuals <- residuals[!is.na(residuals)]
    
    # Apply RINT
    if (length(residuals) > 0) {
      residuals_rint <- RNOmni::RankNorm(residuals)
    } else {
      residuals_rint <- residuals
    }
    
    # Create data frame
    results <- data.frame(
      transformed_residuals = residuals_rint,
      cell_type = cell_type,
      FID = combined_df$FID[!is.na(combined_df[[cell_type]])]
    )
    
    return(results)
  } else {
    return(NULL)
  }
})

# Combine results
pheno_lms <- bind_rows(pheno_lms)

# Pivot to wide format
pivot_df <- pivot_wider(data = pheno_lms,
                        names_from = "cell_type",
                        values_from = "transformed_residuals")

# Add IID column
pivot_df$IID <- pivot_df$FID
pivot_df <- pivot_df[, c("FID", "IID", setdiff(colnames(pivot_df), c("FID", "IID")))]

cat("Matching with WGS samples...\n")
# Load WGS PSAM file
if (!is.null(opt$wgs_psam_file) && file.exists(opt$wgs_psam_file)) {
  psam <- read.table(opt$wgs_psam_file, header = TRUE, sep = "\t")
  wgs_samples <- psam[, 1]
  
  # Match samples
  samples_df <- pivot_df[pivot_df$FID %in% wgs_samples, ]
  
  # Handle samples not directly matched
  samples_left <- wgs_samples[!wgs_samples %in% samples_df$FID]
  samples_left_df <- pivot_df[!pivot_df$FID %in% wgs_samples, ]
  
  # Try matching via biospecimen file
  # Check if file exists and has content (not empty)
  biospecimen_file_valid <- !is.null(opt$biospecimen_file) && 
                            file.exists(opt$biospecimen_file) && 
                            file.info(opt$biospecimen_file)$size > 0
  
  if (biospecimen_file_valid && nrow(samples_left_df) > 0) {
    rosmap_id_matcher <- read.csv(opt$biospecimen_file)
    if (nrow(rosmap_id_matcher) == 0) {
      cat("  WARNING: Biospecimen file is empty, skipping ID matching\n")
      samples_final <- samples_df
    } else {
      # Filter for WGS assay if column exists (e.g., ROSMAP/MSBB have this)
      if ("assay" %in% colnames(rosmap_id_matcher)) {
      rosmap_id_matcher <- rosmap_id_matcher[rosmap_id_matcher$assay == "wholeGenomeSeq", ]
      }
      
      # Filter out excluded samples if column exists
      if ("exclude" %in% colnames(rosmap_id_matcher)) {
        rosmap_id_matcher <- rosmap_id_matcher[rosmap_id_matcher$exclude == FALSE, ]
      }
      
      # Check if we still have rows after filtering
      if (nrow(rosmap_id_matcher) == 0) {
        cat("  WARNING: No valid samples in biospecimen file after filtering\n")
        samples_final <- samples_df
      } else {
        # Cohort-specific ID formatting:
        # - ROSMAP, MSBB, Mayo, GTEx, NABEC: Use specimenID directly (no prefix)
        # - CMC (MSSM, PENN, PITT): Add "0_" prefix to match PSAM format
        study_name <- if (!is.null(opt$study) && !is.na(opt$study)) toupper(opt$study) else "UNKNOWN"
        
        if (grepl("CMC", study_name) || grepl("MSSM|PENN|PITT", study_name)) {
          # CMC cohorts need "0_" prefix
        rosmap_id_matcher[[opt$biospec_col_specimen]] <- paste0("0_", rosmap_id_matcher[[opt$biospec_col_specimen]])
          cat("  Added '0_' prefix to specimenID for CMC cohort\n")
        } else if (grepl("ROSMAP|MSBB|MAYO|GTEX|NABEC", study_name)) {
          # These cohorts use specimenID directly without prefix
          cat("  Using specimenID directly (no prefix) for", study_name, "\n")
        } else {
          # Default: no prefix for unknown cohorts
          cat("  WARNING: Unknown cohort", study_name, "- using specimenID without prefix\n")
        }
      
      # Match via individualID
      FID_ind_matcher <- combined_df[, c('FID', 'individualID'), drop = FALSE]
      samples_left_df_individualID <- merge(FID_ind_matcher, samples_left_df, by = 'FID')
        # Use configurable column names for biospecimen file
        biospec_cols <- c(opt$biospec_col_individual, opt$biospec_col_specimen)
        samples_left_df_specimenID <- merge(rosmap_id_matcher[, biospec_cols, drop = FALSE],
                                          samples_left_df_individualID, by.x = opt$biospec_col_individual, by.y = 'individualID')
      
      # Update FID and IID to use specimen/genotyping ID for WGS matching
      # First save the original FID (Study+projid), specimen ID, and individual ID
      original_fid <- samples_left_df_specimenID$FID  # This is Study+projid from combined_df
      specimen_ids <- samples_left_df_specimenID[[opt$biospec_col_specimen]]
      individual_ids <- samples_left_df_specimenID[[opt$biospec_col_individual]]
      samples_left_df_specimenID$FID <- specimen_ids  # Update to MAP ID for genotyping matching
      samples_left_df_specimenID$IID <- specimen_ids
      samples_left_df_specimenID$individualID_clinical <- original_fid  # Preserve Study+projid for clinical metadata mapping
      # Remove other biospecimen columns (but keep individualID_clinical)
      biospec_cols_to_remove <- setdiff(biospec_cols, opt$biospec_col_individual)
      samples_left_df_specimenID <- samples_left_df_specimenID[, 
                                                                 !colnames(samples_left_df_specimenID) %in% biospec_cols_to_remove]
      
      # Combine matched samples
      # Ensure both dataframes have the same columns before rbind
      if (nrow(samples_df) > 0) {
        # Add individualID_clinical column to samples_df if it doesn't exist
        # For direct matches, clinical FID = phenotype FID (both are Study+projid)
        if (!"individualID_clinical" %in% colnames(samples_df)) {
          samples_df$individualID_clinical <- samples_df$FID
        }
        # Ensure column order matches
        common_cols <- intersect(colnames(samples_df), colnames(samples_left_df_specimenID))
        samples_final <- rbind(samples_df[, common_cols, drop = FALSE], 
                               samples_left_df_specimenID[, common_cols, drop = FALSE])
      } else {
        samples_final <- samples_left_df_specimenID
        }
        
        # Filter to only keep samples that are in WGS PSAM
        samples_final <- samples_final[samples_final$FID %in% wgs_samples, ]
        cat("  After filtering for WGS PSAM samples: ", nrow(samples_final), " samples\n")
      }
    }
  } else {
    samples_final <- samples_df
  }
} else {
  cat("Warning: WGS PSAM file not found. Using all samples.\n")
  samples_final <- pivot_df
}

# Save phenotype file (save to current directory for Nextflow, also copy to output_dir)
cat("Saving phenotype file...\n")
# Remove Individual_ID and individualID_clinical columns (only keep FID, IID, and cell types)
samples_final_output <- samples_final[, !colnames(samples_final) %in% c("Individual_ID", "individualID_clinical"), drop = FALSE]
write.table(samples_final_output, 
            opt$phenotype_output,
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save samples file
samples_with_phenotypes <- samples_final[, c("FID", "IID"), drop = FALSE]
samples_with_phenotypes <- rbind(c("FID", "IID"), samples_with_phenotypes)
write.table(samples_with_phenotypes,
            opt$samples_output,
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

cat("Preparing covariate file...\n")
cat("DEBUG: samples_final columns:", paste(colnames(samples_final), collapse=", "), "\n")
cat("DEBUG: samples_final rows:", nrow(samples_final), "\n")
# Prepare covariate file
if (!is.null(opt$pca_file) && file.exists(opt$pca_file)) {
  pca <- read.csv(opt$pca_file)
  
  # Create FID-specimenID mapper (following reference script line 161, 183-186)
  # For samples directly matched (samples_df): FID is already correct
  # For biospecimen-mapped samples: FID is specimenID
  # We need to create a mapper: phenotype FID -> clinical metadata FID
  
  # Build matcher dataframe from the two-step mapping done above
  # Check if samples_final has the individualID_clinical column (from biospecimen mapping)
  cat("DEBUG: Checking for individualID_clinical column...\n")
  if ("individualID_clinical" %in% colnames(samples_final)) {
    cat("DEBUG: Found individualID_clinical column! Creating matcher from it.\n")
    # Use the preserved individualID_clinical for mapping
    matcher_df_combined <- data.frame(pheno_FID = samples_final$FID,
                                      clinical_FID = samples_final$individualID_clinical)
    cat("DEBUG: First few matcher_df rows:\n")
    print(head(matcher_df_combined, 10))
    cat("DEBUG: First few pca FIDs:", head(pca$FID, 10), "\n")
    cat("DEBUG: First few combined_df FIDs:", head(combined_df$FID, 10), "\n")
  } else if (exists("samples_df") && nrow(samples_df) > 0 && exists("samples_left_df_specimenID") && nrow(samples_left_df_specimenID) > 0) {
    # Direct matches: FID = Study+projid (already matches combined_df)
    matcher_df_1 <- data.frame(pheno_FID = samples_df$FID, 
                               clinical_FID = samples_df$FID)
    # Biospecimen matches: use individualID_clinical if available
    if ("individualID_clinical" %in% colnames(samples_left_df_specimenID)) {
      matcher_df_2 <- data.frame(pheno_FID = samples_left_df_specimenID$FID,
                                 clinical_FID = samples_left_df_specimenID$individualID_clinical)
    } else {
      matcher_df_2 <- data.frame(pheno_FID = samples_left_df_specimenID$FID,
                                 clinical_FID = samples_left_df_specimenID$FID)
    }
    matcher_df_combined <- rbind(matcher_df_1, matcher_df_2)
  } else if (exists("samples_df") && nrow(samples_df) > 0) {
    matcher_df_combined <- data.frame(pheno_FID = samples_df$FID,
                                      clinical_FID = samples_df$FID)
  } else {
    matcher_df_combined <- data.frame(pheno_FID = samples_final$FID,
                                      clinical_FID = samples_final$FID)
  }
  
  # Merge PCA with phenotype FIDs, then add clinical metadata via matcher
  pca_with_mapper <- merge(pca,
                           matcher_df_combined,
                           by.x = "FID",
                           by.y = "pheno_FID",
                           all.x = FALSE, all.y = FALSE)
  
  # Merge with clinical metadata using clinical_FID
  pca_meta_combined <- merge(pca_with_mapper,
                             combined_df[, c("FID", "msex", "age_death", "age_death_sex",
                                           "age_death2", "age_death2_sex"), drop = FALSE],
                             by.x = "clinical_FID",
                             by.y = "FID",
                             all.x = TRUE, all.y = FALSE)
  
  # Clean up: FID/IID from pca_with_mapper are already correct after merge
  # Use the FID.x and IID.x columns from the merge (these are from pca_with_mapper)
  # Rename them to FID and IID
  if ("FID.x" %in% colnames(pca_meta_combined)) {
    pca_meta_combined$FID <- pca_meta_combined$FID.x
    pca_meta_combined$IID <- pca_meta_combined$IID.x
  } else {
    # FID/IID should already be correct from pca_with_mapper in the merge
  }
  # Remove clinical_FID and any duplicate FID/IID columns
  pca_meta_combined <- pca_meta_combined[, !colnames(pca_meta_combined) %in% c("clinical_FID", "FID.x", "FID.y", "IID.x", "IID.y")]
  
  # Reorder columns: FID, IID, PC1-PC10, msex, age_death, age_death_sex, age_death2, age_death2_sex
  pc_cols <- grep("^PC[0-9]+$", colnames(pca_meta_combined), value = TRUE)
  pca_meta_combined <- pca_meta_combined[, c("FID", "IID", pc_cols, "msex", "age_death", 
                                             "age_death_sex", "age_death2", "age_death2_sex")]
  
  # Save covariate file
  write.table(pca_meta_combined,
              opt$covariate_output,
              row.names = FALSE, quote = FALSE, sep = "\t")
} else {
  cat("Warning: PCA file not found. Creating covariate file without PCA.\n")
  # Create covariate file without PCA
  cov_df <- combined_df[, c("FID", "IID", "msex", "age_death", "age_death_sex", 
                           "age_death2", "age_death2_sex"), drop = FALSE]
  cov_df$IID <- cov_df$FID
  write.table(cov_df,
              opt$covariate_output,
              row.names = FALSE, quote = FALSE, sep = "\t")
}

# Also copy to output_dir if different from current directory
if (opt$output_dir != "." && opt$output_dir != getwd()) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(opt$phenotype_output, file.path(opt$output_dir, opt$phenotype_output), overwrite = TRUE)
  file.copy(opt$samples_output, file.path(opt$output_dir, opt$samples_output), overwrite = TRUE)
  file.copy(opt$covariate_output, file.path(opt$output_dir, opt$covariate_output), overwrite = TRUE)
}

cat("Phenotype and covariate preparation completed!\n")
cat("Output files saved to:", opt$output_dir, "\n")
cat("Number of samples with phenotypes:", nrow(samples_final), "\n")

