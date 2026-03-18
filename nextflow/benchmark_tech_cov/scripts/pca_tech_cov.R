#!/usr/bin/env Rscript
# PCA and Technical Covariate Computation
# Processes batch count matrices, computes PCA, and generates z-score normalized data

suppressPackageStartupMessages({
  library(optparse)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(matrixStats)
  library(broom)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(ggsignif)
  library(modelr)
  library(cowplot)
  library(gridExtra)
  library(RColorBrewer)
  library(DESeq2)
  library(factoextra)
  library(data.table)  # For faster CSV reading with fread()
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

# Helper function to validate gene IDs (Ensembl IDs)
validate_gene_ids <- function(gene_ids, context = "", min_ensembl_pct = 0.8) {
  if (is.null(gene_ids) || length(gene_ids) == 0) {
    stop("ERROR [", context, "]: Gene IDs are NULL or empty")
  }
  
  # Check if they're numeric (bad)
  if (all(grepl("^[0-9]+$", head(gene_ids, min(100, length(gene_ids)))))) {
    stop("ERROR [", context, "]: Gene IDs are numeric (1, 2, 3...) instead of Ensembl IDs. ",
         "First 10 IDs: ", paste(head(gene_ids, 10), collapse=", "))
  }
  
  # Check if they look like Ensembl IDs
  ensembl_count <- sum(grepl("^ENSG", gene_ids, ignore.case = TRUE))
  ensembl_pct <- ensembl_count / length(gene_ids)
  
  if (ensembl_pct < min_ensembl_pct) {
    warning("WARNING [", context, "]: Only ", round(ensembl_pct * 100, 1), 
            "% of gene IDs look like Ensembl IDs (expected at least ", round(min_ensembl_pct * 100), "%). ",
            "First 10 IDs: ", paste(head(gene_ids, 10), collapse=", "))
    if (ensembl_pct < 0.5) {
      stop("ERROR [", context, "]: Less than 50% of gene IDs are Ensembl IDs. ",
           "This suggests a serious problem with gene ID extraction.")
    }
  } else {
    cat("  ✓ Validated: ", round(ensembl_pct * 100, 1), "% of gene IDs are Ensembl IDs\n", sep="")
  }
  
  return(TRUE)
}

# Parse command line arguments
option_list <- list(
  make_option(c("--count_matrix_file"), type="character", default=NULL,
              help="Path to single CSV count matrix file (genes as rows, samples as columns)"),
  make_option(c("--metadata_file"), type="character", default=NULL,
              help="Path to combined metadata file (should include all technical covariates)"),
  make_option(c("--tissue_filter"), type="character", default="dorsolateral prefrontal cortex",
              help="Tissue type to filter (default: dorsolateral prefrontal cortex)"),
  make_option(c("--col_specimenID"), type="character", default=NULL,
              help="Column name for specimenID (will auto-detect if not specified)"),
  make_option(c("--col_synapseID"), type="character", default=NULL,
              help="Column name for synapseID (will auto-detect if not specified)"),
  make_option(c("--col_sample_id_for_matching"), type="character", default=NULL,
              help="Column name to use for matching count matrix sample names to metadata (e.g., 'specimenID', 'synapseID', 'id'). If not specified, will auto-detect by checking which column matches sample names."),
  make_option(c("--col_tissue"), type="character", default=NULL,
              help="Column name for tissue (will auto-detect if not specified)"),
  make_option(c("--col_assay"), type="character", default=NULL,
              help="Column name for assay (will auto-detect if not specified)"),
  make_option(c("--col_RIN"), type="character", default=NULL,
              help="Column name for RIN (optional, will auto-detect if not specified)"),
  make_option(c("--col_age"), type="character", default=NULL,
              help="Column name for age_death (will auto-detect if not specified)"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--zscore_output"), type="character", default="zscore_data.RData",
              help="Output filename for z-score data"),
  make_option(c("--metadata_output"), type="character", default="metadata_DLPFC.RData",
              help="Output filename for metadata"),
  make_option(c("--metrics_output"), type="character", default="combined_metrics.csv",
              help="Output filename for combined metrics")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Loading count matrix from CSV file...\n")
# Load single count matrix file
# Use fread() for faster reading of large CSV files (much faster than read_csv for large files)
# fread() with data.table = FALSE returns a regular data.frame, compatible with existing code
# Assume first column is gene identifier, rest are samples
file_size_mb <- file.info(opt$count_matrix_file)$size / (1024^2)
if (file_size_mb > 100) {
  cat("  Large file detected (", round(file_size_mb, 1), " MB), using fread() for faster loading...\n", sep = "")
  count_matrix_raw <- fread(opt$count_matrix_file, data.table = FALSE)
} else {
  cat("  Using read_csv() for smaller file...\n")
  count_matrix_raw <- read_csv(opt$count_matrix_file, show_col_types = FALSE)
  # Convert tibble to data.frame for consistency
  count_matrix_raw <- as.data.frame(count_matrix_raw)
}

# Check for gene identifier column
# First, check for common gene ID column names (gene_id, ensembl_id, etc.)
potential_gene_cols <- c("gene_id", "Gene_ID", "ensembl_id", "Ensembl_ID", 
                         "gene", "Gene", "gene_symbol", "GeneSymbol", "gene_name", "GeneName")
gene_col <- NULL

for (col in potential_gene_cols) {
  if (col %in% colnames(count_matrix_raw)) {
    col_values <- head(count_matrix_raw[[col]], 100)
    # Check if values look like Ensembl IDs or gene symbols
    if (any(grepl("^ENSG", col_values, ignore.case = TRUE)) || 
        any(grepl("^[A-Za-z]", col_values))) {
      gene_col <- col
      cat("Found gene identifier column:", gene_col, "\n")
      cat("  First 5 gene IDs:", paste(head(count_matrix_raw[[gene_col]], 5), collapse=", "), "\n")
      break
    }
  }
}

# If not found, check first column
if (is.null(gene_col)) {
  first_col <- colnames(count_matrix_raw)[1]
  cat("  First column name:", first_col, "\n")
  cat("  First column class:", class(count_matrix_raw[[1]]), "\n")
  cat("  First 5 values in first column:", paste(head(count_matrix_raw[[1]], 5), collapse=", "), "\n")
  
  first_col_values <- head(count_matrix_raw[[1]], 100)
  looks_like_gene_ids <- any(grepl("^ENSG", first_col_values, ignore.case = TRUE)) || 
                         any(grepl("^[A-Za-z]", first_col_values)) ||
                         (!is.numeric(first_col_values) && !all(is.na(first_col_values)))
  
  if (looks_like_gene_ids) {
    gene_col <- first_col
    cat("Using first column '", gene_col, "' as gene identifier\n")
  } else {
    cat("First column does not appear to be gene identifiers\n")
    # Check second column as fallback
    if (ncol(count_matrix_raw) > 1) {
      second_col <- colnames(count_matrix_raw)[2]
      second_col_values <- head(count_matrix_raw[[second_col]], 100)
      if (any(grepl("^ENSG", second_col_values, ignore.case = TRUE))) {
        gene_col <- second_col
        cat("Using second column '", gene_col, "' as gene identifier (contains Ensembl IDs)\n")
      }
    }
  }
}

if (is.null(gene_col)) {
  stop("ERROR: Cannot find gene identifier column with Ensembl IDs. Checked columns: ", 
       paste(colnames(count_matrix_raw)[1:min(5, ncol(count_matrix_raw))], collapse=", "))
}

sample_cols <- setdiff(colnames(count_matrix_raw), gene_col)

# Extract count matrix (samples as columns)
# CRITICAL: Convert to data.frame first (tibbles don't support rownames properly)
combined_rnaseq_mat <- as.data.frame(count_matrix_raw[, sample_cols, drop = FALSE])
# Set rownames from gene identifier column
if (!is.null(gene_col) && gene_col %in% colnames(count_matrix_raw)) {
  gene_ids_from_col <- as.character(count_matrix_raw[[gene_col]])
  rownames(combined_rnaseq_mat) <- gene_ids_from_col
  cat("  Set rownames from column:", gene_col, "\n")
  cat("  First 5 rownames:", paste(head(rownames(combined_rnaseq_mat), 5), collapse=", "), "\n")
  
  # CHECKPOINT 1: Validate gene IDs after initial extraction
  validate_gene_ids(rownames(combined_rnaseq_mat), "After setting rownames from gene_id column")
} else {
  stop("ERROR: Cannot find gene identifier column. Cannot proceed.")
}

# Get sample names from column names
sample_names <- colnames(combined_rnaseq_mat)
cat("Loaded count matrix:", nrow(combined_rnaseq_mat), "genes x", ncol(combined_rnaseq_mat), "samples\n")

cat("Loading metadata...\n")
# Load metadata
# Use fread() for large files, and choose the correct reader for smaller
# delimited text files based on the file extension.
file_size_mb <- file.info(opt$metadata_file)$size / (1024^2)
if (file_size_mb > 10) {
  cat("  Large metadata file detected (", round(file_size_mb, 1), " MB), using fread()...\n", sep = "")
  metadata_sep <- if (grepl("\\.(tsv|txt)$", opt$metadata_file, ignore.case = TRUE)) "\t" else ","
  rosmap_meta <- fread(opt$metadata_file, data.table = FALSE, sep = metadata_sep)
} else {
  if (grepl("\\.(tsv|txt)$", opt$metadata_file, ignore.case = TRUE)) {
    cat("  Using read_tsv() for smaller metadata file...\n")
    rosmap_meta <- read_tsv(opt$metadata_file, show_col_types = FALSE)
  } else {
    cat("  Using read_csv() for smaller metadata file...\n")
    rosmap_meta <- read_csv(opt$metadata_file, show_col_types = FALSE)
  }
  # Convert tibble to data.frame for consistency
  rosmap_meta <- as.data.frame(rosmap_meta)
}

# Optionally filter to geneExpression rows if the metadata dataType column exists
if ("dataType" %in% colnames(rosmap_meta)) {
  cat("Filtering metadata to dataType == 'geneExpression'\n")
  rosmap_meta <- rosmap_meta %>% filter(tolower(dataType) == 'geneexpression')
}

# Detect column names with flexible matching
# If col_sample_id_for_matching is provided, specimenID is optional (for cohorts like GTEx that use different column names)
col_specimenID <- if (!is.null(opt$col_specimenID)) {
  opt$col_specimenID
} else {
  # Make required = FALSE if col_sample_id_for_matching is provided (e.g., GTEx uses SAMPID)
  required_specimenID <- is.null(opt$col_sample_id_for_matching)
  find_column(rosmap_meta, c("specimenID", "SpecimenID", "specimen_id", "sampleID", "SampleID", "sample_id"), required = required_specimenID)
}

col_tissue <- if (!is.null(opt$col_tissue)) {
  opt$col_tissue
} else {
  # Tissue column is optional for some cohorts (GTEx metadata may not have standard tissue column)
  find_column(rosmap_meta, c("tissue", "Tissue", "tissueType", "TissueType", "tissue_type", "SMTSD"), required = FALSE)
}

col_assay <- if (!is.null(opt$col_assay)) {
  opt$col_assay
} else {
  find_column(rosmap_meta, c("assay", "Assay", "assayType", "AssayType", "assay_type"))
}

col_age <- if (!is.null(opt$col_age)) {
  opt$col_age
} else {
  find_column(rosmap_meta, c("age_death", "ageDeath", "age_at_death", "Age_death", "AgeDeath", "Age_at_death"))
}

col_RIN <- if (!is.null(opt$col_RIN)) {
  opt$col_RIN
} else {
  find_column(rosmap_meta, c("RIN", "rin", "RnaIntegrityNumber", "RNA_integrity_number"))
}

cat("Detected column names:\n")
cat("  specimenID:", col_specimenID, "\n")
cat("  tissue:", col_tissue, "\n")
cat("  assay:", col_assay, "\n")
cat("  age:", col_age, "\n")
cat("  RIN:", ifelse(is.null(col_RIN), "NOT FOUND (optional)", col_RIN), "\n")

# Metadata file should already contain all technical covariates
# Filter by assay if column exists
if (!is.null(col_assay) && col_assay %in% colnames(rosmap_meta)) {
  rosmap_meta <- rosmap_meta %>% filter(!!sym(col_assay) == 'rnaSeq')
}

# Find column to use for matching count matrix sample names to metadata
# This is cohort-specific: ROSMAP uses synapseID, Mayo uses specimenID, etc.
if (!is.null(opt$col_sample_id_for_matching)) {
  # User explicitly specified which column to use
  if (!opt$col_sample_id_for_matching %in% colnames(rosmap_meta)) {
    stop("ERROR: Specified column '", opt$col_sample_id_for_matching, "' not found in metadata. Available columns: ", 
         paste(head(colnames(rosmap_meta), 10), collapse=", "), "...")
  }
  col_synapseID <- opt$col_sample_id_for_matching
  cat("  Using specified column for sample matching:", col_synapseID, "\n")
  # Create synapseID column from the specified column for consistency
  rosmap_meta$synapseID <- rosmap_meta[[col_synapseID]]
  col_synapseID <- "synapseID"
} else {
  # Auto-detect: try to find which column matches sample names
  # Priority: 1) Check if specimenID matches (Mayo), 2) Check synapseID/id (ROSMAP), 3) Fallback to specimenID
  if (!is.null(col_specimenID) && col_specimenID %in% colnames(rosmap_meta)) {
    specimen_matches <- sum(rosmap_meta[[col_specimenID]] %in% sample_names, na.rm = TRUE)
    cat("  Checking specimenID matches:", specimen_matches, "out of", length(sample_names), "sample names\n")
    
    if (specimen_matches > 0) {
      cat("  Using specimenID for sample matching (", specimen_matches, " matches found)\n", sep = "")
      rosmap_meta$synapseID <- rosmap_meta[[col_specimenID]]
      col_synapseID <- "synapseID"
    } else {
      # Try to find synapseID column
      col_synapseID_candidate <- find_column(rosmap_meta, c("synapseID", "SynapseID", "synapse_id", "id"))
      if (!is.null(col_synapseID_candidate) && col_synapseID_candidate %in% colnames(rosmap_meta)) {
        synapse_matches <- sum(rosmap_meta[[col_synapseID_candidate]] %in% sample_names, na.rm = TRUE)
        cat("  Checking ", col_synapseID_candidate, " matches:", synapse_matches, "\n", sep = "")
        if (synapse_matches > 0) {
          col_synapseID <- col_synapseID_candidate
          rosmap_meta$synapseID <- rosmap_meta[[col_synapseID]]
          col_synapseID <- "synapseID"
        } else {
          # Fallback to specimenID even if no matches (will error later if truly no matches)
          cat("  WARNING: No matches found, falling back to specimenID\n")
          rosmap_meta$synapseID <- rosmap_meta[[col_specimenID]]
          col_synapseID <- "synapseID"
        }
      } else {
        # Fallback to specimenID
        cat("  No synapseID column found, using specimenID\n")
        rosmap_meta$synapseID <- rosmap_meta[[col_specimenID]]
        col_synapseID <- "synapseID"
      }
    }
  } else {
    # Try to find synapseID column
    col_synapseID <- find_column(rosmap_meta, c("synapseID", "SynapseID", "synapse_id", "id"))
    if (is.null(col_synapseID) || !col_synapseID %in% colnames(rosmap_meta)) {
      stop("Cannot determine sample ID column for matching with count matrices. Please specify --col_sample_id_for_matching")
    }
  }
}

# Filter to samples in count matrices
rosmap_meta <- rosmap_meta[rosmap_meta[[col_synapseID]] %in% sample_names, ]

# Normalize age if column exists
if (!is.null(col_age) && col_age %in% colnames(rosmap_meta)) {
  rosmap_meta$Age_norm <- gsub("90+", "90", rosmap_meta[[col_age]])
  rosmap_meta$Age_norm <- as.numeric(gsub("[+]", "", rosmap_meta$Age_norm))
}

# Save combined metrics (save to current directory for Nextflow, also copy to output_dir)
write_csv(rosmap_meta, opt$metrics_output)
if (opt$output_dir != "." && opt$output_dir != getwd()) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(opt$metrics_output, file.path(opt$output_dir, opt$metrics_output), overwrite = TRUE)
}

# Filter for specified tissue
cat("Filtering for tissue:", opt$tissue_filter, "\n")
use_tissues <- c(opt$tissue_filter)
if (!is.null(col_tissue) && col_tissue %in% colnames(rosmap_meta)) {
  use_tissue_sample_ids <- rosmap_meta %>% 
    filter(!!sym(col_tissue) %in% use_tissues) %>%
    pull(!!sym(col_synapseID))
  
  rosmap_meta_filtered <- rosmap_meta %>% filter(!!sym(col_tissue) %in% use_tissues)
} else {
  cat("Warning: Tissue column not found, using all samples\n")
  use_tissue_sample_ids <- rosmap_meta[[col_synapseID]]
  rosmap_meta_filtered <- rosmap_meta
}

# Some metadata sources (e.g. GVEX manifest) can contain multiple rows per
# subject/sample ID. Keep the first row for each matched sample before setting
# row names so downstream metadata/count alignment remains well-defined.
if (anyDuplicated(use_tissue_sample_ids) > 0) {
  n_dup_ids <- sum(duplicated(use_tissue_sample_ids))
  cat("Warning:", n_dup_ids, "duplicate metadata rows found for sample IDs; keeping first occurrence per ID\n")
  rosmap_meta_filtered <- rosmap_meta_filtered[!duplicated(rosmap_meta_filtered[[col_synapseID]]), , drop = FALSE]
  use_tissue_sample_ids <- rosmap_meta_filtered[[col_synapseID]]
}
rownames(rosmap_meta_filtered) <- rosmap_meta_filtered[[col_synapseID]]

# Filter count matrix to match tissue-filtered samples
cat("Filtering count matrix to tissue-specific samples...\n")
# Find intersection of sample names
common_samples <- intersect(colnames(combined_rnaseq_mat), use_tissue_sample_ids)
if (length(common_samples) == 0) {
  stop("No common samples found between count matrix and metadata. Please check sample ID matching.")
}
cat("Found", length(common_samples), "samples matching tissue filter\n")

# Store rownames before filtering (in case they get lost)
gene_ids_before_filter <- rownames(combined_rnaseq_mat)
combined_rnaseq_mat <- combined_rnaseq_mat[, common_samples, drop = FALSE]
# Ensure rownames are preserved after filtering
if (is.null(rownames(combined_rnaseq_mat)) || length(rownames(combined_rnaseq_mat)) == 0 || 
    all(rownames(combined_rnaseq_mat) == as.character(1:nrow(combined_rnaseq_mat)))) {
  if (length(gene_ids_before_filter) == nrow(combined_rnaseq_mat)) {
    rownames(combined_rnaseq_mat) <- gene_ids_before_filter
    cat("  Restored rownames after filtering. First 5:", paste(head(rownames(combined_rnaseq_mat), 5), collapse=", "), "\n")
  }
}

# CHECKPOINT 2: Validate gene IDs after filtering
validate_gene_ids(rownames(combined_rnaseq_mat), "After filtering to tissue-specific samples")
# Also filter metadata to match
rosmap_meta_filtered <- rosmap_meta_filtered[common_samples, , drop = FALSE]

cat("Creating DESeq2 object...\n")
# CRITICAL: Ensure combined_rnaseq_mat is a data.frame (not tibble) to preserve rownames
# Store current rownames before any conversion
current_rownames <- rownames(combined_rnaseq_mat)
if (inherits(combined_rnaseq_mat, "tbl_df") || inherits(combined_rnaseq_mat, "tbl")) {
  cat("  Converting tibble to data.frame to preserve rownames...\n")
  combined_rnaseq_mat <- as.data.frame(combined_rnaseq_mat)
  # Re-set rownames after conversion (they might be lost)
  if (length(current_rownames) == nrow(combined_rnaseq_mat)) {
    rownames(combined_rnaseq_mat) <- current_rownames
    cat("  Re-set rownames after tibble conversion. First 5:", paste(head(rownames(combined_rnaseq_mat), 5), collapse=", "), "\n")
  } else if (exists("gene_ids_before_filter") && length(gene_ids_before_filter) == nrow(combined_rnaseq_mat)) {
    rownames(combined_rnaseq_mat) <- gene_ids_before_filter
    cat("  Re-set rownames from stored gene_ids_before_filter. First 5:", paste(head(rownames(combined_rnaseq_mat), 5), collapse=", "), "\n")
  } else {
    stop("ERROR: Cannot restore rownames after tibble conversion. current_rownames length: ", length(current_rownames), 
         ", matrix rows: ", nrow(combined_rnaseq_mat))
  }
}

# CRITICAL: Store gene IDs (rownames) BEFORE converting to matrix
gene_ids_count <- rownames(combined_rnaseq_mat)
cat("  Gene IDs from combined_rnaseq_mat length:", length(gene_ids_count), "\n")
if (length(gene_ids_count) > 0) {
  cat("  First 5 gene IDs:", paste(head(gene_ids_count, 5), collapse=", "), "\n")
}

# CHECKPOINT 3: Validate gene IDs before creating DESeq2 object
validate_gene_ids(gene_ids_count, "Before creating DESeq2 object")

# Convert to matrix - ensure rownames are preserved
count_matrix <- as.matrix(combined_rnaseq_mat)
# Explicitly set rownames (as.matrix should preserve them from data.frame, but be safe)
rownames(count_matrix) <- gene_ids_count
cat("  count_matrix rownames length:", length(rownames(count_matrix)), "\n")
cat("  count_matrix first 5 rownames:", paste(head(rownames(count_matrix), 5), collapse=", "), "\n")

# DESeq2 requires numeric integer counts. If the matrix is character, the gene column or
# a header row may have been included as data, or a sample column was read as text.
if (typeof(count_matrix) == "character" || mode(count_matrix) == "character") {
  stop("ERROR: Count matrix contains character/string values. DESeq2 requires integer counts. ",
       "Check that (1) the first column is gene IDs and is used only as rownames (not as data), ",
       "(2) column names are not read as a data row, and (3) all count columns are numeric.")
}

# Remove samples with NA RIN if RIN column exists
if (!is.null(col_RIN) && col_RIN %in% colnames(rosmap_meta_filtered)) {
  cat("Filtering samples with NA RIN...\n")
  rosmap_meta_cleaned <- rosmap_meta_filtered[!is.na(rosmap_meta_filtered[[col_RIN]]), ]
  count_matrix_cleaned <- count_matrix[, !(colnames(count_matrix) %in% 
                                           rosmap_meta_filtered[is.na(rosmap_meta_filtered[[col_RIN]]), ][[col_synapseID]])]
  # Preserve rownames after subsetting
  rownames(count_matrix_cleaned) <- rownames(count_matrix)
  cat("Samples before RIN filter:", ncol(count_matrix), "\n")
  cat("Samples after RIN filter:", ncol(count_matrix_cleaned), "\n")
  count_matrix <- count_matrix_cleaned
} else {
  cat("RIN column not found, using all samples\n")
  rosmap_meta_cleaned <- rosmap_meta_filtered
}

# DESeq2 requires integer counts; CSV read often yields numeric/double (e.g. 1.0, 2.0)
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

# Create DESeq2 object - rownames should be preserved
dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = rosmap_meta_cleaned, 
                              design = ~1)

# Verify rownames are in DESeq2 object
cat("  DESeq2 object rownames length:", length(rownames(dds)), "\n")
if (length(rownames(dds)) == 0) {
  stop("ERROR: DESeq2 object has no rownames. Cannot proceed.")
}
cat("  First 5 DESeq2 rownames:", paste(head(rownames(dds), 5), collapse=", "), "\n")

# CHECKPOINT 4: Validate gene IDs in DESeq2 object
validate_gene_ids(rownames(dds), "In DESeq2 object after creation")

# Filter genes with low variance
cat("Filtering genes with low variance...\n")
dds_hasvar <- dds[rowVars(assay(dds)) > 0.1, ]

# Perform normalization
cat("Performing DESeq2 normalization...\n")
dds_hasvar <- DESeq(dds_hasvar)
rld <- vst(dds_hasvar, blind = FALSE)

# Center genes
cat("Centering and z-score transforming...\n")
# CRITICAL: Get gene IDs from DESeq2 object (rownames should be preserved)
gene_ids <- rownames(dds_hasvar)
if (is.null(gene_ids) || length(gene_ids) == 0) {
  # Try from rld object
  gene_ids <- rownames(rld)
}
if (is.null(gene_ids) || length(gene_ids) == 0) {
  stop("ERROR: Cannot get gene IDs from DESeq2 object. rownames(dds_hasvar) and rownames(rld) are both empty.")
}
cat("  Gene IDs (rownames) length:", length(gene_ids), "\n")
cat("  First 5 gene IDs:", paste(head(gene_ids, 5), collapse=", "), "\n")

rld_assay <- assay(rld)
# Set rownames on assay (it might not have them)
rownames(rld_assay) <- gene_ids
cat("  rld_assay rownames length:", length(rownames(rld_assay)), "\n")

centered_data <- scale(rld_assay, center = TRUE, scale = FALSE)
# Preserve rownames after scale
rownames(centered_data) <- gene_ids

# Z-score transform per sample
zscore_data <- t(scale(t(centered_data)))
# Preserve rownames after transpose and scale
rownames(zscore_data) <- gene_ids

cat("  zscore_data dimensions:", dim(zscore_data), "\n")
cat("  zscore_data rownames length:", length(rownames(zscore_data)), "\n")
cat("  First 5 rownames:", paste(head(rownames(zscore_data), 5), collapse=", "), "\n")

# Save outputs (save to current directory for Nextflow, also copy to output_dir)
cat("Saving outputs...\n")
# Verify rownames before saving
if (is.null(rownames(zscore_data)) || length(rownames(zscore_data)) == 0) {
  stop("ERROR: zscore_data has no rownames before saving. Cannot proceed.")
}

# CHECKPOINT 5: Final validation before saving zscore_data
validate_gene_ids(rownames(zscore_data), "Before saving zscore_data (FINAL CHECKPOINT)")
cat("  ✓ All gene ID validations passed. Saving zscore_data with ", length(rownames(zscore_data)), " genes.\n", sep="")
saveRDS(zscore_data, opt$zscore_output)
saveRDS(rosmap_meta_filtered, opt$metadata_output)

# Also copy to output_dir if different from current directory
if (opt$output_dir != "." && opt$output_dir != getwd()) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(opt$zscore_output, file.path(opt$output_dir, opt$zscore_output), overwrite = TRUE)
  file.copy(opt$metadata_output, file.path(opt$output_dir, opt$metadata_output), overwrite = TRUE)
}

cat("PCA and technical covariate computation completed!\n")
cat("Output files saved to current directory and", opt$output_dir, "\n")

