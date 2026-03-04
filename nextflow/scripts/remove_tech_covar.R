#!/usr/bin/env Rscript
# Remove Batch Effects and Technical Covariates
# Identifies top tech covariates and removes batch effects

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
  library(PCAtools)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--zscore_data"), type="character", default=NULL,
              help="Path to z-score data RData file"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Path to metadata RData file"),
  make_option(c("--top_n_tech_cov"), type="integer", default=20,
              help="Number of top technical covariates to remove (default: 20)"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--corrected_output"), type="character", default="corrected_data.RData",
              help="Output filename for corrected data"),
  make_option(c("--metadata_output"), type="character", default="metadata_cleaned.csv",
              help="Output filename for cleaned metadata"),
  make_option(c("--tech_cov_output"), type="character", default="top_tech_covariates.txt",
              help="Output filename for top tech covariates list")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Loading z-score data and metadata...\n")
# Load data
zscore_data <- readRDS(opt$zscore_data)

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

# CRITICAL: Store gene IDs (rownames) IMMEDIATELY after loading, before any operations
# Check if rownames exist
if (is.null(rownames(zscore_data)) || length(rownames(zscore_data)) == 0) {
  # Try dimnames
  if (!is.null(dimnames(zscore_data)[[1]]) && length(dimnames(zscore_data)[[1]]) > 0) {
    gene_ids_original <- dimnames(zscore_data)[[1]]
    rownames(zscore_data) <- gene_ids_original
  } else {
    stop("ERROR: Loaded zscore_data has no rownames or dimnames. Cannot proceed.")
  }
} else {
  gene_ids_original <- rownames(zscore_data)
}
cat("  Original gene IDs (rownames) length:", length(gene_ids_original), "\n")
cat("  First 5 gene IDs:", paste(head(gene_ids_original, 5), collapse=", "), "\n")

# CHECKPOINT: Validate gene IDs immediately after loading
validate_gene_ids(gene_ids_original, "After loading zscore_data")

# Handle metadata - could be RData or CSV
if (endsWith(opt$metadata, ".RData") || endsWith(opt$metadata, ".rds")) {
  rosmap_meta <- readRDS(opt$metadata)
} else {
  rosmap_meta <- read.csv(opt$metadata)
}

# Prepare metadata for PCA
cat("Preparing metadata...\n")
# Remove tech cov (columns) with any NA values
# Transpose to check each column for completeness, then use to subset columns
cols_without_na <- complete.cases(t(rosmap_meta))
rosmap_meta_cleaned <- rosmap_meta[, cols_without_na, drop = FALSE]

# Remove tech cov with only one value
num_unique_values <- sapply(rosmap_meta_cleaned, function(x) length(unique(x)))
cols_to_keep <- num_unique_values > 1
rosmap_meta_cleaned_rmVar <- rosmap_meta_cleaned[, cols_to_keep]

# Set row names - try multiple ID columns (cohort-specific)
# Priority: synapseID (ROSMAP), id, specimenID (Mayo/MSBB)
if ("synapseID" %in% colnames(rosmap_meta_cleaned_rmVar)) {
  rownames(rosmap_meta_cleaned_rmVar) <- rosmap_meta_cleaned_rmVar$synapseID
  cat("  Set rownames from synapseID column\n")
} else if ("id" %in% colnames(rosmap_meta_cleaned_rmVar)) {
  rownames(rosmap_meta_cleaned_rmVar) <- rosmap_meta_cleaned_rmVar$id
  cat("  Set rownames from id column\n")
} else if ("specimenID" %in% colnames(rosmap_meta_cleaned_rmVar)) {
  rownames(rosmap_meta_cleaned_rmVar) <- rosmap_meta_cleaned_rmVar$specimenID
  cat("  Set rownames from specimenID column\n")
} else {
  stop("ERROR: Cannot set rownames - need synapseID, id, or specimenID column")
}

# Get tech covariate list
tech_covar_list <- colnames(rosmap_meta_cleaned_rmVar)
tech_covar_list <- tech_covar_list[!tech_covar_list %in% c("X", "specimenID", "id", "synapseID", "tissue", "assay", "organ",
                                                             "Started.job.on", "Started.mapping.on", "Finished.on")]

cat("Preparing metadata for PCA analysis...\n")
# Prepare metadata - preserve row names when converting to data.frame
metadata <- as.data.frame(rosmap_meta_cleaned_rmVar)[, tech_covar_list, drop = FALSE]
# Ensure row names are preserved
rownames(metadata) <- rownames(rosmap_meta_cleaned_rmVar)

# Convert percentage columns to numeric
percent_cols <- c("Uniquely.mapped.reads..", "Mismatch.rate.per.base...", 
                  "Deletion.rate.per.base", "Insertion.rate.per.base",
                  "X..of.reads.mapped.to.multiple.loci", 
                  "X..of.reads.mapped.to.too.many.loci",
                  "X..of.reads.unmapped..too.short", 
                  "X..of.reads.unmapped..other")

for (col in percent_cols) {
  if (col %in% colnames(metadata)) {
    metadata[[col]] <- as.numeric(sub("%", "", metadata[[col]], fixed = TRUE))
  }
}

# Coerce to factor or numeric
for (tech_cov in tech_covar_list) {
  if (tech_cov %in% colnames(metadata)) {
    if (is.numeric(metadata[[tech_cov]])) {
      metadata[[tech_cov]] <- as.numeric(metadata[[tech_cov]])
    } else {
      metadata[[tech_cov]] <- factor(metadata[[tech_cov]])
    }
  }
}

# Filter expression data to match metadata
cat("Matching samples between zscore_data and metadata...\n")
cat("  zscore_data samples:", length(colnames(zscore_data)), "\n")
cat("  metadata samples:", length(rownames(metadata)), "\n")

# Find common samples
common_samples <- intersect(colnames(zscore_data), rownames(metadata))
cat("  common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("ERROR: No common samples found between zscore_data and metadata. Check sample IDs.")
}

# Filter both to common samples and ensure same order
zscore_data <- zscore_data[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Verify they match
if (!all(colnames(zscore_data) == rownames(metadata))) {
  stop("ERROR: Sample names still don't match after filtering. This should not happen.")
}

cat("  Successfully matched", length(common_samples), "samples\n")

cat("Performing PCA to identify top technical covariates...\n")
# Perform PCA
p <- PCAtools::pca(zscore_data, metadata = metadata)

# Custom eigencor function
eigencorvalue <- function(pcaobj, components = getComponents(pcaobj, seq_len(10)), 
                          metavars, corUSE = 'pairwise.complete.obs', 
                          corFUN = 'pearson') {
  data <- pcaobj$rotated
  metadata <- pcaobj$metadata
  corvals <- list()
  
  for (i in seq_len(length(components))) {
    if (!is.numeric(data[, components[i]])) {
      warning(components[i], ' is not numeric - please check the source data')
    }
  }
  
  for (i in seq_len(length(metavars))) {
    if (!is.numeric(metadata[, metavars[i]])) {
      warning(metavars[i], ' is not numeric - please check the source data')
    }
  }
  
  xvals <- data.matrix(data[, which(colnames(data) %in% components), drop = FALSE])
  yvals <- metadata[, which(colnames(metadata) %in% metavars), drop = FALSE]
  
  # Coerce non-numeric variables to numeric
  character_columns <- !unlist(lapply(yvals, is.numeric))
  character_columns <- names(which(character_columns))
  for (c in character_columns) {
    yvals[, c] <- as.numeric(as.factor(yvals[, c]))
  }
  
  yvals <- data.matrix(yvals)
  cor_matrix <- cor(xvals, yvals, use = corUSE, method = corFUN)
  
  corvals$x_names <- colnames(xvals)
  corvals$y_names <- colnames(yvals)
  corvals$correlation_values <- cor_matrix
  
  return(corvals)
}

# Calculate correlations
all_corvals <- eigencorvalue(pcaobj = p, metavars = tech_covar_list, corFUN = "pearson")

# Extract top tech covariates
correlation_matrix <- all_corvals$correlation_values
abs_correlation_values <- abs(correlation_matrix)
cor_vector <- as.vector(abs_correlation_values)

# Find top N covariates
top_n_indices <- order(cor_vector, decreasing = TRUE)[1:(opt$top_n_tech_cov * length(all_corvals$x_names))]
all_covals_indices <- order(cor_vector, decreasing = TRUE)

x_names <- all_corvals$x_names
y_names <- all_corvals$y_names

all_covariates <- data.frame(
  x_var = rep(x_names, each = length(y_names)),
  y_var = rep(y_names, times = length(x_names)),
  correlation = cor_vector,
  stringsAsFactors = FALSE
)

all_covariates <- all_covariates[all_covals_indices, ]
all_covariates_sorted <- unique(all_covariates$y_var)
top_tech_cov <- all_covariates_sorted[1:opt$top_n_tech_cov]

# Save top tech covariates (save to current directory for Nextflow, also copy to output_dir)
write.table(top_tech_cov, opt$tech_cov_output, 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
if (opt$output_dir != "." && opt$output_dir != getwd()) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(opt$tech_cov_output, file.path(opt$output_dir, opt$tech_cov_output), overwrite = TRUE)
}

cat("Removing batch effects and technical covariates...\n")
# Prepare metadata for batch removal
metadata_nobatch <- metadata[, !colnames(metadata) %in% c("sequencingBatch", "notes", "libraryPrep"), drop = FALSE]

# Ensure all columns are numeric
non_numeric_columns <- names(which(!sapply(metadata_nobatch, function(x) all(is.numeric(x)))))
for (i in seq_len(length(non_numeric_columns))) {
  if (non_numeric_columns[i] %in% colnames(metadata_nobatch)) {
    metadata_nobatch[, non_numeric_columns[i]] <- as.numeric(as.factor(metadata_nobatch[, non_numeric_columns[i]]))
  }
}

# Select top tech covariates
metadata_nobatch <- metadata_nobatch[, colnames(metadata_nobatch) %in% top_tech_cov, drop = FALSE]

# Remove batch effects
# CRITICAL: Use the original gene_ids we stored at the beginning
# removeBatchEffect may lose rownames, so we'll restore them afterwards
zscore_data_df <- as.data.frame(zscore_data)
# Ensure rownames are set (use original gene_ids)
if (exists("gene_ids_original") && length(gene_ids_original) == nrow(zscore_data_df)) {
  rownames(zscore_data_df) <- gene_ids_original
} else if (!is.null(rownames(zscore_data))) {
  rownames(zscore_data_df) <- rownames(zscore_data)
} else if (!is.null(dimnames(zscore_data)[[1]])) {
  rownames(zscore_data_df) <- dimnames(zscore_data)[[1]]
}
cat("  zscore_data_df dimensions:", dim(zscore_data_df), "\n")
cat("  zscore_data_df rownames length:", length(rownames(zscore_data_df)), "\n")
cat("  First 5 rownames:", paste(head(rownames(zscore_data_df), 5), collapse=", "), "\n")

# Get batch columns if they exist
batch1 <- NULL
batch2 <- NULL
if ("sequencingBatch" %in% colnames(metadata)) {
  batch1 <- as.vector(metadata[, "sequencingBatch"])
}
if ("libraryPrep" %in% colnames(metadata)) {
  batch2 <- as.vector(metadata[, "libraryPrep"])
}

if (any(is.na(zscore_data_df)) == FALSE && any(is.na(metadata_nobatch)) == FALSE) {
  cat("Applying removeBatchEffect...\n")
  if (!is.null(batch1) && !is.null(batch2)) {
    zscore_data_removedBatchEff_cov <- removeBatchEffect(
      x = zscore_data_df,
      batch = batch1,
      batch2 = batch2,
      covariates = metadata_nobatch
    )
  } else if (!is.null(batch1)) {
    zscore_data_removedBatchEff_cov <- removeBatchEffect(
      x = zscore_data_df,
      batch = batch1,
      covariates = metadata_nobatch
    )
  } else {
    zscore_data_removedBatchEff_cov <- removeBatchEffect(
      x = zscore_data_df,
      covariates = metadata_nobatch
    )
  }
} else {
  stop("NA values detected in zscore data or metadata. Cannot proceed.")
}

# Ensure rownames are preserved in the corrected data
# ALWAYS restore from gene_ids_original (stored at the beginning)
cat("Checking rownames in corrected data before saving...\n")
if (exists("gene_ids_original") && length(gene_ids_original) == nrow(zscore_data_removedBatchEff_cov)) {
  cat("  Restoring rownames from original gene_ids_original...\n")
  rownames(zscore_data_removedBatchEff_cov) <- gene_ids_original
  cat("  Rownames restored. Length:", length(rownames(zscore_data_removedBatchEff_cov)), "\n")
  cat("  First 5 rownames:", paste(head(rownames(zscore_data_removedBatchEff_cov), 5), collapse=", "), "\n")
} else if (is.null(rownames(zscore_data_removedBatchEff_cov)) || length(rownames(zscore_data_removedBatchEff_cov)) == 0) {
  cat("  WARNING: rownames are missing and gene_ids_original not available!\n")
  stop("ERROR: Cannot restore rownames. gene_ids_original length:", if(exists("gene_ids_original")) length(gene_ids_original) else "not found", 
       ", data rows:", nrow(zscore_data_removedBatchEff_cov))
} else {
  cat("  Rownames present. Length:", length(rownames(zscore_data_removedBatchEff_cov)), "\n")
  cat("  First 5 rownames:", paste(head(rownames(zscore_data_removedBatchEff_cov), 5), collapse=", "), "\n")
  # Verify they look like gene IDs (Ensembl IDs start with ENSG)
  if (!all(grepl("^ENSG", head(rownames(zscore_data_removedBatchEff_cov), min(10, length(rownames(zscore_data_removedBatchEff_cov))))))) {
    cat("  WARNING: Rownames don't look like Ensembl IDs! Restoring from gene_ids_original...\n")
    if (exists("gene_ids_original") && length(gene_ids_original) == nrow(zscore_data_removedBatchEff_cov)) {
      rownames(zscore_data_removedBatchEff_cov) <- gene_ids_original
      cat("  Rownames restored from original. First 5:", paste(head(rownames(zscore_data_removedBatchEff_cov), 5), collapse=", "), "\n")
    }
  }
}

# Convert to matrix if it's a data.frame (for consistency)
if (is.data.frame(zscore_data_removedBatchEff_cov)) {
  gene_ids_final <- rownames(zscore_data_removedBatchEff_cov)
  zscore_data_removedBatchEff_cov <- as.matrix(zscore_data_removedBatchEff_cov)
  rownames(zscore_data_removedBatchEff_cov) <- gene_ids_final
  cat("  Converted to matrix, rownames preserved.\n")
}

# Save corrected data (save to current directory for Nextflow, also copy to output_dir)
cat("Saving corrected data...\n")
cat("  Final data class:", class(zscore_data_removedBatchEff_cov), "\n")
cat("  Final rownames length:", length(rownames(zscore_data_removedBatchEff_cov)), "\n")

# CHECKPOINT: Final validation before saving corrected_data
validate_gene_ids(rownames(zscore_data_removedBatchEff_cov), "Before saving corrected_data (FINAL CHECKPOINT)")
cat("  ✓ All gene ID validations passed. Saving corrected_data with ", length(rownames(zscore_data_removedBatchEff_cov)), " genes.\n", sep="")
saveRDS(zscore_data_removedBatchEff_cov, opt$corrected_output)

# Save cleaned metadata
write.csv(rosmap_meta_cleaned_rmVar, opt$metadata_output, row.names = FALSE)

# Also copy to output_dir if different from current directory
if (opt$output_dir != "." && opt$output_dir != getwd()) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(opt$corrected_output, file.path(opt$output_dir, opt$corrected_output), overwrite = TRUE)
  file.copy(opt$metadata_output, file.path(opt$output_dir, opt$metadata_output), overwrite = TRUE)
}

# Generate eigencor plot (optional)
cat("Generating eigencor plot...\n")
tryCatch({
  p_removedBatchEff_cov <- PCAtools::pca(zscore_data_removedBatchEff_cov, metadata = metadata)
  eigencor_plot <- eigencorplot(p_removedBatchEff_cov, 
                                metavars = tech_covar_list, 
                                cexLabY = 0.5, 
                                rotLabY = 0.8, 
                                corFUN = "pearson",
                                main = "Correlation of PCs with technical covariates after correction",
                                titleX = "PCs", 
                                titleY = "technical covariates",
                                corMultipleTestCorrection = "hochberg",
                                signifSymbols = c('***', '**', '*', ''),
                                signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                scale = FALSE)
  
  png(file.path(opt$output_dir, "eigencor_plot_removedBatchEff_cov.png"), 
      width = 15, height = 10, units = 'in', res = 300)
  print(eigencor_plot)
  dev.off()
}, error = function(e) {
  cat("Warning: Could not generate eigencor plot:", e$message, "\n")
})

cat("Batch effect and technical covariate removal completed!\n")
cat("Output files saved to:", opt$output_dir, "\n")

