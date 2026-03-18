#!/usr/bin/env Rscript
meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data_input/nimh_hbcc/HBCC_metadata.csv")
# Step 1: drop cols with any NA
cols_without_na <- complete.cases(t(meta))
dropped_na <- setdiff(colnames(meta), colnames(meta)[cols_without_na])
meta_cleaned <- meta[, cols_without_na, drop = FALSE]
# Step 2: drop constant cols (unique==1)
num_unique <- sapply(meta_cleaned, function(x) length(unique(x)))
cols_to_keep <- num_unique > 1
dropped_const <- colnames(meta_cleaned)[!cols_to_keep]
meta_rmVar <- meta_cleaned[, cols_to_keep, drop = FALSE]
# Step 3: exclude IDs/non-covariates
exclude <- c("X", "specimenID", "id", "synapseID", "tissue", "assay", "organ",
             "Started.job.on", "Started.mapping.on", "Finished.on")
tech_cov <- setdiff(colnames(meta_rmVar), exclude)
cat("Tech covariates (HBCC):", length(tech_cov), "\n")
cat(paste(tech_cov, collapse = ", "), "\n")
cat("\nDropped (any NA):", paste(dropped_na, collapse = ", "), "\n")
cat("Dropped (constant):", paste(dropped_const, collapse = ", "), "\n")
