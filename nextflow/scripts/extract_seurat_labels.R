#!/usr/bin/env Rscript
#
# Extract per-cell metadata from a Hodge label-transferred Seurat .rds object
# and derive a broad ~19-subclass call (argmax over the broad
# prediction.score.<Subclass> columns), writing a flat per-cell CSV that
# feeds directly into scripts/compute_sn_proportions.R (unmodified).
#
# Background: these *_HodgeLabelled.rds files carry a discrete fine-grained
# cluster label ("predicted.id", ~100-140 categories, e.g. "L5.6.NP_2") plus
# ~19 broad "prediction.score.<Subclass>" *continuous* columns (Astrocyte,
# L4.IT, SST, VIP, ...) matching the resolution already used elsewhere in
# this repo's snRNA-seq GWAS pipeline. There is no ready-made discrete
# subclass column, so we derive one here via row-wise argmax over the broad
# score columns.
#
# Usage:
#   Rscript scripts/extract_seurat_labels.R \
#       --input_rds  data/Brainscope/LabelTransferred/CMC_HodgeLabelled.rds \
#       --output_csv data_input/sn_brainscope_cmc/percell_annotations.csv

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--input_rds", type = "character",
              help = "Path to the *_HodgeLabelled.rds Seurat object"),
  make_option("--subject_col", type = "character", default = "Donor",
              help = "meta.data column identifying the subject/donor [default: %default]"),
  make_option("--output_csv", type = "character",
              help = "Output path for the per-cell annotation CSV"),
  make_option("--metadata_cols", type = "character",
              default = "Cohort,Biological_Sex,Age_death,Disorder,X1000G_ancestry,Genotype_data,PMI,pH,RIN",
              help = "Comma-separated extra meta.data columns to carry through [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds) || is.null(opt$output_csv)) {
  stop("Both --input_rds and --output_csv are required.")
}

cat("=== extract_seurat_labels.R ===\n")
cat("Input RDS  :", opt$input_rds, "\n")
cat("Subject col:", opt$subject_col, "\n")

cat("\nReading Seurat object (this can take a while for large files)...\n")
t0 <- Sys.time()
obj <- readRDS(opt$input_rds)
cat("  Load time:", round(as.numeric(Sys.time() - t0), 1), "s\n")

if (!inherits(obj, "Seurat")) {
  stop("Expected a Seurat object, got class: ", paste(class(obj), collapse = ","))
}

md <- as.data.table(obj@meta.data, keep.rownames = "cell_barcode")
# Free the (potentially huge) expression matrices ASAP.
rm(obj)
invisible(gc())

cat("meta.data:", nrow(md), "cells x", ncol(md), "columns\n")

stopifnot(opt$subject_col %in% colnames(md))
stopifnot("predicted.id" %in% colnames(md))

# ── Identify broad-subclass prediction score columns ─────────────────────────
# Fine-grained cluster columns look like "prediction.score.L5.6.NP_2" or
# "prediction.score.OPC_2_1.SEAAD" (numeric cluster suffix, optionally
# followed by ".SEAAD"). Broad subclass columns have no such suffix, e.g.
# "prediction.score.Astrocyte".
pred_cols  <- grep("^prediction\\.score\\.", colnames(md), value = TRUE)
fine_cols  <- grep("_[0-9]+(\\.SEAAD)?$", pred_cols, value = TRUE)
broad_cols <- setdiff(pred_cols, c(fine_cols, "prediction.score.max"))

if (length(broad_cols) == 0) {
  stop("No broad prediction.score.<Subclass> columns detected -- check schema for ",
       opt$input_rds)
}

cat("\nBroad subclass columns detected (", length(broad_cols), "):\n")
cat(paste(" ", sort(sub("^prediction\\.score\\.", "", broad_cols))), sep = "\n")

# ── Argmax across broad columns to get a discrete call + its score ──────────
score_mat <- as.matrix(md[, ..broad_cols])
storage.mode(score_mat) <- "double"
all_na   <- rowSums(!is.na(score_mat)) == 0
score_mat[is.na(score_mat)] <- -Inf
best_idx <- max.col(score_mat, ties.method = "first")

broad_labels <- sub("^prediction\\.score\\.", "", broad_cols)
md[, broad_subclass       := broad_labels[best_idx]]
md[, broad_subclass_score := score_mat[cbind(seq_len(nrow(score_mat)), best_idx)]]
if (any(all_na)) {
  md[all_na, `:=`(broad_subclass = NA_character_, broad_subclass_score = NA_real_)]
}

cat("\nBroad subclass distribution:\n")
print(md[, .N, by = broad_subclass][order(-N)])

# ── Assemble output columns ──────────────────────────────────────────────────
meta_cols_requested <- trimws(strsplit(opt$metadata_cols, ",")[[1]])
meta_cols_available <- intersect(meta_cols_requested, colnames(md))
missing_meta <- setdiff(meta_cols_requested, meta_cols_available)
if (length(missing_meta) > 0) {
  cat("\nNote: requested metadata columns not found and will be skipped:",
      paste(missing_meta, collapse = ", "), "\n")
}

id_col <- if ("cell_id" %in% colnames(md)) "cell_id" else "cell_barcode"

keep_cols <- unique(c(opt$subject_col, id_col, "predicted.id",
                      "broad_subclass", "broad_subclass_score",
                      "prediction.score.max", meta_cols_available))
keep_cols <- intersect(keep_cols, colnames(md))

out <- md[, ..keep_cols]

dir.create(dirname(opt$output_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(out, opt$output_csv)

cat("\nOutput written:", opt$output_csv, "\n")
cat("Rows:", nrow(out), " Columns:", ncol(out), "\n")
cat("Unique subjects (", opt$subject_col, "):", uniqueN(out[[opt$subject_col]]), "\n")
cat("Done.\n")
