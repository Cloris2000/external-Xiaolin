#!/usr/bin/env Rscript
#
# Compute per-subject cell type proportions from single-nucleus annotation data.
# Outputs two files consumed by the combined pipeline in precomputed-proportions mode:
#   1. cell_proportions.csv  – wide-format: one row per subject, columns = cell types
#   2. metadata.csv          – one row per subject with clinical/demographic columns

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--input_file", type = "character",
              help = "Path to the per-cell annotation CSV (e.g. Green_sn.csv)"),
  make_option("--subject_col", type = "character", default = "Donor.ID",
              help = "Column identifying the subject/donor [default: %default]"),
  make_option("--celltype_col", type = "character", default = "Subclass",
              help = "Column with cell type annotations [default: %default]"),
  make_option("--confidence_col", type = "character", default = NULL,
              help = "Optional confidence score column for filtering"),
  make_option("--confidence_threshold", type = "double", default = 0.8,
              help = "Minimum confidence to retain a cell [default: %default]"),
  make_option("--source_filter", type = "character", default = NULL,
              help = "Optional: keep only rows where 'Source' == this value"),
  make_option("--crosswalk_file", type = "character", default = NULL,
              help = "Optional CSV crosswalk with columns 'from_id' and 'to_id' to remap subject IDs after aggregation"),
  make_option("--crosswalk_from", type = "character", default = "from_id",
              help = "Column name in crosswalk file for the source ID (matches individualID) [default: %default]"),
  make_option("--crosswalk_to", type = "character", default = "to_id",
              help = "Column name in crosswalk file for the target ID (e.g. CMC-style) [default: %default]"),
  make_option("--metadata_cols", type = "character",
              default = "msex,age_death,Study,projid",
              help = "Comma-separated metadata columns to extract per subject [default: %default]"),
  make_option("--output_proportions", type = "character",
              default = "cell_proportions.csv",
              help = "Output path for proportions CSV [default: %default]"),
  make_option("--output_metadata", type = "character",
              default = "metadata.csv",
              help = "Output path for metadata CSV [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("=== compute_sn_proportions.R ===\n")
cat("Input file :", opt$input_file, "\n")
cat("Subject col:", opt$subject_col, "\n")
cat("Cell type  :", opt$celltype_col, "\n")

# ── Read input ───────────────────────────────────────────────────────────────
cat("\nReading input file (this may take a moment for large files)...\n")
dt <- fread(opt$input_file, header = TRUE)
cat("  Rows read:", nrow(dt), " Columns:", ncol(dt), "\n")

# Drop unnamed row-index column produced by R's write.csv(row.names=TRUE)
if (colnames(dt)[1] %in% c("", "V1") && is.integer(dt[[1]])) {
  cat("  Dropping unnamed row-index column:", colnames(dt)[1], "\n")
  dt[, (colnames(dt)[1]) := NULL]
}

stopifnot(opt$subject_col %in% colnames(dt))
stopifnot(opt$celltype_col %in% colnames(dt))

# ── Filter by source if requested ────────────────────────────────────────────
if (!is.null(opt$source_filter) && "Source" %in% colnames(dt)) {
  n_before <- nrow(dt)
  dt <- dt[Source == opt$source_filter]
  cat("  Source filter '", opt$source_filter, "': kept", nrow(dt), "of", n_before, "rows\n")
}

# ── Filter empty cell type labels ────────────────────────────────────────────
empty_ct <- dt[[opt$celltype_col]] == "" | is.na(dt[[opt$celltype_col]])
if (any(empty_ct)) {
  cat("  Dropping", sum(empty_ct), "cells with empty cell type\n")
  dt <- dt[!empty_ct]
}

# ── Filter by confidence ─────────────────────────────────────────────────────
if (!is.null(opt$confidence_col) && opt$confidence_col %in% colnames(dt)) {
  n_before <- nrow(dt)
  dt <- dt[as.numeric(dt[[opt$confidence_col]]) >= opt$confidence_threshold]
  cat("  Confidence >=", opt$confidence_threshold, ": kept", nrow(dt), "of", n_before, "cells\n")
}

# ── Clean cell type names (make safe for column names and filenames) ─────────
# Replace /, spaces, and hyphens with dots
clean_name <- function(x) gsub("[ /\\-]+", ".", x)
dt[, ct_clean := clean_name(get(opt$celltype_col))]

cat("\nUnique cell types (", uniqueN(dt$ct_clean), "):\n")
cat(paste(" ", sort(unique(dt$ct_clean))), sep = "\n")

# ── Count cells per subject x cell type ──────────────────────────────────────
counts <- dt[, .N, by = c(opt$subject_col, "ct_clean")]

# Compute proportions within each subject
counts[, total := sum(N), by = eval(opt$subject_col)]
counts[, proportion := N / total]

cat("\nSubjects:", uniqueN(counts[[opt$subject_col]]), "\n")

# ── Pivot to wide format ─────────────────────────────────────────────────────
props_wide <- dcast(counts, get(opt$subject_col) ~ ct_clean, value.var = "proportion", fill = 0)
setnames(props_wide, old = names(props_wide)[1], new = "individualID")

cat("Proportions matrix:", nrow(props_wide), "subjects x",
    ncol(props_wide) - 1, "cell types\n")

# ── Extract per-subject metadata ─────────────────────────────────────────────
meta_cols_requested <- trimws(strsplit(opt$metadata_cols, ",")[[1]])
meta_cols_available <- intersect(meta_cols_requested, colnames(dt))
cat("\nMetadata columns requested:", paste(meta_cols_requested, collapse = ", "), "\n")
cat("Metadata columns found    :", paste(meta_cols_available, collapse = ", "), "\n")

if (length(meta_cols_available) > 0) {
  meta <- unique(dt[, c(opt$subject_col, meta_cols_available), with = FALSE])
} else {
  meta <- unique(dt[, opt$subject_col, with = FALSE])
}
setnames(meta, old = opt$subject_col, new = "individualID")
# Deduplicate: keep first row per subject
meta <- meta[!duplicated(individualID)]

cat("Metadata rows:", nrow(meta), "\n")

# ── Optional crosswalk: remap individualID to a different ID system ──────────
if (!is.null(opt$crosswalk_file)) {
  cat("\nApplying crosswalk from:", opt$crosswalk_file, "\n")
  xwalk <- fread(opt$crosswalk_file, header = TRUE)
  from_col <- opt$crosswalk_from
  to_col   <- opt$crosswalk_to
  stopifnot(from_col %in% colnames(xwalk), to_col %in% colnames(xwalk))
  xwalk_map <- setNames(xwalk[[to_col]], xwalk[[from_col]])

  n_before <- nrow(props_wide)
  props_wide[, individualID := xwalk_map[individualID]]
  props_wide <- props_wide[!is.na(individualID)]
  cat("  Proportions: remapped", nrow(props_wide), "of", n_before, "subjects\n")

  meta[, individualID := xwalk_map[individualID]]
  meta <- meta[!is.na(individualID)]
  cat("  Metadata:    remapped", nrow(meta), "subjects\n")
}

# ── Write outputs ────────────────────────────────────────────────────────────
dir.create(dirname(opt$output_proportions), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$output_metadata), showWarnings = FALSE, recursive = TRUE)

fwrite(props_wide, opt$output_proportions)
fwrite(meta, opt$output_metadata)

cat("\nOutputs written:\n")
cat("  Proportions:", opt$output_proportions, "\n")
cat("  Metadata   :", opt$output_metadata, "\n")
cat("Done.\n")
