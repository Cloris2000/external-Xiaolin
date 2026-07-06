#!/usr/bin/env Rscript
# Convert GTEx v10 GCT gene read counts to pipeline-ready CSV
# filtered to a specific tissue (default: Brain - Frontal Cortex BA9)
#
# GCT format:
#   Line 1: "#1.2"
#   Line 2: <num_genes>\t<num_samples>
#   Line 3: header — Name, Description, SAMPID1, SAMPID2, ...
#   Lines 4+: gene rows
#
# Output: gene_id (rows) x SAMPID (columns), CSV format
# matching the existing pipeline's count_matrix_file expectation

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--gct", type = "character",
    help = "Path to GTEx GCT gz file (gene_reads.gct.gz)"),
  make_option("--sample_attrs", type = "character",
    help = "Path to GTEx SampleAttributesDS.txt"),
  make_option("--output", type = "character",
    default = "GTEx_v10_BA9_raw_count_matrix.csv",
    help = "Output CSV path [default: GTEx_v10_BA9_raw_count_matrix.csv]"),
  make_option("--tissue", type = "character",
    default = "Brain - Frontal Cortex (BA9)",
    help = "Tissue to filter on (SMTSD) [default: Brain - Frontal Cortex (BA9)]"),
  make_option("--strip_version", type = "logical", default = TRUE,
    help = "Strip Ensembl version suffix from gene IDs [default: TRUE]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$gct) || is.null(opt$sample_attrs)) {
  stop("--gct and --sample_attrs are required.")
}

cat("=============================================================================\n")
cat("GTEx v10 GCT -> Pipeline CSV Conversion\n")
cat("=============================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Get target SAMPIDs from SampleAttributes
# ---------------------------------------------------------------------------
cat("Loading SampleAttributes:", opt$sample_attrs, "\n")
sample_attrs <- fread(opt$sample_attrs, sep = "\t", data.table = FALSE, quote = "")
ba9 <- sample_attrs[sample_attrs$SMTSD == opt$tissue, "SAMPID"]
cat("  BA9 samples in SampleAttributes:", length(ba9), "\n\n")

# ---------------------------------------------------------------------------
# 2. Read GCT header to find which columns are BA9
# ---------------------------------------------------------------------------
cat("Reading GCT header:", opt$gct, "\n")
con <- gzcon(file(opt$gct, "rb"))
header_lines <- readLines(con, n = 3)
close(con)

# Line 3 is the column header
col_names <- strsplit(header_lines[3], "\t")[[1]]
cat("  Total columns in GCT:", length(col_names), "\n")

# First two columns are Name and Description
sample_cols <- col_names[3:length(col_names)]
cat("  Total samples in GCT:", length(sample_cols), "\n")

# Find which sample columns are BA9
ba9_cols <- intersect(sample_cols, ba9)
cat("  BA9 samples found in GCT:", length(ba9_cols), "\n")

missing <- setdiff(ba9, sample_cols)
if (length(missing) > 0) {
  cat("  WARNING:", length(missing), "BA9 SAMPIDs from SampleAttributes not in GCT\n")
  cat("  First 5 missing:", paste(head(missing, 5), collapse = ", "), "\n")
}

# Column indices to keep (1-based): Name col + BA9 sample cols
name_idx <- 1  # "Name" = gene_id
ba9_idx <- which(col_names %in% ba9_cols)
keep_idx <- c(name_idx, ba9_idx)
cat("  Keeping", length(keep_idx), "columns (1 gene_id + ", length(ba9_idx), "samples)\n\n")

# ---------------------------------------------------------------------------
# 3. Read the full GCT, selecting only needed columns
# ---------------------------------------------------------------------------
cat("Reading GCT data (this may take a minute for a ~900MB file)...\n")

# fread with select= is very efficient — reads only needed columns
gct_data <- fread(
  cmd = paste("zcat", shQuote(opt$gct), "| tail -n +3"),
  sep = "\t",
  header = TRUE,
  select = keep_idx,
  data.table = FALSE,
  showProgress = TRUE
)

cat("  Read", nrow(gct_data), "genes x", ncol(gct_data) - 1, "samples\n")

# ---------------------------------------------------------------------------
# 4. Clean up gene IDs
# ---------------------------------------------------------------------------
# Rename first column to gene_id
colnames(gct_data)[1] <- "gene_id"

# Strip Ensembl version suffix (e.g. ENSG00000223972.5 -> ENSG00000223972)
if (opt$strip_version) {
  gct_data$gene_id <- sub("\\.[0-9]+$", "", gct_data$gene_id)
  cat("  Stripped Ensembl version suffixes from gene IDs\n")
}

# Remove duplicate gene IDs (keep first occurrence)
n_dup <- sum(duplicated(gct_data$gene_id))
if (n_dup > 0) {
  cat("  Removing", n_dup, "duplicate gene IDs\n")
  gct_data <- gct_data[!duplicated(gct_data$gene_id), ]
}

cat("  Final matrix:", nrow(gct_data), "genes x", ncol(gct_data) - 1, "samples\n\n")

# ---------------------------------------------------------------------------
# 5. Write output
# ---------------------------------------------------------------------------
out_dir <- dirname(opt$output)
if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)

cat("Writing output to:", opt$output, "\n")
fwrite(gct_data, opt$output, sep = ",", quote = FALSE, na = "NA")

# Quick sanity check
sample_totals <- colSums(gct_data[, -1], na.rm = TRUE)
cat("\n=============================================================================\n")
cat("Done.\n")
cat("  File:    ", opt$output, "\n")
cat("  Genes:   ", nrow(gct_data), "\n")
cat("  Samples: ", ncol(gct_data) - 1, "\n")
cat("  Tissue:  ", opt$tissue, "\n")
cat("  Count range per sample: [", min(sample_totals), ",", max(sample_totals), "]\n")
cat("=============================================================================\n\n")
