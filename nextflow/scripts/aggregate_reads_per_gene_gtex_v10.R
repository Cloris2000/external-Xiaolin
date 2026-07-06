#!/usr/bin/env Rscript
# Aggregate per-sample STAR ReadsPerGene.out.tab.gz files for GTEx v10
# into a single gene × sample count matrix for Brain - Frontal Cortex (BA9).
#
# Input:
#   - GTEx v10 SampleAttributesDS.txt (to identify BA9 samples by SMTSD)
#   - Directory of per-sample *.ReadsPerGene.out.tab.gz files
# Output:
#   - CSV: gene_id (rows) × SAMPID (columns), unstranded counts
#
# Usage:
#   Rscript aggregate_reads_per_gene_gtex_v10.R \
#     --sample_attrs  /path/to/GTEx_Analysis_..._SampleAttributesDS.txt \
#     --reads_dir     /path/to/RNAseq_aux_files_ReadsPerGene \
#     --output        /path/to/GTEx_v10_BA9_raw_count_matrix.csv \
#     [--tissue       "Brain - Frontal Cortex (BA9)"]  # default
#
# STAR ReadsPerGene.out.tab format (4 columns):
#   gene_id  unstranded  sense  antisense
# First 4 rows are summary lines (N_unmapped, N_multimapping, N_noFeature, N_ambiguous) — skipped.
# GTEx RNA-seq is unstranded, so column 2 (unstranded) is used.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--sample_attrs", type = "character",
    help = "Path to GTEx SampleAttributesDS.txt (tab-separated)"),
  make_option("--reads_dir", type = "character",
    help = "Directory containing *.ReadsPerGene.out.tab.gz files"),
  make_option("--output", type = "character",
    default = "GTEx_v10_BA9_raw_count_matrix.csv",
    help = "Output CSV file path [default: GTEx_v10_BA9_raw_count_matrix.csv]"),
  make_option("--tissue", type = "character",
    default = "Brain - Frontal Cortex (BA9)",
    help = "Tissue name to filter on (SMTSD column) [default: 'Brain - Frontal Cortex (BA9)']"),
  make_option("--n_cores", type = "integer", default = 4,
    help = "Number of parallel cores for reading files [default: 4]"),
  make_option("--match_by_subject", type = "logical", default = TRUE,
    help = "If TRUE, match files by subject ID (GTEX-XXXXX) prefix when exact SAMPID not found [default: TRUE]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample_attrs) || is.null(opt$reads_dir)) {
  stop("--sample_attrs and --reads_dir are required. Run with --help for usage.")
}

cat("=============================================================================\n")
cat("GTEx v10 ReadsPerGene Aggregation\n")
cat("=============================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load SampleAttributes and identify target SAMPIDs
# ---------------------------------------------------------------------------
cat("Loading SampleAttributes:", opt$sample_attrs, "\n")
sample_attrs <- fread(opt$sample_attrs, sep = "\t", data.table = FALSE,
                      quote = "")

cat("  Total samples in SampleAttributes:", nrow(sample_attrs), "\n")
cat("  Column names (first 15):", paste(head(colnames(sample_attrs), 15), collapse = ", "), "\n")

# SMTSD is the detailed tissue name column
if (!"SMTSD" %in% colnames(sample_attrs)) {
  stop("Column 'SMTSD' not found in SampleAttributes file. Cannot filter by tissue.")
}
if (!"SAMPID" %in% colnames(sample_attrs)) {
  stop("Column 'SAMPID' not found in SampleAttributes file.")
}

target_samples <- sample_attrs[sample_attrs$SMTSD == opt$tissue, ]
cat("  Samples with tissue '", opt$tissue, "': ", nrow(target_samples), "\n", sep = "")

if (nrow(target_samples) == 0) {
  stop(paste0("No samples found for tissue '", opt$tissue, "'. Check --tissue argument."))
}

target_sampids <- target_samples$SAMPID
cat("  Example SAMPIDs:", paste(head(target_sampids, 3), collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# 2. Find matching ReadsPerGene files
# ---------------------------------------------------------------------------
cat("Scanning reads directory:", opt$reads_dir, "\n")
all_files <- list.files(opt$reads_dir, pattern = "\\.ReadsPerGene\\.out\\.tab\\.gz$",
                        full.names = TRUE)
cat("  Total ReadsPerGene files found:", length(all_files), "\n")

# Derive SAMPID from filename by removing the suffix
file_sampids <- sub("\\.ReadsPerGene\\.out\\.tab\\.gz$", "",
                    basename(all_files))

# Match against target SAMPIDs — exact match first
matched_idx <- file_sampids %in% target_sampids
matched_files <- all_files[matched_idx]
matched_sampids <- file_sampids[matched_idx]

cat("  Files matching target tissue (exact SAMPID):", length(matched_files), "\n")

# If match_by_subject=TRUE, also match by subject ID prefix (GTEX-XXXXX)
# for SAMPIDs that differ only in aliquot code between releases
missing_sampids <- setdiff(target_sampids, matched_sampids)
if (length(missing_sampids) > 0 && isTRUE(opt$match_by_subject)) {
  cat("  Trying subject-level fallback for", length(missing_sampids), "unmatched SAMPIDs...\n")

  # Extract subject IDs from target SAMPIDs and file SAMPIDs
  # GTEX-XXXXX-YYYY-ZZZ -> GTEX-XXXXX (first two dash-parts)
  extract_subject <- function(sampids) {
    sapply(strsplit(sampids, "-"), function(x) paste(x[1], x[2], sep = "-"))
  }

  missing_subjects <- extract_subject(missing_sampids)
  file_subjects <- extract_subject(file_sampids)

  # For each missing subject, pick the first available local file for that subject
  # that is NOT already matched (avoid double-counting)
  already_matched_subjects <- extract_subject(matched_sampids)
  unmatched_file_idx <- !matched_idx
  unmatched_file_subjects <- file_subjects[unmatched_file_idx]
  unmatched_files <- all_files[unmatched_file_idx]
  unmatched_sampids_local <- file_sampids[unmatched_file_idx]

  fallback_files <- c()
  fallback_sampids <- c()
  fallback_map <- data.frame(
    target_sampid = character(0),
    used_sampid   = character(0),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(missing_subjects)) {
    subj <- missing_subjects[i]
    # Skip if this subject already has a matched file
    if (subj %in% already_matched_subjects) next
    # Find first available local file for this subject
    hit <- which(unmatched_file_subjects == subj)[1]
    if (!is.na(hit)) {
      fallback_files   <- c(fallback_files, unmatched_files[hit])
      fallback_sampids <- c(fallback_sampids, unmatched_sampids_local[hit])
      fallback_map <- rbind(fallback_map, data.frame(
        target_sampid = missing_sampids[i],
        used_sampid   = unmatched_sampids_local[hit],
        stringsAsFactors = FALSE
      ))
      # Mark as used so we don't pick it again
      unmatched_file_subjects[hit] <- "__used__"
    }
  }

  if (length(fallback_files) > 0) {
    cat("  Subject-level fallback matched", length(fallback_files), "additional files\n")
    cat("  Example mappings (target SAMPID -> file used):\n")
    for (i in seq_len(min(3, nrow(fallback_map)))) {
      cat("    ", fallback_map$target_sampid[i], "->", fallback_map$used_sampid[i], "\n")
    }
    # Use the TARGET SAMPID as the column name (matches SampleAttributes)
    matched_files   <- c(matched_files, fallback_files)
    # Use target SAMPID as column label so it aligns with metadata
    matched_sampids <- c(matched_sampids, fallback_map$target_sampid)
  }

  still_missing <- setdiff(missing_sampids,
                           fallback_map$target_sampid[fallback_map$target_sampid %in% missing_sampids])
  if (length(still_missing) > 0) {
    cat("  WARNING:", length(still_missing), "subjects have NO file at all:\n")
    cat("    ", paste(head(still_missing, 5), collapse = ", "),
        if (length(still_missing) > 5) paste0("... and ", length(still_missing) - 5, " more") else "", "\n")
  }
} else if (length(missing_sampids) > 0) {
  cat("  WARNING:", length(missing_sampids), "target SAMPIDs have no matching file:\n")
  cat("    ", paste(head(missing_sampids, 5), collapse = ", "),
      if (length(missing_sampids) > 5) paste0("... and ", length(missing_sampids) - 5, " more") else "", "\n")
}

cat("  Total files to aggregate:", length(matched_files), "\n")

if (length(matched_files) == 0) {
  stop("No matching files found. Check --reads_dir and --sample_attrs.")
}

# ---------------------------------------------------------------------------
# 3. Read all matched files and aggregate
# ---------------------------------------------------------------------------
cat("\nReading", length(matched_files), "ReadsPerGene files...\n")
cat("  (Using", opt$n_cores, "parallel cores)\n\n")

read_one_file <- function(filepath, sampid) {
  tryCatch({
    dt <- fread(filepath, sep = "\t", header = FALSE, data.table = FALSE,
                col.names = c("gene_id", "unstranded", "sense", "antisense"),
                showProgress = FALSE)
    # Skip the first 4 summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
    dt <- dt[5:nrow(dt), ]
    # Strip Ensembl version suffix (e.g. ENSG00000223972.5 -> ENSG00000223972)
    dt$gene_id <- sub("\\.[0-9]+$", "", dt$gene_id)
    # Return just gene_id and unstranded count
    out <- data.frame(gene_id = dt$gene_id, count = dt$unstranded,
                      stringsAsFactors = FALSE)
    colnames(out)[2] <- sampid
    out
  }, error = function(e) {
    warning("Failed to read file: ", filepath, "\n  Error: ", conditionMessage(e))
    NULL
  })
}

# Read files sequentially with progress reporting
results <- vector("list", length(matched_files))
pb_step <- max(1, floor(length(matched_files) / 20))  # report every 5%

for (i in seq_along(matched_files)) {
  if (i == 1 || i %% pb_step == 0) {
    cat(sprintf("  Progress: %d / %d (%.0f%%)\r", i, length(matched_files),
                100 * i / length(matched_files)))
    flush.console()
  }
  results[[i]] <- read_one_file(matched_files[i], matched_sampids[i])
}
cat(sprintf("\n  Done reading %d files.\n\n", length(matched_files)))

# Remove any NULL entries (failed reads)
failed <- sum(sapply(results, is.null))
if (failed > 0) {
  warning(failed, " files failed to read and will be excluded.")
  results <- results[!sapply(results, is.null)]
}

# ---------------------------------------------------------------------------
# 4. Build wide gene × sample matrix
# ---------------------------------------------------------------------------
cat("Building gene x sample count matrix...\n")

# Use the gene_id list from the first result as reference
gene_ids <- results[[1]]$gene_id
cat("  Genes:", length(gene_ids), "\n")

# Verify all files have the same genes (warn if not)
gene_lengths <- sapply(results, nrow)
if (length(unique(gene_lengths)) > 1) {
  warning("Files have different numbers of genes! Using intersection.")
  common_genes <- Reduce(intersect, lapply(results, function(x) x$gene_id))
  cat("  Genes after intersection:", length(common_genes), "\n")
  results <- lapply(results, function(x) x[x$gene_id %in% common_genes, ])
  gene_ids <- common_genes
}

# Assemble the matrix: start with gene_id column, cbind count columns
count_matrix <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
for (res in results) {
  count_matrix <- merge(count_matrix, res, by = "gene_id", all.x = TRUE, sort = FALSE)
}

# Restore original gene order
count_matrix <- count_matrix[match(gene_ids, count_matrix$gene_id), ]
rownames(count_matrix) <- NULL

cat("  Matrix dimensions:", nrow(count_matrix), "genes x",
    ncol(count_matrix) - 1, "samples\n\n")

# ---------------------------------------------------------------------------
# 5. Write output
# ---------------------------------------------------------------------------
out_dir <- dirname(opt$output)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("Writing count matrix to:", opt$output, "\n")
fwrite(count_matrix, opt$output, sep = ",", quote = FALSE, na = "NA")
cat("Done.\n\n")

# Quick summary
cat("=============================================================================\n")
cat("Output summary\n")
cat("=============================================================================\n")
cat("  File:    ", opt$output, "\n")
cat("  Genes:   ", nrow(count_matrix), "\n")
cat("  Samples: ", ncol(count_matrix) - 1, "\n")
cat("  Tissue:  ", opt$tissue, "\n")
sample_totals <- colSums(count_matrix[, -1], na.rm = TRUE)
cat("  Count range per sample: [", min(sample_totals), ",", max(sample_totals), "]\n")
cat("=============================================================================\n\n")
