#!/usr/bin/env Rscript
# Prepare GTEx v10 sample metadata for the Nextflow GWAS pipeline
#
# Merges:
#   - GTEx_Analysis_2022-06-06_v10_Annotations_SampleAttributesDS.txt
#     (sample-level: SAMPID, SMTSD tissue, SMRIN RIN score, SMTSISCH ischemic time, etc.)
#   - GTEx_Analysis_2022-06-06_v10_Annotations_SubjectPhenotypesDS.txt
#     (subject-level: SUBJID, SEX, AGE, COHORT, etc.)
#
# Outputs:
#   - GTEx_v10_BA9_sample_metadata.csv  — one row per BA9 sample (SAMPID)
#     Contains: SAMPID, subject_id, SEX, AGE, SMRIN, SMTSISCH, SMNABTCHT, + other QC columns
#
# Usage:
#   Rscript prepare_gtex_v10_metadata.R \
#     --sample_attrs   /path/to/GTEx_..._SampleAttributesDS.txt \
#     --subject_pheno  /path/to/GTEx_..._SubjectPhenotypesDS.txt \
#     --output         /path/to/GTEx_v10_BA9_sample_metadata.csv \
#     [--tissue        "Brain - Frontal Cortex (BA9)"]

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--sample_attrs", type = "character",
    help = "Path to GTEx v10 SampleAttributesDS.txt"),
  make_option("--subject_pheno", type = "character",
    help = "Path to GTEx v10 SubjectPhenotypesDS.txt"),
  make_option("--output", type = "character",
    default = "GTEx_v10_BA9_sample_metadata.csv",
    help = "Output CSV file path [default: GTEx_v10_BA9_sample_metadata.csv]"),
  make_option("--tissue", type = "character",
    default = "Brain - Frontal Cortex (BA9)",
    help = "Tissue to filter on (SMTSD column) [default: 'Brain - Frontal Cortex (BA9)']"),
  make_option("--count_matrix", type = "character", default = NULL,
    help = "Optional: path to count matrix CSV. If provided, metadata rows are keyed on the exact SAMPIDs present as columns in the count matrix (ensures perfect alignment).")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample_attrs) || is.null(opt$subject_pheno)) {
  stop("--sample_attrs and --subject_pheno are required. Run with --help for usage.")
}

cat("=============================================================================\n")
cat("GTEx v10 Metadata Preparation\n")
cat("=============================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load and filter SampleAttributes
# ---------------------------------------------------------------------------
cat("Loading SampleAttributes:", opt$sample_attrs, "\n")
sample_attrs <- fread(opt$sample_attrs, sep = "\t", data.table = FALSE,
                      quote = "")
cat("  Total rows:", nrow(sample_attrs), "\n")
cat("  Columns:", paste(colnames(sample_attrs), collapse = ", "), "\n\n")

# Filter to target tissue
ba9 <- sample_attrs[sample_attrs$SMTSD == opt$tissue, ]
cat("  Rows for '", opt$tissue, "': ", nrow(ba9), "\n", sep = "")

if (nrow(ba9) == 0) {
  stop(paste0("No rows found for tissue '", opt$tissue, "'."))
}

# Extract subject_id from SAMPID: GTEX-XXXXX-YYYY-... -> GTEX-XXXXX
ba9$subject_id <- sapply(strsplit(ba9$SAMPID, "-"),
                          function(x) paste(x[1], x[2], sep = "-"))
cat("  Example subject_ids:", paste(head(unique(ba9$subject_id), 3), collapse = ", "), "\n")
cat("  Unique subjects:", length(unique(ba9$subject_id)), "\n\n")

# ---------------------------------------------------------------------------
# 1b. If count matrix provided, build a SAMPID roster from its column names
#     and use subject-level fallback so metadata rows align exactly with the
#     count matrix columns (same logic as aggregate_reads_per_gene_gtex_v10.R)
# ---------------------------------------------------------------------------
if (!is.null(opt$count_matrix)) {
  cat("Reading count matrix column names from:", opt$count_matrix, "\n")
  cm_header <- scan(opt$count_matrix, what = "", nlines = 1, sep = ",", quiet = TRUE)
  matrix_sampids <- cm_header[cm_header != "gene_id" & cm_header != ""]
  cat("  SAMPIDs in count matrix:", length(matrix_sampids), "\n")

  extract_subject <- function(s) {
    sapply(strsplit(s, "-"), function(x) paste(x[1], x[2], sep = "-"))
  }

  matrix_subjects <- extract_subject(matrix_sampids)

  # For each matrix SAMPID, find the best matching row in ba9 SampleAttributes
  # Prefer exact SAMPID match; fall back to same subject_id + brain tissue code
  cat("  Matching count matrix SAMPIDs to SampleAttributes rows...\n")
  rows_list <- lapply(seq_along(matrix_sampids), function(i) {
    sid  <- matrix_sampids[i]
    subj <- matrix_subjects[i]
    # Exact match
    hit <- ba9[ba9$SAMPID == sid, ]
    if (nrow(hit) > 0) return(hit[1, ])
    # Subject-level fallback: any BA9 row for this subject
    hit <- ba9[ba9$subject_id == subj, ]
    if (nrow(hit) > 0) {
      row <- hit[1, ]
      row$SAMPID <- sid   # use the matrix SAMPID as the key
      return(row)
    }
    # No match — build a minimal skeleton row so the sample is still represented
    # (SEX/AGE will be NA; these samples will be excluded downstream by the pipeline)
    skeleton <- ba9[1, ]
    skeleton[1, ] <- NA
    skeleton$SAMPID     <- sid
    skeleton$subject_id <- subj
    return(skeleton)
  })
  ba9 <- do.call(rbind, rows_list)
  ba9$SAMPID <- matrix_sampids   # guarantee column order matches matrix
  cat("  Rows after count-matrix alignment:", nrow(ba9), "\n")
  cat("  Rows with SEX:", sum(!is.na(ba9$SEX)), "\n\n")
}

# ---------------------------------------------------------------------------
# 2. Load SubjectPhenotypes
# ---------------------------------------------------------------------------
cat("Loading SubjectPhenotypes:", opt$subject_pheno, "\n")
subject_pheno <- fread(opt$subject_pheno, sep = "\t", data.table = FALSE,
                       quote = "")
cat("  Total rows:", nrow(subject_pheno), "\n")
cat("  Columns:", paste(colnames(subject_pheno), collapse = ", "), "\n\n")

# Validate required columns
for (col in c("SUBJID", "SEX", "AGE")) {
  if (!col %in% colnames(subject_pheno)) {
    stop(paste0("Required column '", col, "' not found in SubjectPhenotypes file."))
  }
}

# SEX coding in GTEx: 1 = Male, 2 = Female
# pheno_prep.R handles "1"/"2" coding already via tolower/ifelse logic
cat("  SEX distribution in SubjectPhenotypes:\n")
print(table(subject_pheno$SEX, useNA = "ifany"))
cat("  AGE range:", min(subject_pheno$AGE, na.rm = TRUE), "to",
    max(subject_pheno$AGE, na.rm = TRUE), "\n\n")

# ---------------------------------------------------------------------------
# 3. Merge sample attributes with subject phenotypes
# ---------------------------------------------------------------------------
cat("Merging sample attributes with subject phenotypes...\n")

merged <- left_join(ba9, subject_pheno,
                    by = c("subject_id" = "SUBJID"),
                    suffix = c("", ".subject"))

cat("  Merged rows:", nrow(merged), "\n")
cat("  Samples with SEX:", sum(!is.na(merged$SEX)), "\n")
cat("  Samples with AGE:", sum(!is.na(merged$AGE)), "\n")

# Check for samples without subject phenotype match
n_missing_pheno <- sum(is.na(merged$SEX))
if (n_missing_pheno > 0) {
  cat("  WARNING:", n_missing_pheno, "samples have no subject phenotype match.\n")
  missing_subjects <- unique(merged$subject_id[is.na(merged$SEX)])
  cat("  Missing subject_ids:", paste(head(missing_subjects, 5), collapse = ", "),
      if (length(missing_subjects) > 5) paste0("... (", length(missing_subjects), " total)") else "", "\n")
}

# ---------------------------------------------------------------------------
# 4. Select and rename output columns
# ---------------------------------------------------------------------------
# Core columns always included
core_cols <- c("SAMPID", "subject_id", "SEX", "AGE", "COHORT")

# QC columns from SampleAttributes that are useful as technical covariates
# SMRIN:     RNA Integrity Number
# SMTSISCH:  Ischemic time (seconds)
# SMNABTCHT: Library prep batch
# SMCENTER:  Collection center
# SMGEBTCHT: Genotyping batch
# SMMTRLTP:  Material type
qc_cols <- intersect(
  c("SMRIN", "SMTSISCH", "SMNABTCHT", "SMCENTER", "SMGEBTCHT", "SMMTRLTP",
    "SMNABTCHD", "SMNABTCH", "SMEXPID"),
  colnames(merged)
)

cat("\n  Including QC columns:", paste(qc_cols, collapse = ", "), "\n")

keep_cols <- c(core_cols, qc_cols)
keep_cols <- keep_cols[keep_cols %in% colnames(merged)]

output_df <- merged[, keep_cols, drop = FALSE]

# Remove duplicated SAMPIDs (keep first occurrence)
n_dup <- sum(duplicated(output_df$SAMPID))
if (n_dup > 0) {
  cat("  Removing", n_dup, "duplicate SAMPID rows.\n")
  output_df <- output_df[!duplicated(output_df$SAMPID), ]
}

cat("\n  Final output rows:", nrow(output_df), "\n")
cat("  Final output columns:", paste(colnames(output_df), collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# 5. Write output
# ---------------------------------------------------------------------------
out_dir <- dirname(opt$output)
if (!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}

cat("Writing metadata to:", opt$output, "\n")
fwrite(output_df, opt$output, sep = ",", quote = TRUE, na = "NA")

cat("\n=============================================================================\n")
cat("Metadata preparation complete.\n")
cat("=============================================================================\n")
cat("  File:           ", opt$output, "\n")
cat("  Samples (rows): ", nrow(output_df), "\n")
cat("  Columns:        ", paste(colnames(output_df), collapse = ", "), "\n")
cat("\n")
cat("Pipeline config parameters to use:\n")
cat("  col_sample_id_for_matching = 'SAMPID'\n")
cat("  col_individualID           = 'subject_id'\n")
cat("  col_msex                   = 'SEX'\n")
cat("  col_age                    = 'AGE'\n")
cat("  tissue_filter              = '", opt$tissue, "'\n", sep = "")
cat("=============================================================================\n\n")
