#!/usr/bin/env Rscript
# Manhattan Plot Generator for METAL Meta-Analysis Results
# Generates Manhattan and QQ plots for meta-analysis output with gene annotation

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(qqman)
  library(topr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--input"), type="character", default=NULL,
              help="Input METAL meta-analysis .tbl file", metavar="FILE"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Output prefix for plots", metavar="PREFIX"),
  make_option(c("--cell_type"), type="character", default="Unknown",
              help="Cell type name for plot title", metavar="STRING"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory for plots", metavar="DIR")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate inputs
if (is.null(opt$input)) {
  stop("Input file (--input) is required")
}

if (!file.exists(opt$input)) {
  stop(paste("Input file not found:", opt$input))
}

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Set output prefix
if (is.null(opt$output_prefix)) {
  opt$output_prefix <- gsub("\\.tbl$", "", basename(opt$input))
}

cat("Reading meta-analysis results...\n")
cat(paste("Input file:", opt$input, "\n"))

# Read METAL output using fread (much faster for large files)
# Format with ANALYZE HETEROGENEITY:
#   MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq Effect StdErr
#   P-value Direction HetISq HetChiSq HetDf HetPVal
raw_tbl <- fread(opt$input, header = TRUE, stringsAsFactors = FALSE, showProgress = TRUE) %>%
  as.data.frame()

cat("Columns in METAL output:", paste(colnames(raw_tbl), collapse=", "), "\n")

meta_data <- raw_tbl %>%
  separate(MarkerName, into = c("CHROM", "POS", "REF", "ALT"), sep = ":", remove = FALSE) %>%
  rename(ID = MarkerName, P = `P-value`) %>%
  mutate(
    CHROM = gsub("chr", "", CHROM),
    CHROM = as.numeric(CHROM),
    POS   = as.numeric(POS),
    P     = as.numeric(P)
  ) %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1)

# Carry heterogeneity columns when present (ANALYZE HETEROGENEITY output)
het_cols <- intersect(c("HetISq", "HetChiSq", "HetDf", "HetPVal"), colnames(meta_data))
has_heterogeneity <- length(het_cols) > 0
if (has_heterogeneity) {
  cat("Heterogeneity columns detected:", paste(het_cols, collapse=", "), "\n")
  meta_data <- meta_data %>%
    mutate(across(all_of(het_cols), as.numeric))
} else {
  cat("No heterogeneity columns in this file (expected HetISq, HetChiSq, HetDf, HetPVal).\n")
}

plot_data <- meta_data %>% dplyr::select(ID, CHROM, POS, P)

cat(paste("Loaded", nrow(meta_data), "variants\n"))
cat(paste("After filtering:", nrow(plot_data), "valid variants\n"))

if (nrow(plot_data) == 0) {
  stop("No valid variants after filtering!")
}

# Calculate lambda (genomic inflation factor)
cat("Calculating genomic inflation factor (lambda)...\n")
chi2 <- qchisq(1 - plot_data$P, 1)
lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)
cat(paste("Lambda =", round(lambda, 4), "\n"))

# Count significant associations
sig_5e8 <- sum(plot_data$P < 5e-8, na.rm = TRUE)
sig_1e5 <- sum(plot_data$P < 1e-5, na.rm = TRUE)
min_p <- min(plot_data$P, na.rm = TRUE)

cat("\nSummary Statistics:\n")
cat(paste("  Total variants:", nrow(plot_data), "\n"))
cat(paste("  Genome-wide significant (p < 5e-8):", sig_5e8, "\n"))
cat(paste("  Suggestive (p < 1e-5):", sig_1e5, "\n"))
cat(paste("  Minimum p-value:", format(min_p, scientific = TRUE), "\n"))
cat(paste("  Genomic inflation (lambda):", round(lambda, 4), "\n"))

# Heterogeneity summary (from ANALYZE HETEROGENEITY output)
if (has_heterogeneity) {
  cat("\nHeterogeneity Summary (from METAL ANALYZE HETEROGENEITY):\n")
  het_summary <- meta_data %>%
    filter(!is.na(HetISq)) %>%
    summarise(
      n_variants      = n(),
      median_HetISq   = round(median(HetISq, na.rm = TRUE), 2),
      pct_ISq_gt25    = round(100 * mean(HetISq > 25,  na.rm = TRUE), 1),
      pct_ISq_gt50    = round(100 * mean(HetISq > 50,  na.rm = TRUE), 1),
      pct_ISq_gt75    = round(100 * mean(HetISq > 75,  na.rm = TRUE), 1),
      n_sig_HetPVal   = sum(HetPVal < 0.05, na.rm = TRUE)
    )
  print(het_summary)

  # High-heterogeneity top hits table
  high_het_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_high_het_variants.txt"))
  high_het <- meta_data %>%
    filter(P < 1e-5, !is.na(HetISq), HetISq > 50) %>%
    dplyr::select(ID, CHROM, POS, P, HetISq, HetChiSq, HetDf, HetPVal) %>%
    arrange(P)
  if (nrow(high_het) > 0) {
    write.table(high_het, file = high_het_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(paste("  High-heterogeneity top variants (I² > 50%):", nrow(high_het), "written to", high_het_file, "\n"))
  }
}

# Generate Manhattan plot with gene annotation using topr
cat("\nGenerating annotated Manhattan plot using topr...\n")
manhattan_annotated_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_manhattan_annotated.png"))

# Create snpset for suggestive hits to highlight
suggestive_snps <- plot_data %>%
  filter(P < 1e-5) %>%
  pull(ID)

cat(paste("  Found", length(suggestive_snps), "suggestive SNPs to highlight\n"))

# Create annotated Manhattan plot with gene annotation
tryCatch({
  p_manhattan <- manhattanExtra(
    df = plot_data,
    genome_wide_thresh = 5e-8,
    suggestive_thresh = 1e-5,
    annotate = 1e-5,
    build = 37,  # GRCh37/hg19
    label_size = 4,
    title = opt$cell_type
  )

  ggsave(
    filename = manhattan_annotated_file,
    plot = p_manhattan,
    width = 16,
    height = 6,
    dpi = 300
  )
  cat(paste("  Saved annotated plot:", manhattan_annotated_file, "\n"))
}, error = function(e) {
  cat(paste("  Warning: Could not create annotated plot with topr:", e$message, "\n"))
  cat("  Error details:", toString(e), "\n")
})

# Also generate basic Manhattan plot with qqman (as backup)
cat("Generating basic Manhattan plot with qqman...\n")
manhattan_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_manhattan_basic.png"))

qqman_data <- plot_data %>%
  rename(CHR = CHROM, BP = POS, SNP = ID)

png(file = manhattan_file, width = 4800, height = 1800, res = 300)
qqman::manhattan(qqman_data,
                 suggestiveline = -log10(1e-5),
                 genomewideline = -log10(5e-8),
                 cex = 0.5,
                 cex.axis = 0.8)
title(main = paste("Meta-Analysis Manhattan Plot:", opt$cell_type))
dev.off()
cat(paste("  Saved basic plot:", manhattan_file, "\n"))

# Generate QQ plot
cat("Generating QQ plot...\n")
qq_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_qq.png"))
png(file = qq_file, width = 1400, height = 1400, res = 300)
qqman::qq(plot_data$P,
   main = paste("Meta-Analysis QQ Plot:", opt$cell_type),
   sub = paste("Lambda =", round(lambda, 3)))
dev.off()
cat(paste("  Saved:", qq_file, "\n"))

# Save summary statistics (includes heterogeneity summary when available)
summary_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_summary.txt"))
summary_stats <- data.frame(
  CellType        = opt$cell_type,
  N_Variants      = nrow(plot_data),
  Lambda          = round(lambda, 4),
  Min_P           = format(min_p, scientific = TRUE),
  Significant_5e8 = sig_5e8,
  Suggestive_1e5  = sig_1e5
)
if (has_heterogeneity) {
  summary_stats$Median_HetISq      <- round(median(meta_data$HetISq, na.rm = TRUE), 2)
  summary_stats$Pct_ISq_gt25       <- round(100 * mean(meta_data$HetISq > 25, na.rm = TRUE), 1)
  summary_stats$Pct_ISq_gt50       <- round(100 * mean(meta_data$HetISq > 50, na.rm = TRUE), 1)
  summary_stats$Pct_ISq_gt75       <- round(100 * mean(meta_data$HetISq > 75, na.rm = TRUE), 1)
  summary_stats$N_sig_HetPVal_0.05 <- sum(meta_data$HetPVal < 0.05, na.rm = TRUE)
}
write.table(summary_stats, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("  Saved:", summary_file, "\n"))

# Find top hits (p < 1e-5)
if (sig_1e5 > 0) {
  top_hits_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_top_hits.txt"))
  top_hits_base <- meta_data %>%
    filter(P < 1e-5) %>%
    arrange(P) %>%
    head(100)
  # Attach heterogeneity columns when present
  if (has_heterogeneity) {
    top_hits <- top_hits_base %>%
      dplyr::select(ID, CHROM, POS, P, HetISq, HetChiSq, HetDf, HetPVal)
  } else {
    top_hits <- top_hits_base %>% dplyr::select(ID, CHROM, POS, P)
  }
  write.table(top_hits, file = top_hits_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("  Saved top hits:", top_hits_file, "\n"))
  
  # Also show top 10 in console
  cat("\nTop 10 variants:\n")
  print(head(top_hits, 10))
}

cat("\nVisualization complete!\n")
