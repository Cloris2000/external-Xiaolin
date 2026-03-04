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
# Format: MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq Effect StdErr P-value Direction
meta_data <- fread(opt$input, header = TRUE, stringsAsFactors = FALSE, showProgress = TRUE) %>%
  as.data.frame() %>%
  separate(MarkerName, into = c("CHROM", "POS", "REF", "ALT"), sep = ":", remove = FALSE) %>%
  rename(
    ID = MarkerName,
    P = `P-value`
  ) %>%
  mutate(
    CHROM = gsub("chr", "", CHROM),  # Remove "chr" prefix if present
    CHROM = as.numeric(CHROM),
    POS = as.numeric(POS),
    P = as.numeric(P)
  ) %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1) %>%
  dplyr::select(ID, CHROM, POS, P)

cat(paste("Loaded", nrow(meta_data), "variants\n"))
cat("Columns:", paste(colnames(meta_data), collapse=", "), "\n")

cat(paste("After filtering:", nrow(meta_data), "valid variants\n"))

if (nrow(meta_data) == 0) {
  stop("No valid variants after filtering!")
}

# Calculate lambda (genomic inflation factor)
cat("Calculating genomic inflation factor (lambda)...\n")
chi2 <- qchisq(1 - meta_data$P, 1)
lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)
cat(paste("Lambda =", round(lambda, 4), "\n"))

# Count significant associations
sig_5e8 <- sum(meta_data$P < 5e-8, na.rm = TRUE)
sig_1e5 <- sum(meta_data$P < 1e-5, na.rm = TRUE)
min_p <- min(meta_data$P, na.rm = TRUE)

cat("\nSummary Statistics:\n")
cat(paste("  Total variants:", nrow(meta_data), "\n"))
cat(paste("  Genome-wide significant (p < 5e-8):", sig_5e8, "\n"))
cat(paste("  Suggestive (p < 1e-5):", sig_1e5, "\n"))
cat(paste("  Minimum p-value:", format(min_p, scientific = TRUE), "\n"))
cat(paste("  Genomic inflation (lambda):", round(lambda, 4), "\n"))

# Generate Manhattan plot with gene annotation using topr
cat("\nGenerating annotated Manhattan plot using topr...\n")
manhattan_annotated_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_manhattan_annotated.png"))

# Create snpset for suggestive hits to highlight
suggestive_snps <- meta_data %>%
  filter(P < 1e-5) %>%
  pull(ID)

cat(paste("  Found", length(suggestive_snps), "suggestive SNPs to highlight\n"))

# Create annotated Manhattan plot with gene annotation
tryCatch({
  p_manhattan <- manhattanExtra(
    df = meta_data,
    genome_wide_thresh = 5e-8,
    suggestive_thresh = 1e-5,
    annotate = 1e-5,  # Annotate SNPs with p < 1e-5
    build = 37,  # GRCh37/hg19
    label_size = 4,
    title = opt$cell_type
  )
  
  # Save with same dimensions as reference script
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

# Prepare data for qqman (needs CHR, BP, SNP, P)
qqman_data <- meta_data %>%
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
qqman::qq(meta_data$P,
   main = paste("Meta-Analysis QQ Plot:", opt$cell_type),
   sub = paste("Lambda =", round(lambda, 3)))
dev.off()
cat(paste("  Saved:", qq_file, "\n"))

# Save summary statistics
summary_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_summary.txt"))
summary_stats <- data.frame(
  CellType = opt$cell_type,
  N_Variants = nrow(meta_data),
  Lambda = round(lambda, 4),
  Min_P = format(min_p, scientific = TRUE),
  Significant_5e8 = sig_5e8,
  Suggestive_1e5 = sig_1e5
)
write.table(summary_stats, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("  Saved:", summary_file, "\n"))

# Find top hits (p < 1e-5)
if (sig_1e5 > 0) {
  top_hits_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_top_hits.txt"))
  top_hits <- meta_data %>%
    filter(P < 1e-5) %>%
    arrange(P) %>%
    head(100)
  write.table(top_hits, file = top_hits_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("  Saved top hits:", top_hits_file, "\n"))
  
  # Also show top 10 in console
  cat("\nTop 10 variants:\n")
  print(head(top_hits, 10))
}

cat("\nVisualization complete!\n")
