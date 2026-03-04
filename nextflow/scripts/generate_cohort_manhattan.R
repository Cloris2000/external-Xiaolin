#!/usr/bin/env Rscript
# Manhattan Plot Generator for Per-Cohort REGENIE raw_p Results
# Reads .regenie.raw_p format - generates annotated Manhattan plot (suggestive SNPs with genes)

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(qqman)
  library(topr)
})

option_list <- list(
  make_option(c("--input"), type="character", default=NULL,
              help="Input .regenie.raw_p file", metavar="FILE"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Output prefix for plots", metavar="PREFIX"),
  make_option(c("--cell_type"), type="character", default="Unknown",
              help="Cell type name for plot title", metavar="STRING"),
  make_option(c("--cohort"), type="character", default=NULL,
              help="Cohort name for plot title", metavar="STRING"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory for plots", metavar="DIR"),
  make_option(c("--downsample_n"), type="integer", default=100000,
              help="Max SNPs to sample from non-suggestive (P>=1e-5) for plotting", metavar="N")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) stop("Input file (--input) is required")
if (!file.exists(opt$input)) stop(paste("Input file not found:", opt$input))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
if (is.null(opt$output_prefix)) opt$output_prefix <- gsub("\\.regenie\\.raw_p$", "", basename(opt$input))

title_label <- if (!is.null(opt$cohort)) paste0(opt$cell_type, " (", opt$cohort, ")") else opt$cell_type

cat("Reading REGENIE raw_p results...\n")
cat(paste("Input file:", opt$input, "\n"))

# Read raw_p format: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA P
gwas_data <- fread(opt$input, header = TRUE, stringsAsFactors = FALSE, showProgress = TRUE) %>%
  as.data.frame() %>%
  rename(POS = GENPOS) %>%
  mutate(
    CHROM = as.numeric(gsub("chr", "", CHROM)),
    POS = as.numeric(POS),
    P = as.numeric(P)
  ) %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1) %>%
  dplyr::select(ID, CHROM, POS, P)

cat(paste("Loaded", nrow(gwas_data), "variants\n"))
if (nrow(gwas_data) == 0) stop("No valid variants after filtering!")

# Lambda
chi2 <- qchisq(1 - gwas_data$P, 1)
lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)
sig_5e8 <- sum(gwas_data$P < 5e-8, na.rm = TRUE)
sig_1e5 <- sum(gwas_data$P < 1e-5, na.rm = TRUE)
min_p <- min(gwas_data$P, na.rm = TRUE)

cat("\nSummary Statistics:\n")
cat(paste("  Total variants:", nrow(gwas_data), "\n"))
cat(paste("  Genome-wide significant (p < 5e-8):", sig_5e8, "\n"))
cat(paste("  Suggestive (p < 1e-5):", sig_1e5, "\n"))
cat(paste("  Minimum p-value:", format(min_p, scientific = TRUE), "\n"))
cat(paste("  Lambda:", round(lambda, 4), "\n"))

# Downsample for plotting: keep all suggestive (P < 1e-5), sample from rest
important <- gwas_data %>% filter(P < 1e-5)
other <- gwas_data %>% filter(P >= 1e-5)
if (nrow(other) > opt$downsample_n) {
  set.seed(42)
  other_sampled <- other %>% sample_n(opt$downsample_n)
  gwas_plot <- bind_rows(important, other_sampled) %>% arrange(CHROM, POS)
  cat(paste("\nDownsampled for plotting:", nrow(gwas_data), "->", nrow(gwas_plot),
            "variants (kept all", nrow(important), "suggestive, sampled", opt$downsample_n, "from rest)\n"))
} else {
  gwas_plot <- gwas_data
  cat(paste("\nNo downsampling needed (", nrow(gwas_data), "variants)\n"))
}

# Annotated Manhattan (topr) - only suggestive SNPs annotated, use downsampled data
manhattan_annotated_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_manhattan_annotated.png"))
tryCatch({
  p_manhattan <- manhattanExtra(
    df = gwas_plot,
    genome_wide_thresh = 5e-8,
    suggestive_thresh = 1e-5,
    annotate = 1e-5,
    build = 37,
    label_size = 4,
    title = title_label
  )
  ggsave(manhattan_annotated_file, plot = p_manhattan, width = 16, height = 6, dpi = 300)
  cat(paste("  Saved annotated plot:", manhattan_annotated_file, "\n"))
}, error = function(e) {
  cat(paste("  Warning: Could not create annotated plot:", e$message, "\n"))
})

# # Basic Manhattan (qqman)
# manhattan_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_manhattan_basic.png"))
# qqman_data <- gwas_plot %>% rename(CHR = CHROM, BP = POS, SNP = ID)
# png(file = manhattan_file, width = 4800, height = 1800, res = 300)
# qqman::manhattan(qqman_data,
#                 suggestiveline = -log10(1e-5),
#                 genomewideline = -log10(5e-8),
#                 cex = 0.5,
#                 cex.axis = 0.8)
# title(main = paste("Manhattan Plot:", title_label))
# dev.off()
# cat(paste("  Saved basic plot:", manhattan_file, "\n"))

# # QQ plot
# qq_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_qq.png"))
# png(file = qq_file, width = 1400, height = 1400, res = 300)
# qqman::qq(gwas_data$P, main = paste("QQ Plot:", title_label), sub = paste("Lambda =", round(lambda, 3)))
# dev.off()
# cat(paste("  Saved:", qq_file, "\n"))

# # Summary
# summary_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_summary.txt"))
# write.table(data.frame(
#   Cohort = if (!is.null(opt$cohort)) opt$cohort else ".",
#   CellType = opt$cell_type,
#   N_Variants = nrow(gwas_data),
#   N_Plotted = nrow(gwas_plot),
#   Lambda = round(lambda, 4),
#   Min_P = format(min_p, scientific = TRUE),
#   Significant_5e8 = sig_5e8,
#   Suggestive_1e5 = sig_1e5
# ), file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
# cat(paste("  Saved:", summary_file, "\n"))

# if (sig_1e5 > 0) {
#   top_hits_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_top_hits.txt"))
#   top_hits <- gwas_data %>% filter(P < 1e-5) %>% arrange(P) %>% head(100)
#   write.table(top_hits, file = top_hits_file, sep = "\t", quote = FALSE, row.names = FALSE)
#   cat(paste("  Saved top hits:", top_hits_file, "\n"))
# }

cat("\nDone.\n")
