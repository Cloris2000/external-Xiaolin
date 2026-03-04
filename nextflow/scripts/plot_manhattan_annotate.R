#!/usr/bin/env Rscript
# Manhattan plot and gene annotation for GWAS results
# 
# MODE 1 - Generate Manhattan plot:
#   Usage: Rscript plot_manhattan_annotate.R <gwas_file> <output_prefix> <cell_type> [full]
#   Add 'full' as 4th argument to plot all variants (no downsampling, slower)
#
# MODE 2 - Combine two Manhattan plots vertically:
#   Usage: Rscript plot_manhattan_annotate.R combine <plot1.png> <plot2.png> <output.png> [label1] [label2]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(qqman)
  library(topr)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if in COMBINE mode
if (length(args) >= 1 && tolower(args[1]) == "combine") {
  # ============================================================
  # COMBINE MODE - Stack two Manhattan plots vertically
  # ============================================================
  
  if (length(args) < 4) {
    stop("COMBINE mode usage: Rscript plot_manhattan_annotate.R combine <plot1.png> <plot2.png> <output.png> [label1] [label2]")
  }
  
  # No R packages needed - will use ImageMagick command line tool
  
  plot1_file <- args[2]
  plot2_file <- args[3]
  output_file <- args[4]
  label1 <- if (length(args) >= 5) args[5] else "A"
  label2 <- if (length(args) >= 6) args[6] else "B"
  
  cat("============================================================\n")
  cat("COMBINE MODE: Stacking Manhattan Plots Vertically\n")
  cat("============================================================\n")
  cat("Top panel:   ", plot1_file, "\n")
  cat("Bottom panel:", plot2_file, "\n")
  cat("Output:      ", output_file, "\n")
  cat("Labels:      ", label1, "|", label2, "\n\n")
  
  # Check files exist
  if (!file.exists(plot1_file)) stop("Plot 1 not found: ", plot1_file)
  if (!file.exists(plot2_file)) stop("Plot 2 not found: ", plot2_file)
  
  # Find ImageMagick convert command
  convert_cmd <- NULL
  
  # Try common locations
  possible_paths <- c(
    Sys.which("convert"),  # In PATH
    "/nethome/kcni/xzhou/.anaconda3/envs/test/bin/convert",  # test env
    "/nethome/kcni/xzhou/.anaconda3/bin/convert",  # base env
    "/usr/bin/convert"  # system
  )
  
  for (path in possible_paths) {
    if (path != "" && file.exists(path)) {
      convert_cmd <- path
      break
    }
  }
  
  if (is.null(convert_cmd)) {
    stop("ImageMagick 'convert' command not found. Please install ImageMagick:\n",
         "  conda install -c conda-forge imagemagick")
  }
  
  cat("Using ImageMagick at:", convert_cmd, "\n")
  
  cat("Combining images...\n")
  
  # Create temporary files with labels added
  temp1 <- tempfile(fileext = ".png")
  temp2 <- tempfile(fileext = ".png")
  
  # Add label to first image (use default font, more compatible)
  label1_cmd <- sprintf('"%s" "%s" -gravity NorthWest -pointsize 80 -fill black -weight 700 -annotate +50+50 "%s" "%s"',
                        convert_cmd, plot1_file, label1, temp1)
  cat("Adding label to plot 1...\n")
  system(label1_cmd)
  
  # Add label to second image  
  label2_cmd <- sprintf('"%s" "%s" -gravity NorthWest -pointsize 80 -fill black -weight 700 -annotate +50+50 "%s" "%s"',
                        convert_cmd, plot2_file, label2, temp2)
  cat("Adding label to plot 2...\n")
  system(label2_cmd)
  
  # Stack images vertically
  cat("Stacking images vertically...\n")
  combine_cmd <- sprintf('"%s" "%s" "%s" -append "%s"', convert_cmd, temp1, temp2, output_file)
  system(combine_cmd)
  
  # Clean up temp files
  unlink(c(temp1, temp2))
  
  cat("\n============================================================\n")
  cat("Combined plot saved:", output_file, "\n")
  if (file.exists(output_file)) {
    cat("File size:", round(file.size(output_file) / 1024^2, 2), "MB\n")
  }
  cat("============================================================\n")
  
  quit(save = "no", status = 0)
}

# ============================================================
# PLOT MODE - Generate Manhattan plot from GWAS data
# ============================================================

if (length(args) < 3) {
  stop("Usage: Rscript plot_manhattan_annotate.R <gwas_file> <output_prefix> <cell_type> [full]")
}

gwas_file <- args[1]
output_prefix <- args[2]
cell_type <- args[3]
use_full_data <- length(args) >= 4 && tolower(args[4]) == "full"

cat("============================================================\n")
cat("GWAS Manhattan Plot & Gene Annotation\n")
cat("============================================================\n")
cat("Cell type:", cell_type, "\n")
cat("Input file:", gwas_file, "\n")
cat("Output prefix:", output_prefix, "\n\n")

# Read GWAS results
cat("Reading GWAS results...\n")
gwas <- fread(gwas_file, data.table = FALSE)
cat("Total variants:", nrow(gwas), "\n")

# Clean data and prepare for plotting
# IMPORTANT: manhattanExtra requires "chr1", "chr2" format (with "chr" prefix)
gwas <- gwas %>%
  mutate(
    CHROM = paste0("chr", CHROM),  # Add "chr" prefix for manhattanExtra
    POS = as.numeric(GENPOS)
  ) %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1) %>%
  select(ID, CHROM, POS, P, BETA, everything())

cat("Variants after QC:", nrow(gwas), "\n")

# Lambda (genomic inflation)
chi2 <- qchisq(1 - gwas$P, 1)
lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)
cat("Lambda =", round(lambda, 4), "\n")

# Identify significant SNPs
sig_threshold <- 5e-8
suggestive_threshold <- 1e-5
sig_snps <- gwas %>% filter(P < sig_threshold)
suggestive_snps <- gwas %>% filter(P < suggestive_threshold & P >= sig_threshold)

cat("\nSignificance summary:\n")
cat("  Genome-wide significant (P < 5e-8):", nrow(sig_snps), "\n")
cat("  Suggestive (P < 1e-5):", nrow(suggestive_snps), "\n")

# Get top 20 SNPs for annotation
top_snps <- gwas %>%
  arrange(P) %>%
  head(20)

cat("\nTop 20 SNPs:\n")
print(top_snps %>% select(CHROM, POS, ID, P, BETA))

# SPEED OPTIMIZATION: Downsample variants for faster plotting (unless 'full' mode)
if (use_full_data) {
  cat("\nUsing FULL dataset (all variants, no downsampling)...\n")
  gwas_plot <- gwas
  cat("Plotting with ALL", nrow(gwas_plot), "variants\n")
} else {
  cat("\nDownsampling variants for faster plotting...\n")
  # Keep all significant/suggestive + random sample of the rest
  set.seed(42)
  important_snps <- gwas %>% filter(P < 1e-5)  # Keep all suggestive SNPs (p < 1e-5)
  other_snps <- gwas %>% filter(P >= 1e-5)
  
  # Sample 100K variants from non-significant (or all if fewer)
  if (nrow(other_snps) > 100000) {
    other_snps_sampled <- other_snps %>% sample_n(100000)
    gwas_plot <- bind_rows(important_snps, other_snps_sampled)
  } else {
    gwas_plot <- gwas
  }
  
  cat("Plotting with", nrow(gwas_plot), "variants (", nrow(important_snps), "important +", 
      nrow(gwas_plot) - nrow(important_snps), "sampled)\n")
}

# Generate annotated Manhattan plot using topr (automatically colors significant/suggestive loci)
cat("\nGenerating annotated Manhattan plot using topr...\n")
manhattan_annotated_file <- paste0(output_prefix, "_manhattan_annotated.png")

# Prepare data for topr (needs CHROM, POS, ID, P) - use downsampled for speed
topr_data <- gwas_plot %>%
  select(ID, CHROM, POS, P)

cat("Generating Manhattan plot with automatic coloring for", nrow(topr_data), "variants\n")

tryCatch({
  # manhattanExtra automatically highlights genome-wide sig and suggestive loci
  # Explicitly set thresholds to ensure proper coloring
  p_manhattan_annotated <- manhattanExtra(
    topr_data,
    genome_wide_thresh = 5e-8,
    suggestive_thresh = 1e-5,
    flank_size = 1e6,
    annotate = 1e-5,
    build = 37
  )
  
  ggsave(
    filename = manhattan_annotated_file,
    plot = p_manhattan_annotated,
    width = 16,
    height = 6,
    dpi = 300
  )
  cat("Annotated Manhattan plot saved:", manhattan_annotated_file, "\n")
}, error = function(e) {
  cat("Warning: Could not create annotated plot with topr:", e$message, "\n")
  cat("Error details:", toString(e), "\n")
})

# Generate basic Manhattan plot with qqman (faster, same as generate_meta_manhattan.R)
cat("Generating basic Manhattan plot with qqman...\n")
manhattan_file <- paste0(output_prefix, "_manhattan_basic.png")

# Use downsampled data for faster plotting
# Convert CHROM back to numeric for qqman (it expects numeric chromosomes)
qqman_data <- topr_data %>%
  mutate(CHR = as.numeric(gsub("chr", "", CHROM))) %>%
  rename(BP = POS, SNP = ID) %>%
  select(SNP, CHR, BP, P)

png(file = manhattan_file, width = 4800, height = 1800, res = 300)
qqman::manhattan(qqman_data,
                 suggestiveline = -log10(suggestive_threshold),
                 genomewideline = -log10(sig_threshold),
                 cex = 0.5,
                 cex.axis = 0.8,
                 main = paste("Manhattan Plot:", cell_type))
dev.off()
cat("Basic Manhattan plot saved:", manhattan_file, "\n")

# Generate QQ plot with ALL data (need full distribution for accurate QQ)
cat("Generating QQ plot...\n")
qq_file <- paste0(output_prefix, "_qq.png")
png(file = qq_file, width = 1400, height = 1400, res = 300)
qqman::qq(gwas$P,  # Use full dataset for QQ plot
   main = paste("QQ Plot:", cell_type),
   sub = paste("Lambda =", round(lambda, 3)))
dev.off()
cat("QQ plot saved:", qq_file, "\n")

# Save summary statistics
summary_file <- paste0(output_prefix, "_summary.txt")
summary_stats <- data.frame(
  CellType = cell_type,
  N_Samples = unique(gwas$N)[1],
  N_Variants = nrow(gwas),
  Lambda = round(lambda, 4),
  Min_P = format(min(gwas$P), scientific = TRUE),
  Significant_5e8 = nrow(sig_snps),
  Suggestive_1e5 = nrow(suggestive_snps)
)
write.table(summary_stats, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary statistics saved:", summary_file, "\n")

# Save top hits (p < 1e-5)
if (nrow(suggestive_snps) > 0) {
  top_hits_file <- paste0(output_prefix, "_top_snps.txt")
  top_hits <- gwas %>%
    filter(P < 1e-5) %>%
    arrange(P) %>%
    select(CHROM, POS, ID, P, BETA, SE) %>%
    head(100)
  write.table(top_hits, file = top_hits_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Top hits saved:", top_hits_file, "\n")
}

# Print summary
cat("\n============================================================\n")
cat("Summary Statistics\n")
cat("============================================================\n")
cat("Cell type:", cell_type, "\n")
cat("N samples:", unique(gwas$N)[1], "\n")
cat("Total variants:", nrow(gwas), "\n")
cat("Lambda (genomic inflation):", round(lambda, 4), "\n")
cat("Minimum P-value:", format(min(gwas$P), scientific = TRUE), "\n")
cat("Genome-wide significant (P < 5e-8):", nrow(sig_snps), "\n")
cat("Suggestive (P < 1e-5):", nrow(suggestive_snps), "\n")
cat("\nOutput files created:\n")
cat("  -", manhattan_file, "\n")
cat("  -", qq_file, "\n")
cat("  -", summary_file, "\n")
if (exists("manhattan_annotated_file") && file.exists(manhattan_annotated_file)) {
  cat("  -", manhattan_annotated_file, "\n")
}
if (exists("top_hits_file") && file.exists(top_hits_file)) {
  cat("  -", top_hits_file, "\n")
}
cat("============================================================\n")

cat("\nVisualization complete!\n")
