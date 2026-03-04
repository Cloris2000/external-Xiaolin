#!/usr/bin/env Rscript
# GWAS Results Visualization and Report Generation
# Generates Manhattan plots, QQ plots, and summary statistics for GWAS results

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(qqman)
  library(gridExtra)
  library(knitr)
  library(rmarkdown)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--results_dir"), type="character", default=NULL,
              help="Directory containing Regenie Step 2 results"),
  make_option(c("--meta_dir"), type="character", default=NULL,
              help="Directory containing meta-analysis results"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory for visualizations"),
  make_option(c("--cohort"), type="character", default="ROSMAP",
              help="Cohort name")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Define cell types
cell_types <- c(
  "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.ET", "L5.6.IT.Car3",
  "L5.6.NP", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC",
  "Oligodendrocyte", "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC"
)

# Function to read Regenie results
read_regenie_results <- function(file_path) {
  # Regenie output format: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
  # Skip header and read data
  data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
  
  # Convert LOG10P to P-value
  data$P <- 10^(-data$LOG10P)
  
  # Ensure required columns exist
  required_cols <- c("CHROM", "GENPOS", "ID", "P", "LOG10P")
  if (!all(required_cols %in% colnames(data))) {
    stop(paste("Missing required columns in", file_path))
  }
  
  return(data)
}

# Function to generate Manhattan plot
generate_manhattan <- function(data, cell_type, output_path) {
  # Prepare data for qqman
  gwas_data <- data.frame(
    CHR = data$CHROM,
    BP = data$GENPOS,
    SNP = data$ID,
    P = data$P
  )
  
  # Remove any NA values
  gwas_data <- gwas_data[complete.cases(gwas_data), ]
  
  if (nrow(gwas_data) == 0) {
    warning(paste("No valid data for", cell_type))
    return(NULL)
  }
  
  # Generate Manhattan plot
  png(file = output_path, width = 2400, height = 1200, res = 300)
  manhattan(gwas_data, 
            main = paste("Manhattan Plot:", cell_type, "-", opt$cohort),
            suggestiveline = -log10(5e-8),
            genomewideline = -log10(5e-8),
            cex = 0.6)
  dev.off()
  
  return(output_path)
}

# Function to generate QQ plot
generate_qq <- function(data, cell_type, output_path) {
  # Prepare data
  p_values <- data$P[!is.na(data$P) & data$P > 0 & data$P <= 1]
  
  if (length(p_values) == 0) {
    warning(paste("No valid p-values for", cell_type))
    return(NULL)
  }
  
  # Calculate lambda (genomic inflation factor)
  chi2 <- qchisq(1 - p_values, 1)
  lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)
  
  # Generate QQ plot
  png(file = output_path, width = 1200, height = 1200, res = 300)
  qq(p_values, 
     main = paste("QQ Plot:", cell_type, "-", opt$cohort),
     sub = paste("Lambda =", round(lambda, 3)))
  dev.off()
  
  return(list(path = output_path, lambda = lambda))
}

# Function to generate summary statistics
generate_summary_stats <- function(data, cell_type) {
  p_values <- data$P[!is.na(data$P) & data$P > 0 & data$P <= 1]
  
  if (length(p_values) == 0) {
    return(NULL)
  }
  
  # Calculate lambda
  chi2 <- qchisq(1 - p_values, 1)
  lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, 1)
  
  # Count significant associations
  sig_5e8 <- sum(p_values < 5e-8, na.rm = TRUE)
  sig_5e6 <- sum(p_values < 5e-6, na.rm = TRUE)
  sig_1e5 <- sum(p_values < 1e-5, na.rm = TRUE)
  
  # Min p-value
  min_p <- min(p_values, na.rm = TRUE)
  
  return(data.frame(
    CellType = cell_type,
    N_Variants = length(p_values),
    Lambda = round(lambda, 4),
    Min_P = format(min_p, scientific = TRUE),
    Significant_5e8 = sig_5e8,
    Significant_5e6 = sig_5e6,
    Significant_1e5 = sig_1e5
  ))
}

# Process each cell type
summary_stats_list <- list()
lambda_values <- numeric()

cat("Processing GWAS results...\n")
for (cell_type in cell_types) {
  cat(paste("Processing", cell_type, "...\n"))
  
  # Construct file path
  # Format: {cohort}_{cell_type}_step2_{cell_type}.regenie
  regenie_file <- file.path(opt$results_dir, 
                            paste0(opt$cohort, "_", cell_type, "_step2_", cell_type, ".regenie"))
  
  if (!file.exists(regenie_file)) {
    warning(paste("File not found:", regenie_file))
    next
  }
  
  # Read results
  gwas_data <- read_regenie_results(regenie_file)
  
  # Generate Manhattan plot
  manhattan_path <- file.path(opt$output_dir, paste0(opt$cohort, "_", cell_type, "_manhattan.png"))
  generate_manhattan(gwas_data, cell_type, manhattan_path)
  
  # Generate QQ plot
  qq_path <- file.path(opt$output_dir, paste0(opt$cohort, "_", cell_type, "_qq.png"))
  qq_result <- generate_qq(gwas_data, cell_type, qq_path)
  
  if (!is.null(qq_result)) {
    lambda_values <- c(lambda_values, qq_result$lambda)
  }
  
  # Generate summary statistics
  stats <- generate_summary_stats(gwas_data, cell_type)
  if (!is.null(stats)) {
    summary_stats_list[[cell_type]] <- stats
  }
}

# Combine summary statistics
if (length(summary_stats_list) > 0) {
  summary_df <- do.call(rbind, summary_stats_list)
  
  # Save summary statistics
  summary_file <- file.path(opt$output_dir, paste0(opt$cohort, "_summary_statistics.txt"))
  write.table(summary_df, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\nSummary Statistics:\n")
  print(summary_df)
  
  # Generate summary plot
  png(file = file.path(opt$output_dir, paste0(opt$cohort, "_lambda_distribution.png")), 
      width = 1200, height = 800, res = 300)
  hist(lambda_values, 
       main = paste("Genomic Inflation Factor (Lambda) Distribution -", opt$cohort),
       xlab = "Lambda",
       breaks = 20,
       col = "steelblue")
  abline(v = 1.0, col = "red", lty = 2, lwd = 2)
  dev.off()
  
  # Generate significant associations bar plot
  png(file = file.path(opt$output_dir, paste0(opt$cohort, "_significant_associations.png")), 
      width = 2000, height = 1200, res = 300)
  par(mar = c(10, 5, 4, 2))
  barplot(summary_df$Significant_5e8, 
          names.arg = summary_df$CellType,
          las = 2,
          main = paste("Number of Genome-wide Significant Associations (p < 5e-8) -", opt$cohort),
          ylab = "Number of Associations",
          col = "steelblue")
  dev.off()
}

cat("\nVisualization complete!\n")
cat(paste("Output directory:", opt$output_dir, "\n"))

