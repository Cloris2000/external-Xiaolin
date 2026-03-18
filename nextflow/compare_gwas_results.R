#!/usr/bin/env Rscript
# Compare GWAS Results: Previous Analysis vs Nextflow Pipeline
# Generates Manhattan plots and correlation statistics
# Uses topr package for manhattanExtra plots (consistent with previous analysis)

suppressPackageStartupMessages({
  library(tidyverse)
  library(topr)
  library(gridExtra)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--old_dir"), type="character", 
              default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2",
              help="Directory with previous analysis results (default: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2)"),
  make_option(c("--new_dir"), type="character",
              default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step2",
              help="Directory with Nextflow pipeline results"),
  make_option(c("--output_dir"), type="character",
              default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/comparison",
              help="Output directory for comparison plots"),
  make_option(c("--cell_types"), type="character",
              default="Astrocyte,Microglia,Oligodendrocyte,OPC,Endothelial,Pericyte,VLMC,IT,L4.IT,L5.ET,L5.6.IT.Car3,L5.6.NP,L6.CT,L6b,LAMP5,PAX6,PVALB,SST,VIP",
              help="Comma-separated list of cell types to compare"),
  make_option(c("--maf_threshold"), type="numeric",
              default=0.01,
              help="Minor allele frequency threshold (default: 0.01)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Parse cell types
cell_types <- strsplit(opt$cell_types, ",")[[1]]

cat("=================================================\n")
cat("GWAS Results Comparison\n")
cat("=================================================\n")
cat("Previous analysis:", opt$old_dir, "\n")
cat("Nextflow pipeline:", opt$new_dir, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Cell types:", length(cell_types), "\n")
cat("MAF threshold:", opt$maf_threshold, "\n\n")

# Function to read regenie file and prepare for manhattanExtra
read_and_prepare_regenie <- function(file_path, maf_threshold = 0.01) {
  cat("  Reading:", basename(file_path), "...")
  
  tryCatch({
    # Read regenie file (whitespace separated)
    data <- read.csv(file_path, sep = "")
    
    cat(" (", nrow(data), "variants)")
    
    # Filter for MAF threshold (consistent with reference script)
    if("A1FREQ" %in% colnames(data)) {
      data <- data %>% filter(A1FREQ > maf_threshold & A1FREQ <= (1 - maf_threshold))
      cat(" -> MAF filtered (", nrow(data), "variants)")
    }
    
    # Convert LOG10P to p-value (consistent with reference script)
    if("LOG10P" %in% colnames(data)) {
      data$P <- 10^(-1 * data$LOG10P)
    } else if(!"P" %in% colnames(data)) {
      stop("No LOG10P or P column found")
    }
    
    # Prepare data for manhattanExtra (needs: ID, CHROM, POS, P)
    df_plot <- data %>%
      dplyr::select(
        ID = ID,
        CHROM = CHROM,
        POS = GENPOS,
        P = P
      ) %>%
      mutate(
        CHROM = as.numeric(CHROM),
        POS = as.numeric(POS),
        P = as.numeric(P)
      ) %>%
      filter(!is.na(CHROM) & !is.na(POS) & !is.na(P)) %>%
      filter(CHROM >= 1 & CHROM <= 22)  # Autosomes only
    
    cat(" -> Final (", nrow(df_plot), "variants)\n")
    return(df_plot)
    
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
    return(NULL)
  })
}

# Function to create Manhattan plot using topr (consistent with reference script)
create_manhattan_plot <- function(data, title) {
  tryCatch({
    plot <- manhattanExtra(
      df = data,
      genome_wide_thresh = 5e-08,
      suggestive_thresh = 1e-05,
      annotate = 1e-05,
      build = 37,
      label_size = 4,
      title = title
    )
    return(plot)
  }, error = function(e) {
    cat("  ERROR creating Manhattan plot:", e$message, "\n")
    return(NULL)
  })
}

# Main comparison loop
comparison_stats <- data.frame()

for(cell_type in cell_types) {
  cat("\n-------------------------------------------------\n")
  cat("Processing:", cell_type, "\n")
  cat("-------------------------------------------------\n")
  
  # Find old file (previous analysis)
  # Previous naming: ROSMAP_WGS_step2_update_<CellType>.regenie
  # Escape special characters in cell type for regex
  cell_type_escaped <- gsub("\\.", "\\\\.", cell_type)
  old_pattern1 <- paste0(".*_", cell_type_escaped, "\\.regenie$")
  old_pattern2 <- paste0(".*", cell_type_escaped, ".*\\.regenie$")
  
  all_files <- list.files(opt$old_dir, pattern = "\\.regenie$", full.names = TRUE, recursive = FALSE)
  old_files1 <- all_files[grepl(old_pattern1, basename(all_files), ignore.case = TRUE)]
  old_files2 <- all_files[grepl(old_pattern2, basename(all_files), ignore.case = TRUE)]
  old_files <- unique(c(old_files1, old_files2))
  old_file <- old_files[!grepl("\\.raw_p$", old_files)][1]
  
  # Debug: show what we're looking for
  if(is.na(old_file)) {
    cat("  DEBUG: Looking for pattern matching '", cell_type, "' in", opt$old_dir, "\n")
    cat("  DEBUG: Found files:", paste(basename(all_files), collapse=", "), "\n")
  }
  
  # Find new file (Nextflow pipeline)
  # New naming: ROSMAP_<CellType>_step2_<CellType>.regenie
  new_pattern <- paste0("*", cell_type, "*.regenie")
  new_files <- list.files(opt$new_dir, pattern = gsub("\\.", "\\\\.", new_pattern), full.names = TRUE, recursive = FALSE)
  new_file <- new_files[!grepl("\\.raw_p$", new_files)][1]
  
  if(is.na(old_file) || is.na(new_file)) {
    cat("  WARNING: Missing files for", cell_type, "\n")
    cat("    Previous:", if(is.na(old_file)) "NOT FOUND" else basename(old_file), "\n")
    cat("    Nextflow:", if(is.na(new_file)) "NOT FOUND" else basename(new_file), "\n")
    next
  }
  
  cat("  Files found:\n")
  cat("    Previous:", basename(old_file), "\n")
  cat("    Nextflow:", basename(new_file), "\n\n")
  
  # Read and prepare data
  old_data <- read_and_prepare_regenie(old_file, opt$maf_threshold)
  new_data <- read_and_prepare_regenie(new_file, opt$maf_threshold)
  
  if(is.null(old_data) || is.null(new_data)) {
    cat("  ERROR: Could not read data files\n")
    next
  }
  
  # Create Manhattan plots using topr (consistent with reference script)
  cat("  Creating Manhattan plots...\n")
  p1 <- create_manhattan_plot(old_data, paste0(cell_type, " - Previous Analysis"))
  p2 <- create_manhattan_plot(new_data, paste0(cell_type, " - Nextflow Pipeline"))
  
  if(is.null(p1) || is.null(p2)) {
    cat("  ERROR: Could not create Manhattan plots\n")
    next
  }
  
  # Save combined plot
  pdf_file <- file.path(opt$output_dir, paste0(cell_type, "_manhattan_comparison.pdf"))
  pdf(pdf_file, width = 10, height = 10)
  grid.arrange(p1, p2, ncol = 1)
  dev.off()
  
  cat("  Saved Manhattan plots:", basename(pdf_file), "\n")
  
  # Calculate correlation statistics
  cat("  Calculating correlation statistics...\n")
  
  # Merge by variant ID
  merged <- merge(old_data, new_data, by = "ID", suffixes = c("_old", "_new"))
  
  if(nrow(merged) > 0) {
    # Calculate -log10(P) for both
    merged$neglog10p_old <- -log10(merged$P_old)
    merged$neglog10p_new <- -log10(merged$P_new)
    
    # Remove infinite values
    merged <- merged %>% filter(is.finite(neglog10p_old) & is.finite(neglog10p_new))
    
    if(nrow(merged) > 100) {
      cor_pearson <- cor(merged$neglog10p_old, merged$neglog10p_new, method = "pearson")
      cor_spearman <- cor(merged$neglog10p_old, merged$neglog10p_new, method = "spearman")
      
      cat("    Overlapping variants:", nrow(merged), "\n")
      cat("    Pearson correlation:", round(cor_pearson, 4), "\n")
      cat("    Spearman correlation:", round(cor_spearman, 4), "\n")
      
      # Create scatter plot
      scatter_plot <- ggplot(merged, aes(x = neglog10p_old, y = neglog10p_new)) +
        geom_point(alpha = 0.3, size = 0.5) +
        geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
        labs(title = paste0(cell_type, " - P-value Correlation"),
             x = "-log10(P) Previous Analysis",
             y = "-log10(P) Nextflow Pipeline",
             subtitle = paste0("Pearson r = ", round(cor_pearson, 3), 
                              ", Spearman ρ = ", round(cor_spearman, 3), 
                              " (n = ", format(nrow(merged), big.mark = ","), ")")) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))
      
      corr_file <- file.path(opt$output_dir, paste0(cell_type, "_correlation.pdf"))
      ggsave(corr_file, scatter_plot, width = 8, height = 8)
      cat("  Saved correlation plot:", basename(corr_file), "\n")
      
      # Save stats
      comparison_stats <- rbind(comparison_stats, data.frame(
        cell_type = cell_type,
        n_variants = nrow(merged),
        pearson_r = cor_pearson,
        spearman_rho = cor_spearman,
        old_file = basename(old_file),
        new_file = basename(new_file)
      ))
    } else {
      cat("    WARNING: Too few overlapping variants (", nrow(merged), ") for correlation\n")
    }
  } else {
    cat("    WARNING: No overlapping variants found\n")
  }
}

# Save comparison statistics
if(nrow(comparison_stats) > 0) {
  stats_file <- file.path(opt$output_dir, "comparison_statistics.csv")
  write.csv(comparison_stats, stats_file, row.names = FALSE)
  
  cat("\n=================================================\n")
  cat("Summary Statistics\n")
  cat("=================================================\n")
  print(comparison_stats)
  cat("\nSaved comparison statistics to:", basename(stats_file), "\n")
}

cat("\n=================================================\n")
cat("Comparison complete!\n")
cat("=================================================\n")
cat("Output directory:", opt$output_dir, "\n")
cat("- Manhattan comparison plots: *_manhattan_comparison.pdf\n")
cat("- Correlation plots: *_correlation.pdf\n")
cat("- Summary statistics: comparison_statistics.csv\n\n")

