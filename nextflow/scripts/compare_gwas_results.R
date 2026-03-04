#!/usr/bin/env Rscript
# Compare GWAS results between Nextflow pipeline and individual scripts
# Calculate total significant/suggestive variants and overlaps

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("============================================================\n")
cat("GWAS Results Comparison: Nextflow Pipeline vs Individual Scripts\n")
cat("============================================================\n\n")

# File paths
nf_base <- "results/ROSMAP/regenie_step2"
wgs_base <- "manhattan_plots"

cell_types <- c("VIP", "L5.ET", "SST")

# Initialize results table
results <- data.frame()

for (ct in cell_types) {
  cat("\n", rep("=", 60), "\n", sep="")
  cat("Processing:", ct, "\n")
  cat(rep("=", 60), "\n", sep="")
  
  # Read Nextflow pipeline results
  nf_file <- file.path(nf_base, paste0("ROSMAP_", ct, "_step2.regenie.raw_p"))
  cat("  Nextflow file:", nf_file, "\n")
  nf <- fread(nf_file, data.table = FALSE)
  
  # Read WGS individual scripts results
  wgs_file <- file.path(wgs_base, paste0("ROSMAP_WGS_", ct, ".raw_p"))
  cat("  WGS file:     ", wgs_file, "\n\n")
  wgs <- fread(wgs_file, data.table = FALSE)
  
  # Calculate significance counts
  nf_sig <- sum(nf$P < 5e-8, na.rm = TRUE)
  nf_sug <- sum(nf$P >= 5e-8 & nf$P < 1e-5, na.rm = TRUE)
  nf_total_hits <- sum(nf$P < 1e-5, na.rm = TRUE)
  
  wgs_sig <- sum(wgs$P < 5e-8, na.rm = TRUE)
  wgs_sug <- sum(wgs$P >= 5e-8 & wgs$P < 1e-5, na.rm = TRUE)
  wgs_total_hits <- sum(wgs$P < 1e-5, na.rm = TRUE)
  
  cat("  Nextflow Pipeline:\n")
  cat("    Genome-wide significant (P < 5e-8):", nf_sig, "\n")
  cat("    Suggestive (5e-8 ≤ P < 1e-5):      ", nf_sug, "\n")
  cat("    Total hits (P < 1e-5):             ", nf_total_hits, "\n\n")
  
  cat("  Individual Scripts (WGS):\n")
  cat("    Genome-wide significant (P < 5e-8):", wgs_sig, "\n")
  cat("    Suggestive (5e-8 ≤ P < 1e-5):      ", wgs_sug, "\n")
  cat("    Total hits (P < 1e-5):             ", wgs_total_hits, "\n\n")
  
  # Find overlapping SNPs (by position)
  # Create position keys
  nf$pos_key <- paste(nf$CHROM, nf$GENPOS, sep=":")
  wgs$pos_key <- paste(wgs$CHROM, wgs$GENPOS, sep=":")
  
  # Get suggestive/significant SNPs from each
  nf_hits <- nf %>% filter(P < 1e-5) %>% select(pos_key, P, BETA)
  wgs_hits <- wgs %>% filter(P < 1e-5) %>% select(pos_key, P, BETA)
  
  # Find overlaps
  overlap <- inner_join(nf_hits, wgs_hits, by = "pos_key", suffix = c("_NF", "_WGS"))
  n_overlap <- nrow(overlap)
  
  # Calculate concordance
  if (n_overlap > 0) {
    # Check if effect directions match
    concordant <- sum(sign(overlap$BETA_NF) == sign(overlap$BETA_WGS))
    concordance_pct <- round(100 * concordant / n_overlap, 1)
  } else {
    concordant <- 0
    concordance_pct <- NA
  }
  
  cat("  Overlap Analysis:\n")
  cat("    Shared suggestive/significant SNPs:", n_overlap, "\n")
  cat("    Concordant effect direction:       ", concordant, "(", concordance_pct, "%)\n\n")
  
  # Calculate replication rate
  rep_rate_nf <- round(100 * n_overlap / max(nf_total_hits, 1), 1)
  rep_rate_wgs <- round(100 * n_overlap / max(wgs_total_hits, 1), 1)
  
  cat("  Replication Rate:\n")
  cat("    NF hits replicated in WGS:  ", rep_rate_nf, "%\n")
  cat("    WGS hits replicated in NF:  ", rep_rate_wgs, "%\n")
  
  # Store results
  results <- rbind(results, data.frame(
    CellType = ct,
    NF_Significant = nf_sig,
    NF_Suggestive = nf_sug,
    NF_Total = nf_total_hits,
    WGS_Significant = wgs_sig,
    WGS_Suggestive = wgs_sug,
    WGS_Total = wgs_total_hits,
    Shared_SNPs = n_overlap,
    Concordant = concordant,
    Concordance_Pct = concordance_pct,
    NF_Replication_Pct = rep_rate_nf,
    WGS_Replication_Pct = rep_rate_wgs
  ))
}

cat("\n", rep("=", 60), "\n", sep="")
cat("SUMMARY TABLE\n")
cat(rep("=", 60), "\n\n", sep="")

print(results, row.names = FALSE)

# Save summary table
output_file <- "manhattan_plots/GWAS_comparison_summary.txt"
write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n\nSummary table saved to:", output_file, "\n")

cat("\n", rep("=", 60), "\n", sep="")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n")
