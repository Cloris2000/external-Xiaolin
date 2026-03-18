#!/usr/bin/env Rscript

# Comprehensive Diagnostic Script for GWAS Result Discrepancies
# Purpose: Identify why Pearson r = 0.78 instead of ~0.99

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

# Setup
old_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2"
new_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step2"
output_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/diagnostics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cell_type <- "Astrocyte"

cat("\n=================================================\n")
cat("GWAS DISCREPANCY DIAGNOSTIC REPORT\n")
cat("=================================================\n")
cat("Cell Type:", cell_type, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =================================================
# STEP 1: Load and Compare Basic Statistics
# =================================================
cat("\n[STEP 1] Loading GWAS Results...\n")

old_file <- file.path(old_dir, paste0("ROSMAP_WGS_step2_update_", cell_type, ".regenie"))
new_file <- file.path(new_dir, paste0("ROSMAP_", cell_type, "_step2_", cell_type, ".regenie"))

old_data <- fread(old_file, data.table = FALSE)
new_data <- fread(new_file, data.table = FALSE)

cat("  Previous analysis: ", nrow(old_data), "variants\n")
cat("  Nextflow pipeline: ", nrow(new_data), "variants\n")

# Check column names
cat("\n  Column names comparison:\n")
cat("    Previous:", paste(names(old_data), collapse=", "), "\n")
cat("    Nextflow:", paste(names(new_data), collapse=", "), "\n")

# =================================================
# STEP 2: Sample Size Check
# =================================================
cat("\n[STEP 2] Sample Size Comparison...\n")

# Extract N from results
old_n <- old_data$N[1]
new_n <- new_data$N[1]

cat("  Previous analysis N:", old_n, "\n")
cat("  Nextflow pipeline N:", new_n, "\n")
cat("  Difference:", new_n - old_n, "\n")

if (old_n != new_n) {
  cat("  ⚠️  WARNING: Sample sizes differ! This will affect p-values.\n")
}

# =================================================
# STEP 3: Variant Overlap Analysis
# =================================================
cat("\n[STEP 3] Variant Overlap Analysis...\n")

# Create variant IDs
old_data$variant_id <- with(old_data, paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep=":"))
new_data$variant_id <- with(new_data, paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep=":"))

old_variants <- unique(old_data$variant_id)
new_variants <- unique(new_data$variant_id)

common_variants <- intersect(old_variants, new_variants)
old_only <- setdiff(old_variants, new_variants)
new_only <- setdiff(new_variants, old_variants)

cat("  Previous-only variants:", length(old_only), "\n")
cat("  Nextflow-only variants:", length(new_only), "\n")
cat("  Common variants:", length(common_variants), "\n")
cat("  Overlap percentage:", round(100 * length(common_variants) / length(old_variants), 2), "%\n")

if (length(old_only) > 0 || length(new_only) > 0) {
  cat("  ⚠️  WARNING: Variant sets differ!\n")
  
  # Save non-overlapping variants
  if (length(old_only) > 0) {
    write.table(data.frame(variant_id = old_only[1:min(100, length(old_only))]),
                file.path(output_dir, "variants_previous_only.txt"),
                row.names = FALSE, quote = FALSE)
  }
  if (length(new_only) > 0) {
    write.table(data.frame(variant_id = new_only[1:min(100, length(new_only))]),
                file.path(output_dir, "variants_nextflow_only.txt"),
                row.names = FALSE, quote = FALSE)
  }
}

# =================================================
# STEP 4: MAF Distribution Comparison
# =================================================
cat("\n[STEP 4] MAF Distribution Analysis...\n")

# Filter to common variants for fair comparison
old_common <- old_data %>% filter(variant_id %in% common_variants)
new_common <- new_data %>% filter(variant_id %in% common_variants)

cat("  Previous MAF range:", range(old_common$A1FREQ), "\n")
cat("  Nextflow MAF range:", range(new_common$A1FREQ), "\n")

# MAF comparison plot
pdf(file.path(output_dir, "maf_comparison.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))
hist(old_common$A1FREQ, breaks = 50, main = "Previous Analysis MAF", 
     xlab = "MAF", col = "lightblue")
hist(new_common$A1FREQ, breaks = 50, main = "Nextflow Pipeline MAF", 
     xlab = "MAF", col = "lightcoral")
dev.off()

cat("  ✓ Saved: maf_comparison.pdf\n")

# =================================================
# STEP 5: Merge and Compare P-values
# =================================================
cat("\n[STEP 5] P-value Correlation Analysis...\n")

# Merge datasets
merged <- merge(old_common, new_common, by = "variant_id", suffixes = c("_old", "_new"))
cat("  Merged variants:", nrow(merged), "\n")

# Calculate p-values
merged$p_old <- 10^(-merged$LOG10P_old)
merged$p_new <- 10^(-merged$LOG10P_new)

# Overall correlation
pearson_r <- cor(merged$LOG10P_old, merged$LOG10P_new, method = "pearson")
spearman_r <- cor(merged$LOG10P_old, merged$LOG10P_new, method = "spearman")

cat("  Pearson correlation (LOG10P):", round(pearson_r, 4), "\n")
cat("  Spearman correlation (LOG10P):", round(spearman_r, 4), "\n")

# =================================================
# STEP 6: Effect Size (Beta) Comparison
# =================================================
cat("\n[STEP 6] Effect Size (Beta) Comparison...\n")

beta_cor <- cor(merged$BETA_old, merged$BETA_new, method = "pearson")
cat("  Beta correlation:", round(beta_cor, 4), "\n")

# Check for sign flips
sign_flips <- sum(sign(merged$BETA_old) != sign(merged$BETA_new))
cat("  Sign flips:", sign_flips, "(", round(100*sign_flips/nrow(merged), 2), "%)\n")

if (sign_flips > nrow(merged) * 0.01) {
  cat("  ⚠️  WARNING: >1% sign flips suggest allele coding differences!\n")
}

# Beta scatter plot
pdf(file.path(output_dir, "beta_comparison.pdf"), width = 8, height = 8)
plot(merged$BETA_old, merged$BETA_new, 
     pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3),
     xlab = "Previous Beta", ylab = "Nextflow Beta",
     main = sprintf("Effect Size Comparison (r = %.3f)", beta_cor))
abline(0, 1, col = "red", lwd = 2)
abline(h = 0, v = 0, col = "gray", lty = 2)
dev.off()

cat("  ✓ Saved: beta_comparison.pdf\n")

# =================================================
# STEP 7: Standard Error Comparison
# =================================================
cat("\n[STEP 7] Standard Error Comparison...\n")

se_cor <- cor(merged$SE_old, merged$SE_new, method = "pearson")
cat("  SE correlation:", round(se_cor, 4), "\n")

se_ratio <- mean(merged$SE_new / merged$SE_old)
cat("  Mean SE ratio (new/old):", round(se_ratio, 4), "\n")

if (abs(se_ratio - 1) > 0.1) {
  cat("  ⚠️  WARNING: SE ratio deviates >10% from 1.0!\n")
  cat("     This suggests different effective sample sizes or model specifications.\n")
}

# =================================================
# STEP 8: Stratified Analysis by MAF
# =================================================
cat("\n[STEP 8] Stratified Correlation by MAF...\n")

merged$maf_bin <- cut(merged$A1FREQ_old, 
                      breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.5),
                      labels = c("0-1%", "1-5%", "5-10%", "10-20%", "20-50%"))

maf_strat <- merged %>%
  group_by(maf_bin) %>%
  summarise(
    n_variants = n(),
    pearson_r = cor(LOG10P_old, LOG10P_new, method = "pearson"),
    spearman_r = cor(LOG10P_old, LOG10P_new, method = "spearman"),
    mean_se_ratio = mean(SE_new / SE_old)
  )

print(maf_strat)

write.csv(maf_strat, file.path(output_dir, "maf_stratified_correlation.csv"), 
          row.names = FALSE)

# =================================================
# STEP 9: Identify Outlier Variants
# =================================================
cat("\n[STEP 9] Identifying Outlier Variants...\n")

# Calculate residuals from y=x line
merged$log10p_diff <- merged$LOG10P_new - merged$LOG10P_old
merged$abs_diff <- abs(merged$log10p_diff)

# Top 100 discrepant variants
outliers <- merged %>%
  arrange(desc(abs_diff)) %>%
  head(100) %>%
  select(variant_id, CHROM_old, GENPOS_old, A1FREQ_old, 
         LOG10P_old, LOG10P_new, log10p_diff, BETA_old, BETA_new,
         SE_old, SE_new, N_old, N_new)

write.csv(outliers, file.path(output_dir, "top_100_outliers.csv"), row.names = FALSE)
cat("  ✓ Saved: top_100_outliers.csv\n")

# Check if outliers are enriched in specific regions
outlier_chroms <- table(outliers$CHROM_old)
cat("\n  Outlier distribution by chromosome:\n")
print(outlier_chroms)

# =================================================
# STEP 10: QQ Plot Comparison
# =================================================
cat("\n[STEP 10] Generating QQ Plots...\n")

# Prepare data
merged_sorted_old <- merged %>% arrange(p_old)
merged_sorted_new <- merged %>% arrange(p_new)

merged_sorted_old$expected_p <- (1:nrow(merged_sorted_old)) / (nrow(merged_sorted_old) + 1)
merged_sorted_new$expected_p <- (1:nrow(merged_sorted_new)) / (nrow(merged_sorted_new) + 1)

pdf(file.path(output_dir, "qq_plot_comparison.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# Previous analysis QQ plot
plot(-log10(merged_sorted_old$expected_p), -log10(merged_sorted_old$p_old),
     pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3),
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
     main = "Previous Analysis QQ Plot")
abline(0, 1, col = "red", lwd = 2)

# Nextflow QQ plot
plot(-log10(merged_sorted_new$expected_p), -log10(merged_sorted_new$p_new),
     pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.3),
     xlab = "Expected -log10(P)", ylab = "Observed -log10(P)",
     main = "Nextflow Pipeline QQ Plot")
abline(0, 1, col = "red", lwd = 2)

dev.off()
cat("  ✓ Saved: qq_plot_comparison.pdf\n")

# =================================================
# STEP 11: Bland-Altman Plot
# =================================================
cat("\n[STEP 11] Bland-Altman Analysis...\n")

merged$log10p_mean <- (merged$LOG10P_old + merged$LOG10P_new) / 2
merged$log10p_diff <- merged$LOG10P_new - merged$LOG10P_old

mean_diff <- mean(merged$log10p_diff)
sd_diff <- sd(merged$log10p_diff)

pdf(file.path(output_dir, "bland_altman_plot.pdf"), width = 10, height = 8)
plot(merged$log10p_mean, merged$log10p_diff,
     pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.2),
     xlab = "Mean -log10(P) [(Old + New) / 2]",
     ylab = "Difference -log10(P) [New - Old]",
     main = "Bland-Altman Plot")
abline(h = mean_diff, col = "blue", lwd = 2)
abline(h = mean_diff + 1.96 * sd_diff, col = "red", lty = 2, lwd = 2)
abline(h = mean_diff - 1.96 * sd_diff, col = "red", lty = 2, lwd = 2)
legend("topright", 
       legend = c("Mean difference", "±1.96 SD"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)
dev.off()

cat("  Mean difference (new - old):", round(mean_diff, 4), "\n")
cat("  SD of differences:", round(sd_diff, 4), "\n")
cat("  ✓ Saved: bland_altman_plot.pdf\n")

# =================================================
# STEP 12: INFO Score / Imputation Quality Check
# =================================================
cat("\n[STEP 12] Imputation Quality Check...\n")

if ("INFO" %in% names(merged)) {
  info_diff <- merged$INFO_new - merged$INFO_old
  cat("  Mean INFO difference:", round(mean(info_diff), 4), "\n")
  
  if (abs(mean(info_diff)) > 0.01) {
    cat("  ⚠️  WARNING: INFO scores differ between analyses!\n")
  }
} else {
  cat("  INFO column not available in results.\n")
}

# =================================================
# STEP 13: Check for Systematic Bias
# =================================================
cat("\n[STEP 13] Testing for Systematic Bias...\n")

# Linear regression: new ~ old
lm_fit <- lm(LOG10P_new ~ LOG10P_old, data = merged)
intercept <- coef(lm_fit)[1]
slope <- coef(lm_fit)[2]

cat("  Linear model: LOG10P_new = ", round(intercept, 4), " + ", 
    round(slope, 4), " × LOG10P_old\n", sep = "")

if (abs(intercept) > 0.1) {
  cat("  ⚠️  WARNING: Non-zero intercept suggests systematic bias!\n")
}

if (abs(slope - 1) > 0.05) {
  cat("  ⚠️  WARNING: Slope ≠ 1 suggests scaling difference!\n")
}

# =================================================
# STEP 14: Summary Report
# =================================================
cat("\n=================================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("=================================================\n\n")

cat("KEY FINDINGS:\n")
cat("1. Sample sizes: ", 
    ifelse(old_n == new_n, "✓ Identical", "✗ Different"), "\n", sep = "")
cat("2. Variant overlap: ", 
    round(100 * length(common_variants) / length(old_variants), 1), "%\n", sep = "")
cat("3. P-value correlation: ", round(pearson_r, 3), " (Pearson)\n", sep = "")
cat("4. Effect size correlation: ", round(beta_cor, 3), "\n", sep = "")
cat("5. SE ratio (new/old): ", round(se_ratio, 3), "\n", sep = "")
cat("6. Sign flips: ", round(100*sign_flips/nrow(merged), 2), "%\n", sep = "")
cat("7. Systematic bias: ", 
    ifelse(abs(intercept) < 0.1 & abs(slope - 1) < 0.05, "✓ None detected", "✗ Present"), 
    "\n", sep = "")

cat("\nPOSSIBLE EXPLANATIONS FOR r = 0.78:\n")

explanations <- c()
if (old_n != new_n) {
  explanations <- c(explanations, "  • Different sample sizes affecting p-values")
}
if (length(common_variants) < 0.99 * length(old_variants)) {
  explanations <- c(explanations, "  • Different variant sets (QC filters)")
}
if (sign_flips > nrow(merged) * 0.01) {
  explanations <- c(explanations, "  • Allele coding differences (sign flips)")
}
if (abs(se_ratio - 1) > 0.1) {
  explanations <- c(explanations, "  • Different model specifications or covariates")
}
if (abs(slope - 1) > 0.05) {
  explanations <- c(explanations, "  • Systematic scaling difference")
}

if (length(explanations) > 0) {
  cat(paste(explanations, collapse = "\n"), "\n")
} else {
  cat("  • No obvious major differences detected\n")
  cat("  • Discrepancy may be due to:\n")
  cat("    - Minor numerical differences in algorithms\n")
  cat("    - Different covariate values\n")
  cat("    - Different phenotype normalizations\n")
}

cat("\nOUTPUT FILES GENERATED:\n")
cat("  Directory:", output_dir, "\n")
cat("  - maf_comparison.pdf\n")
cat("  - beta_comparison.pdf\n")
cat("  - qq_plot_comparison.pdf\n")
cat("  - bland_altman_plot.pdf\n")
cat("  - maf_stratified_correlation.csv\n")
cat("  - top_100_outliers.csv\n")
if (length(old_only) > 0) cat("  - variants_previous_only.txt\n")
if (length(new_only) > 0) cat("  - variants_nextflow_only.txt\n")

cat("\n=================================================\n")
cat("Diagnostic analysis complete!\n")
cat("=================================================\n\n")

