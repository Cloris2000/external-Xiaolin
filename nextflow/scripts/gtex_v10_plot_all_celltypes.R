#!/usr/bin/env Rscript
# GTEx v10 GWAS – batch plots + summary table for all 19 cell types
#
# For each cell type this script:
#   1. Calls generate_cohort_manhattan.R (reused pipeline script) which produces:
#        - <ct>_manhattan_annotated.png  (topr, gene-annotated)
#   2. Generates a QQ plot via qqman (with lambda GC annotation)
#   3. Computes genomic inflation factor (lambda GC)
#   4. Computes simple SNP heritability estimate:
#        h2_simple = (mean(chi2) - 1) * M / N
#      where M = # SNPs, N = median sample size
#   5. Writes a single combined summary TSV
#
# Usage (run from the nextflow project root):
#   Rscript scripts/gtex_v10_plot_all_celltypes.R

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(qqman)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
script_dir  <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/scripts"
input_dir   <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/GTEx_v10/regenie_step2"
output_dir  <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/GTEx_v10/GWAS_plots"
cohort_script <- file.path(script_dir, "generate_cohort_manhattan.R")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Discover cell types ───────────────────────────────────────────────────────
raw_p_files <- list.files(input_dir, pattern = "\\.regenie\\.raw_p$", full.names = TRUE)
if (length(raw_p_files) == 0) stop("No .regenie.raw_p files found in: ", input_dir)

cell_types <- sub(".*GTEx_v10_(.+)_step2\\.regenie\\.raw_p$", "\\1", basename(raw_p_files))
names(raw_p_files) <- cell_types
cat(sprintf("Found %d cell types: %s\n\n", length(cell_types), paste(sort(cell_types), collapse = ", ")))

# ── Helper: h2 simple estimate ────────────────────────────────────────────────
compute_h2_simple <- function(chi2_vec, n_vec) {
    chi2_vec <- chi2_vec[!is.na(chi2_vec) & is.finite(chi2_vec)]
    if (length(chi2_vec) < 100) return(NA_real_)
    M <- length(chi2_vec)
    N <- median(n_vec, na.rm = TRUE)
    round(max((mean(chi2_vec) - 1) * M / N, 0), 6)
}

# ── Main loop ─────────────────────────────────────────────────────────────────
summary_rows <- list()

for (ct in sort(cell_types)) {
    f <- raw_p_files[[ct]]
    cat(sprintf("==============================\n[%d/%d] %s\n", 
                which(sort(cell_types) == ct), length(cell_types), ct))

    dt <- tryCatch(
        fread(f, showProgress = FALSE),
        error = function(e) { warning("Failed to read ", f, ": ", e$message); NULL }
    )
    if (is.null(dt) || nrow(dt) == 0) { cat("  SKIP: empty file\n"); next }

    # Ensure P column exists
    if (!"P" %in% names(dt) && "LOG10P" %in% names(dt)) {
        dt[, P := 10^(-LOG10P)]
    }

    pvals  <- dt$P[!is.na(dt$P) & dt$P > 0 & dt$P <= 1]
    chi2   <- qchisq(1 - pvals, df = 1)
    lambda <- round(median(chi2, na.rm = TRUE) / qchisq(0.5, 1), 4)
    h2     <- compute_h2_simple(chi2, dt$N)
    n_snps <- nrow(dt)
    n_sig  <- sum(pvals < 5e-8, na.rm = TRUE)
    n_sug  <- sum(pvals < 1e-5, na.rm = TRUE)
    n_med  <- as.integer(median(dt$N, na.rm = TRUE))

    cat(sprintf("  SNPs: %d | N_median: %d | lambda GC: %.4f | h2_simple: %.6f | GW-sig: %d | Suggestive: %d\n",
                n_snps, n_med, lambda, ifelse(is.na(h2), 0, h2), n_sig, n_sug))

    # ── 1. Manhattan plot (via existing pipeline script) ──────────────────────
    manhattan_prefix <- paste0("GTEx_v10_", ct)
    cat("  Generating Manhattan plot...\n")
    manhattan_cmd <- sprintf(
        'Rscript "%s" --input "%s" --output_prefix "%s" --cell_type "%s" --cohort "GTEx_v10" --output_dir "%s"',
        cohort_script, f, manhattan_prefix, ct, output_dir
    )
    ret <- system(manhattan_cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
    if (ret != 0) cat(sprintf("  WARNING: Manhattan script exited with code %d\n", ret))

    # ── 2. QQ plot (qqman) ────────────────────────────────────────────────────
    qq_file <- file.path(output_dir, paste0("GTEx_v10_", ct, "_qq.png"))
    cat("  Generating QQ plot...\n")
    tryCatch({
        png(file = qq_file, width = 1400, height = 1400, res = 150)
        qqman::qq(pvals,
                  main = paste0("QQ Plot – GTEx v10 / ", ct),
                  sub  = paste0("lambda GC = ", lambda))
        dev.off()
        cat(sprintf("  Saved QQ: %s\n", basename(qq_file)))
    }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        cat(sprintf("  WARNING: QQ plot failed: %s\n", e$message))
    })

    summary_rows[[ct]] <- data.frame(
        cell_type          = ct,
        n_snps             = n_snps,
        n_samples_median   = n_med,
        lambda_gc          = lambda,
        h2_simple          = ifelse(is.na(h2), 0, h2),
        n_genome_wide_sig  = n_sig,
        n_suggestive       = n_sug,
        stringsAsFactors   = FALSE
    )
}

# ── Combined summary table ────────────────────────────────────────────────────
summary_df  <- bind_rows(summary_rows) %>% arrange(cell_type)
summary_out <- file.path(output_dir, "gtex_v10_gwas_summary.tsv")
fwrite(summary_df, summary_out, sep = "\t")

cat("\n========================================\n")
cat("All cell types done.\n")
cat("Summary table:", summary_out, "\n")
cat("========================================\n")
print(as.data.frame(summary_df), row.names = FALSE)
cat("========================================\n")
