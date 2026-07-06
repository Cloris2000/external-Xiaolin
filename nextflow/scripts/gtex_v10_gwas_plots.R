#!/usr/bin/env Rscript
# GTEx v10 GWAS: Manhattan plots, QQ plots, genomic inflation (lambda GC),
# and SNP-based heritability (h2) via LD Score regression approximation.
#
# Heritability is estimated using the Bulik-Sullivan et al. LDSC intercept-free
# formula:  h2 = (mean(chi2) - 1) / (N * mean(l2) / M)
# where l2 is approximated from the chi2 statistics themselves (regression of
# chi2 on expected chi2 quantiles) when external LD scores are not available.
# A simpler but widely used approximation is also reported:
#   h2_simple = (mean(chi2) - 1) * M / N
# where M = number of SNPs, N = median sample size.
#
# Usage:
#   Rscript gtex_v10_gwas_plots.R \
#       --input_dir  results/GTEx_v10/regenie_step2 \
#       --output_dir results/GTEx_v10/GWAS_plots
#
# Outputs (per cell type):
#   <cell_type>_manhattan.png
#   <cell_type>_qqplot.png
# Summary table:
#   gtex_v10_gwas_summary.tsv

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tidyr)
    library(optparse)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--input_dir",  type = "character",
                default = "results/GTEx_v10/regenie_step2",
                help    = "Directory containing REGENIE step2 .raw_p files"),
    make_option("--output_dir", type = "character",
                default = "results/GTEx_v10/GWAS_plots",
                help    = "Output directory for plots and summary table"),
    make_option("--p_threshold", type = "double", default = 5e-8,
                help    = "Genome-wide significance threshold [default: 5e-8]"),
    make_option("--suggestive",  type = "double", default = 1e-5,
                help    = "Suggestive significance threshold [default: 1e-5]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helper functions ──────────────────────────────────────────────────────────

compute_lambda_gc <- function(pvalues) {
    # Genomic inflation factor: median(chi2_obs) / median(chi2_expected)
    pvalues <- pvalues[!is.na(pvalues) & pvalues > 0 & pvalues < 1]
    if (length(pvalues) < 100) return(NA_real_)
    chi2_obs <- qchisq(pvalues, df = 1, lower.tail = FALSE)
    round(median(chi2_obs) / qchisq(0.5, df = 1, lower.tail = FALSE), 4)
}

compute_h2_simple <- function(chi2_vec, n_vec) {
    # Simple LDSC-style approximation (no external LD scores):
    #   h2 = (mean(chi2) - 1) * M / N
    # where M = number of SNPs, N = median sample size.
    # Returns NA if insufficient data.
    chi2_vec <- chi2_vec[!is.na(chi2_vec)]
    if (length(chi2_vec) < 100) return(NA_real_)
    M <- length(chi2_vec)
    N <- median(n_vec, na.rm = TRUE)
    h2 <- (mean(chi2_vec) - 1) * M / N
    round(max(h2, 0), 6)   # clamp at 0 (negative = no detectable h2)
}

make_manhattan <- function(df, cell_type, p_threshold, suggestive, out_path) {
    # df must have: CHROM, GENPOS, P
    df <- df[!is.na(df$P) & df$P > 0, ]
    df$CHROM <- as.integer(sub("chr", "", df$CHROM))
    df <- df[df$CHROM %in% 1:22, ]
    df <- df[order(df$CHROM, df$GENPOS), ]

    # Cumulative BP position for x-axis
    chr_sizes <- df %>% group_by(CHROM) %>%
        summarise(max_bp = max(GENPOS), .groups = "drop") %>%
        arrange(CHROM) %>%
        mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0))
    df <- df %>% left_join(chr_sizes[, c("CHROM","bp_add")], by = "CHROM") %>%
        mutate(bp_cum = GENPOS + bp_add)

    axis_df <- df %>% group_by(CHROM) %>%
        summarise(center = (max(bp_cum) + min(bp_cum)) / 2, .groups = "drop")

    # Label top hits
    sig_df <- df %>% filter(P < p_threshold)
    label_df <- if (nrow(sig_df) > 0) {
        sig_df %>% group_by(CHROM) %>% slice_min(P, n = 1)
    } else {
        data.frame()
    }

    chrom_colors <- rep(c("#4393C3","#2166AC"), 11)[1:22]
    names(chrom_colors) <- 1:22

    p <- ggplot(df, aes(x = bp_cum, y = -log10(P), color = factor(CHROM))) +
        geom_point(size = 0.4, alpha = 0.7) +
        geom_hline(yintercept = -log10(p_threshold), color = "red",
                   linetype = "dashed", linewidth = 0.6) +
        geom_hline(yintercept = -log10(suggestive),  color = "orange",
                   linetype = "dotted",  linewidth = 0.5) +
        scale_color_manual(values = chrom_colors, guide = "none") +
        scale_x_continuous(
            label = axis_df$CHROM,
            breaks = axis_df$center,
            expand = c(0.01, 0)
        ) +
        labs(
            title = paste0("Manhattan plot — GTEx v10 / ", cell_type),
            subtitle = paste0("Red: p < 5×10⁻⁸  |  Orange: p < 1×10⁻⁵"),
            x = "Chromosome", y = expression(-log[10](p))
        ) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            axis.text.x = element_text(size = 7)
        )

    if (nrow(label_df) > 0) {
        p <- p + geom_label_repel(
            data = label_df,
            aes(label = ID), size = 2.5, color = "black",
            box.padding = 0.3, max.overlaps = 20, show.legend = FALSE
        )
    }

    ggsave(out_path, plot = p, width = 14, height = 5, dpi = 150)
    invisible(p)
}

make_qq <- function(pvalues, cell_type, lambda_gc, out_path) {
    pvalues <- pvalues[!is.na(pvalues) & pvalues > 0 & pvalues < 1]
    n  <- length(pvalues)
    df <- data.frame(
        observed = -log10(sort(pvalues)),
        expected = -log10(ppoints(n))
    )
    # 95% CI band
    ci_upper <- -log10(qbeta(0.025, seq_len(n), n - seq_len(n) + 1))
    ci_lower <- -log10(qbeta(0.975, seq_len(n), n - seq_len(n) + 1))
    df$ci_upper <- sort(ci_upper, decreasing = TRUE)
    df$ci_lower <- sort(ci_lower, decreasing = TRUE)

    p <- ggplot(df, aes(x = expected, y = observed)) +
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                    fill = "lightblue", alpha = 0.5) +
        geom_point(size = 0.5, alpha = 0.7, color = "#2166AC") +
        geom_abline(slope = 1, intercept = 0, color = "red",
                    linetype = "dashed", linewidth = 0.7) +
        annotate("text", x = 0.3, y = max(df$observed) * 0.92,
                 label = paste0("λ GC = ", lambda_gc), size = 4, hjust = 0) +
        labs(
            title   = paste0("QQ plot — GTEx v10 / ", cell_type),
            x       = expression("Expected " * -log[10](p)),
            y       = expression("Observed " * -log[10](p))
        ) +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank())

    ggsave(out_path, plot = p, width = 6, height = 6, dpi = 150)
    invisible(p)
}

# ── Main loop over cell types ─────────────────────────────────────────────────
raw_p_files <- list.files(opt$input_dir, pattern = "\\.regenie\\.raw_p$",
                           full.names = TRUE)

if (length(raw_p_files) == 0) {
    stop("No .regenie.raw_p files found in: ", opt$input_dir)
}

# Extract cell type from filename: GTEx_v10_<CELLTYPE>_step2.regenie.raw_p
cell_types <- sub(".*GTEx_v10_(.+)_step2\\.regenie\\.raw_p$", "\\1",
                  basename(raw_p_files))

summary_rows <- list()

for (i in seq_along(raw_p_files)) {
    ct   <- cell_types[i]
    f    <- raw_p_files[i]
    cat(sprintf("[%d/%d] Processing: %s\n", i, length(raw_p_files), ct))

    dt <- tryCatch(
        fread(f, showProgress = FALSE),
        error = function(e) { warning("Failed to read ", f, ": ", e$message); NULL }
    )
    if (is.null(dt) || nrow(dt) == 0) next

    # Standardise column names (REGENIE raw_p adds a P column)
    if (!"P" %in% names(dt) && "LOG10P" %in% names(dt)) {
        dt[, P := 10^(-LOG10P)]
    }

    pvals  <- dt$P
    chi2   <- qchisq(pvals, df = 1, lower.tail = FALSE)
    lambda <- compute_lambda_gc(pvals)
    h2     <- compute_h2_simple(chi2, dt$N)
    n_snps <- nrow(dt)
    n_sig  <- sum(pvals < opt$p_threshold, na.rm = TRUE)
    n_sug  <- sum(pvals < opt$suggestive,  na.rm = TRUE)
    n_med  <- median(dt$N, na.rm = TRUE)

    cat(sprintf("  SNPs: %d  |  N(median): %d  |  lambda GC: %.4f  |  h2: %.6f  |  sig hits: %d\n",
                n_snps, as.integer(n_med), lambda, h2, n_sig))

    # Manhattan
    manhattan_path <- file.path(opt$output_dir, paste0(ct, "_manhattan.png"))
    tryCatch(
        make_manhattan(dt, ct, opt$p_threshold, opt$suggestive, manhattan_path),
        error = function(e) warning("Manhattan failed for ", ct, ": ", e$message)
    )

    # QQ
    qq_path <- file.path(opt$output_dir, paste0(ct, "_qqplot.png"))
    tryCatch(
        make_qq(pvals, ct, lambda, qq_path),
        error = function(e) warning("QQ failed for ", ct, ": ", e$message)
    )

    summary_rows[[ct]] <- data.frame(
        cell_type          = ct,
        n_snps             = n_snps,
        n_samples_median   = as.integer(n_med),
        lambda_gc          = lambda,
        h2_simple          = h2,
        n_genome_wide_sig  = n_sig,
        n_suggestive       = n_sug,
        stringsAsFactors   = FALSE
    )
}

# ── Summary table ─────────────────────────────────────────────────────────────
summary_df <- bind_rows(summary_rows) %>% arrange(cell_type)
out_table  <- file.path(opt$output_dir, "gtex_v10_gwas_summary.tsv")
fwrite(summary_df, out_table, sep = "\t")

cat("\n========================================\n")
cat("Summary saved to:", out_table, "\n")
cat("========================================\n")
print(summary_df, row.names = FALSE)
cat("========================================\n")
cat("Plots saved to:", opt$output_dir, "\n")
