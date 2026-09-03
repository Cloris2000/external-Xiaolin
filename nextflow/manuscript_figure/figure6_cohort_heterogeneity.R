#!/usr/bin/env Rscript
# Figure 6: Cross-cohort consistency of top CTP-GWAS loci
# Three forest plots (VIP × TMEM106B, L5.6.IT.Car3 × CACNA1C, Microglia × chr6)
# each annotated with Cochran's Q and I² heterogeneity statistics.
#
# Outputs:
#   manuscript_figure/figure6_cohort_heterogeneity.png
#   manuscript_figure/figure6_cohort_heterogeneity.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

OUT_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"

# Colours (consistent with other figures)
CLR_EUR  <- "#2166AC"   # European
CLR_AFR  <- "#D95F02"   # Mixed ancestries (predominantly AFR)
CLR_META <- "#B2182B"   # Meta-analysis diamond

FIG_DPI <- 200

# ─────────────────────────────────────────────────────────────────────────────
# Helper: compute Cochran's Q and I²
# ─────────────────────────────────────────────────────────────────────────────
compute_heterogeneity <- function(beta, se) {
  wi   <- 1 / se^2
  bhat <- sum(wi * beta) / sum(wi)
  Q    <- sum(wi * (beta - bhat)^2)
  df   <- length(beta) - 1L
  I2   <- max(0, (Q - df) / Q) * 100
  Qp   <- pchisq(Q, df = df, lower.tail = FALSE)
  list(Q = Q, df = df, I2 = I2, Qp = Qp)
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: make one forest plot with I² annotation
# ─────────────────────────────────────────────────────────────────────────────
make_forest <- function(cohorts_df, meta_df, cohort_order, panel_letter,
                        title_str, x_label, sep_yintercept = NULL) {

  df <- bind_rows(cohorts_df, meta_df)
  df$label <- factor(df$label, levels = rev(cohort_order))

  df <- df %>%
    mutate(
      ci_lo    = beta - 1.96 * se,
      ci_hi    = beta + 1.96 * se,
      plabel   = ifelse(pval < 0.001,
                        formatC(pval, format = "e", digits = 2),
                        formatC(pval, format = "f", digits = 3)),
      nlabel   = ifelse(is.na(n), "", paste0("N=", formatC(n, format = "d", big.mark = ","))),
      anc_label = case_when(
        ancestry == "EUR"  ~ "European",
        ancestry == "AFR"  ~ "Mixed ancestries",
        ancestry == "META" ~ "Meta-analysis"
      )
    )
  df$anc_label <- factor(df$anc_label,
                          levels = c("European", "Mixed ancestries", "Meta-analysis"))

  # Heterogeneity (cohorts only)
  coh <- cohorts_df
  het <- compute_heterogeneity(coh$beta, coh$se)
  het_label <- sprintf("I² = %.1f%%", het$I2)

  x_range <- range(c(df$ci_lo, df$ci_hi), na.rm = TRUE)
  x_left  <- x_range[1] - diff(x_range) * 0.05
  x_right <- x_range[2] + diff(x_range) * 0.05

  p <- ggplot(df, aes(y = label, x = beta, colour = anc_label)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.5)

  if (!is.null(sep_yintercept)) {
    p <- p + geom_hline(yintercept = sep_yintercept,
                        linetype = "dotted", colour = "grey65", linewidth = 0.5)
  }

  p <- p +
    geom_errorbarh(
      aes(xmin = ci_lo, xmax = ci_hi, linewidth = is_meta),
      height = 0.28
    ) +
    geom_point(aes(shape = is_meta, size = is_meta)) +

    # N labels (left)
    geom_text(aes(x = x_left - 0.02, label = nlabel),
              hjust = 1, size = 3.0, colour = "grey40") +

    # P-value labels (right)
    geom_text(aes(x = x_right + 0.02, label = paste0("P=", plabel)),
              hjust = 0, size = 3.0, colour = "grey20") +

    # I² annotation (top-right inside panel)
    annotate("text",
             x = Inf, y = Inf,
             label = het_label,
             hjust = 1.05, vjust = 1.3,
             size = 3.2, colour = "grey30", fontface = "italic") +

    scale_colour_manual(
      name   = "Ancestry",
      values = c("European"        = CLR_EUR,
                 "Mixed ancestries" = CLR_AFR,
                 "Meta-analysis"   = CLR_META),
      guide  = guide_legend(override.aes = list(shape = 16, size = 3,
                                                linewidth = 0.6))
    ) +
    scale_linewidth_manual(values = c("FALSE" = 0.7, "TRUE" = 1.4), guide = "none") +
    scale_shape_manual(values    = c("FALSE" = 16,  "TRUE" = 18),   guide = "none") +
    scale_size_manual(values     = c("FALSE" = 2.5, "TRUE" = 4.2),  guide = "none") +

    scale_x_continuous(expand = expansion(add = c(0.30, 0.44))) +

    labs(
      title = paste0(panel_letter, ". ", title_str),
      x     = x_label,
      y     = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 12, hjust = 0),
      axis.text.y        = element_text(size = 10, colour = "black"),
      axis.text.x        = element_text(size = 9,  colour = "black"),
      axis.title.x       = element_text(size = 10),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      legend.position    = "none",
      plot.margin        = margin(8, 8, 8, 8)
    )

  list(plot = p, het = het)
}

# ─────────────────────────────────────────────────────────────────────────────
# Panel A — VIP × TMEM106B (chr7:12284378:G:A, effect allele = A)
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Building Panel A (VIP × TMEM106B) ---\n")
vip_cohorts <- data.frame(
  label    = c("ROSMAP", "Mayo", "MSBB", "NABEC",
               "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
               "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650"),
  beta     = c( 0.318855,  0.154605,  0.229885,  0.0346592,
                0.0990067, 0.127878, -0.0901743,  0.00345954,
                0.0280202, -0.0826726, -0.193582),
  se       = c( 0.0451451, 0.0571266, 0.0764044, 0.0886242,
                0.0810084, 0.0894102,  0.09388,   0.067703,
                0.0851864,  0.0949025,  0.0960267),
  n        = c(795L, 257L, 252L, 210L,
               242L,  92L, 161L, 394L,
               202L,  79L,  97L),
  pval     = c(1.63042e-12, 0.00680268, 0.00262289, 0.695738,
               0.22164,     0.152649,   0.33679,    0.959247,
               0.74221,     0.383683,   0.0438087),
  ancestry = c(rep("EUR", 8), rep("AFR", 3)),
  is_meta  = FALSE,
  stringsAsFactors = FALSE
)

vip_meta <- data.frame(
  label = "Meta-analysis", beta = 0.1178, se = 0.0220,
  n = NA_integer_, pval = 8.831e-08, ancestry = "META",
  is_meta = TRUE, stringsAsFactors = FALSE
)

vip_order <- c("ROSMAP", "Mayo", "MSBB", "NABEC",
               "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
               "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
               "Meta-analysis")

res_A <- make_forest(
  vip_cohorts, vip_meta, vip_order,
  panel_letter    = "A",
  title_str       = "VIP × TMEM106B (chr7:12284378:G:A)",
  x_label         = "Beta (95% CI), effect allele A",
  sep_yintercept  = 3.5
)
pA <- res_A$plot
cat(sprintf("  VIP: I² = %.1f%%, Q p = %.4f\n", res_A$het$I2, res_A$het$Qp))

# ─────────────────────────────────────────────────────────────────────────────
# Panel B — L5.6.IT.Car3 × CACNA1C (chr12:2324042:T:C, effect allele = T)
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Building Panel B (L5.6.IT.Car3 × CACNA1C) ---\n")
car3_cohorts <- data.frame(
  label    = c("ROSMAP", "Mayo", "MSBB", "NABEC",
               "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
               "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650"),
  beta     = c( 0.0423,  0.1314,  0.1710, -0.0316,
                0.1854,  0.0163, -0.0203,  0.1114,
                0.2233,  0.1080,  0.1091),
  se       = c( 0.04636, 0.07910, 0.07265, 0.09907,
                0.07003, 0.11999, 0.09023, 0.06855,
                0.07263, 0.07202, 0.08284),
  n        = c(796L, 257L, 252L, 210L,
               242L,  92L, 161L, 394L,
               202L,  79L,  97L),
  pval     = c(0.3616, 0.0968, 0.01856, 0.7499,
               0.00812, 0.8916, 0.8216, 0.1040,
               0.002110, 0.1336, 0.1879),
  ancestry = c(rep("EUR", 8), rep("AFR", 3)),
  is_meta  = FALSE,
  stringsAsFactors = FALSE
)

car3_meta <- data.frame(
  label = "Meta-analysis", beta = 0.1018, se = 0.0221,
  n = NA_integer_, pval = 3.959e-06, ancestry = "META",
  is_meta = TRUE, stringsAsFactors = FALSE
)

car3_order <- c("ROSMAP", "Mayo", "MSBB", "NABEC",
                "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
                "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
                "Meta-analysis")

res_B <- make_forest(
  car3_cohorts, car3_meta, car3_order,
  panel_letter    = "B",
  title_str       = "L5.6.IT.Car3 × CACNA1C (chr12:2324042:T:C)",
  x_label         = "Beta (95% CI), effect allele T",
  sep_yintercept  = 3.5
)
pB <- res_B$plot
cat(sprintf("  Car3: I² = %.1f%%, Q p = %.4f\n", res_B$het$I2, res_B$het$Qp))

# ─────────────────────────────────────────────────────────────────────────────
# Panel C — Microglia × PRKN region (chr6:164862615:T:C, effect allele = T)
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Building Panel C (Microglia × chr6) ---\n")
micro_cohorts <- data.frame(
  label    = c("ROSMAP", "Mayo", "MSBB", "NABEC",
               "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
               "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650"),
  beta     = c(-0.00811,  0.08979,  0.07872,  0.02090,
                0.10286,  0.12062,  0.02824,  0.14401,
                0.12192,  0.15638,  0.03106),
  se       = c( 0.04728,  0.07317,  0.07070,  0.07921,
                0.08206,  0.11202,  0.08396,  0.07213,
                0.09186,  0.11578,  0.10052),
  n        = c(825L, 256L, 225L, 210L,
               245L,  94L, 166L, 394L,
               202L,  79L,  97L),
  pval     = c(0.8631, 0.2198, 0.2655, 0.7913,
               0.2098, 0.2815, 0.7358, 0.0458,
               0.1843, 0.1769, 0.7572),
  ancestry = c("EUR", rep("EUR", 6), rep("AFR", 3), "EUR"),
  is_meta  = FALSE,
  stringsAsFactors = FALSE
)

micro_meta <- data.frame(
  label = "Meta-analysis", beta = 0.0619, se = 0.0204,
  n = NA_integer_, pval = 0.002393, ancestry = "META",
  is_meta = TRUE, stringsAsFactors = FALSE
)

micro_order <- c("ROSMAP", "Mayo", "MSBB", "NABEC",
                 "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
                 "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
                 "Meta-analysis")

res_C <- make_forest(
  micro_cohorts, micro_meta, micro_order,
  panel_letter    = "C",
  title_str       = "Microglia × PRKN region (chr6:164862615:T:C)",
  x_label         = "Beta (95% CI), effect allele T",
  sep_yintercept  = 3.5
)
pC <- res_C$plot
cat(sprintf("  Microglia: I² = %.1f%%, Q p = %.4f\n", res_C$het$I2, res_C$het$Qp))

# ─────────────────────────────────────────────────────────────────────────────
# Shared legend
# ─────────────────────────────────────────────────────────────────────────────
legend_plot <- ggplot(
  data.frame(
    x   = c(0, 0, 0),
    y   = c(1, 2, 3),
    grp = factor(c("European", "Mixed ancestries", "Meta-analysis"),
                 levels = c("European", "Mixed ancestries", "Meta-analysis"))
  ),
  aes(x = x, y = y, colour = grp)
) +
  geom_point(size = 3, shape = 16) +
  scale_colour_manual(
    name   = "Ancestry / stratum",
    values = c("European"         = CLR_EUR,
               "Mixed ancestries" = CLR_AFR,
               "Meta-analysis"    = CLR_META),
    guide  = guide_legend(nrow = 1,
                          override.aes = list(shape = 16, size = 3))
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(1.0, "lines"),
    legend.margin    = margin(0, 0, 0, 0)
  )

shared_leg <- get_legend(legend_plot)

# ─────────────────────────────────────────────────────────────────────────────
# Assemble Figure 6
# Three forest plots side by side with shared legend underneath
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Assembling Figure 6 ---\n")

forest_row <- plot_grid(pA, pB, pC,
                        nrow       = 1,
                        rel_widths = c(1, 1, 1),
                        align      = "h",
                        axis       = "tb")

fig6 <- plot_grid(
  forest_row,
  shared_leg,
  ncol        = 1,
  rel_heights = c(1, 0.08)
)

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Saving Figure 6 ---\n")
out_png <- file.path(OUT_DIR, "figure6_cohort_heterogeneity.png")
out_pdf <- file.path(OUT_DIR, "figure6_cohort_heterogeneity.pdf")

cowplot::save_plot(out_png, fig6,
                   base_width  = 18,
                   base_height = 9,
                   dpi         = FIG_DPI,
                   bg          = "white")

cowplot::save_plot(out_pdf, fig6,
                   base_width  = 18,
                   base_height = 9,
                   bg          = "white")

cat("\n=== Figure 6 outputs ===\n")
cat("  PNG:", out_png, "\n")
cat("  PDF:", out_pdf, "\n")
cat("Done.\n")
