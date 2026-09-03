#!/usr/bin/env Rscript
# Forest plot for chr6:164862615 locus in Microglia 15-cohort meta-analysis
# SNP: chr6:164862615:T:C  (METAL Allele1 = T; effect allele = T throughout)
# Present in 11/15 cohorts; AMP_AD_Mayo, AMP_AD_Rush, GTEx_v10, ROSMAP_array: SNP absent in cohort GWAS
# REGENIE output: ALLELE1 = C; betas are flipped to T-allele convention.
# Meta-analysis Direction: --+++-++++++++?  (effect allele T, Effect = +0.0619)

library(ggplot2)
library(dplyr)

# ── per-cohort data (REGENIE ALLELE1=C, flipped → T-allele betas) ─────────────
cohorts <- data.frame(
  label    = c("CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
               "MSBB", "Mayo", "NABEC",
               "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
               "ROSMAP"),
  beta     = c( 0.10286,  0.12062,  0.02824,  0.14401,
                0.07872,  0.08979,  0.02090,
                0.12192,  0.15638,  0.03106,
               -0.00811),
  se       = c( 0.08206,  0.11202,  0.08396,  0.07213,
                0.07070,  0.07317,  0.07921,
                0.09186,  0.11578,  0.10052,
                0.04728),
  n        = c(245L, 94L, 166L, 394L,
               225L, 256L, 210L,
               202L,  79L,  97L,
               825L),
  pval     = c(0.2098, 0.2815, 0.7358, 0.0458,
               0.2655, 0.2198, 0.7913,
               0.1843, 0.1769, 0.7572,
               0.8631),
  ancestry = c(rep("EUR", 7), rep("AFR", 3), "EUR"),
  is_meta  = FALSE,
  stringsAsFactors = FALSE
)

# ── meta-analysis (METAL, 11 visible cohorts + others, Allele1 = T) ───────────
meta <- data.frame(
  label    = "Meta-analysis",
  beta     =  0.0619,
  se       =  0.0204,
  n        = NA_integer_,
  pval     = 0.002393,
  ancestry = "META",
  is_meta  = TRUE,
  stringsAsFactors = FALSE
)

df <- bind_rows(cohorts, meta)

cohort_order <- c("CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
                  "MSBB", "Mayo", "NABEC",
                  "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
                  "ROSMAP",
                  "Meta-analysis")
df$label <- factor(df$label, levels = rev(cohort_order))

df <- df %>%
  mutate(
    ci_lo     = beta - 1.96 * se,
    ci_hi     = beta + 1.96 * se,
    plabel    = ifelse(pval < 0.001,
                       formatC(pval, format = "e", digits = 2),
                       formatC(pval, format = "f", digits = 3)),
    nlabel    = ifelse(is.na(n), "", paste0("N = ", formatC(n, format = "d", big.mark = ","))),
    anc_label = dplyr::case_when(
      ancestry == "EUR"  ~ "European",
      ancestry == "AFR"  ~ "Mixed ancestries",
      ancestry == "META" ~ "Meta-analysis"
    )
  )
df$anc_label <- factor(df$anc_label,
                       levels = c("European", "Mixed ancestries", "Meta-analysis"))

# ── colours ───────────────────────────────────────────────────────────────────
clr_eur  <- "#2166AC"
clr_afr  <- "#D95F02"
clr_meta <- "#B2182B"

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(y = label, x = beta, colour = anc_label)) +

  geom_hline(yintercept = 3.5, linetype = "dotted", colour = "grey65", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +

  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi, linewidth = is_meta),
    height = 0.28
  ) +
  geom_point(aes(shape = is_meta, size = is_meta)) +

  geom_text(
    aes(x = min(ci_lo, na.rm = TRUE) - 0.02, label = nlabel),
    hjust = 1, size = 3.2, colour = "grey35"
  ) +
  geom_text(
    aes(x = max(ci_hi, na.rm = TRUE) + 0.02, label = paste0("P = ", plabel)),
    hjust = 0, size = 3.2, colour = "grey20"
  ) +

  scale_colour_manual(
    name   = "Ancestry",
    values = c("European"         = clr_eur,
               "Mixed ancestries" = clr_afr,
               "Meta-analysis"    = clr_meta),
    guide  = guide_legend(override.aes = list(shape = 16, size = 3, linewidth = 0.6))
  ) +
  scale_linewidth_manual(values = c("FALSE" = 0.7, "TRUE" = 1.3), guide = "none") +
  scale_shape_manual(values    = c("FALSE" = 16, "TRUE" = 18),    guide = "none") +
  scale_size_manual(values     = c("FALSE" = 2.5, "TRUE" = 4.0),  guide = "none") +

  labs(
    title = "chr6:164862615:T:C — PRKN locus (Microglia)",
    x     = "Beta (95% CI), effect allele T",
    y     = NULL
  ) +

  scale_x_continuous(expand = expansion(add = c(0.24, 0.40))) +

  theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y        = element_text(size = 10, colour = "black"),
    axis.text.x        = element_text(size = 9,  colour = "black"),
    axis.title.x       = element_text(size = 10),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "none",
    legend.background  = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(0.8, "lines"),
    plot.margin        = margin(8, 100, 8, 8)
  )

# ── save ──────────────────────────────────────────────────────────────────────
out_dir <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/meta_analysis_15cohorts/plots/forest"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "Microglia_PRKN_chr6_164862615_forest.png")
ggsave(out, plot = p, width = 6, height = 8, dpi = 200, bg = "white")
message("Saved: ", out)
