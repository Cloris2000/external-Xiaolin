#!/usr/bin/env Rscript
# Forest plot for CACNA1C locus in L5.6.IT.Car3 15-cohort meta-analysis
# SNP: chr12:2324042:T:C  (METAL Allele1 = T; effect allele = T throughout)
# Present in 11/15 cohorts (Direction: ??++-?+++-++++?)
# Missing from: AMP_AD_Mayo, AMP_AD_Rush (not in meta input), GTEx_v10 (?),
#   ROSMAP_array (SNP array coverage).
# METAL Effect = +0.1018 for T allele; all per-cohort betas are flipped from
#   REGENIE output (which uses C=ALLELE1) so that every entry = beta for T allele.
# NIMH HBCC cohorts are predominantly African-American.

library(ggplot2)
library(dplyr)

# ---------- per-cohort data (REGENIE ALLELE1=C, flipped → T-allele betas) ----
# ancestry: "EUR" = European, "AFR" = Mixed ancestries (predominantly AFR)
cohorts <- data.frame(
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
               0.00812, 0.8916,  0.8216,  0.1040,
               0.002110, 0.1336,  0.1879),
  ancestry = c(rep("EUR", 8), rep("AFR", 3)),
  is_meta  = FALSE,
  stringsAsFactors = FALSE
)

# ---------- meta-analysis (METAL, 11 cohorts, Allele1 = T) -------------------
meta <- data.frame(
  label    = "Meta-analysis",
  beta     =  0.1018,
  se       =  0.0221,
  n        = NA_integer_,
  pval     = 3.959e-06,
  ancestry = "META",
  is_meta  = TRUE,
  stringsAsFactors = FALSE
)

df <- bind_rows(cohorts, meta)

cohort_order <- c("ROSMAP", "Mayo", "MSBB", "NABEC",
                  "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
                  "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
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

# ---------- colours ----------------------------------------------------------
clr_eur  <- "#2166AC"
clr_afr  <- "#D95F02"
clr_meta <- "#B2182B"

# ---------- plot -------------------------------------------------------------
p <- ggplot(df, aes(y = label, x = beta, colour = anc_label)) +

  # dotted separator between EUR and AFR cohort groups
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
    guide  = guide_legend(override.aes = list(shape = 16, size = 3,
                                              linewidth = 0.6))
  ) +
  scale_linewidth_manual(values = c("FALSE" = 0.7, "TRUE" = 1.3), guide = "none") +
  scale_shape_manual(values    = c("FALSE" = 16, "TRUE" = 18),    guide = "none") +
  scale_size_manual(values     = c("FALSE" = 2.5, "TRUE" = 4.0),  guide = "none") +

  labs(
    title = "chr12:2324042:T:C — CACNA1C (L5.6.IT.Car3)",
    x     = "Beta (95% CI), effect allele T",
    y     = NULL
  ) +

  scale_x_continuous(expand = expansion(add = c(0.26, 0.42))) +

  theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y        = element_text(size = 10, colour = "black"),
    axis.text.x        = element_text(size = 9,  colour = "black"),
    axis.title.x       = element_text(size = 10),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "none",
    legend.background  = element_rect(fill = "white", colour = "grey80",
                                      linewidth = 0.3),
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(0.8, "lines"),
    plot.margin        = margin(8, 100, 8, 8)
  )

# ---------- save -------------------------------------------------------------
out_dir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts/plots/forest"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "L5.6.IT.Car3_CACNA1C_chr12_2324042_forest.png")
ggsave(out, plot = p, width = 6, height = 8, dpi = 200, bg = "white")
message("Saved: ", out)
