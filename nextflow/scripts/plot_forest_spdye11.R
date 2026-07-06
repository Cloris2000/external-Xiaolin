#!/usr/bin/env Rscript
# Forest plot for top SNP near SPDYE11 in Oligodendrocyte 15-cohort meta-analysis
# SNP: chr7:73043561:G:A (effect allele = A)
# 9 of 15 cohorts contributed; all show consistent negative direction.
# Meta Effect from METAL is for A allele (minor allele).

library(ggplot2)
library(dplyr)

# ---------- per-cohort data (from REGENIE, effect allele = A) ----------
cohorts <- data.frame(
  label   = c("CMC_MSSM", "CMC_PENN", "CMC_PITT", "GVEX",
              "MSBB", "Mayo", "NABEC", "NIMH_HBCC_Omni5M", "ROSMAP"),
  beta    = c(-0.459541, -0.364685, -0.36499, -0.32809,
              -0.220253, -0.21288, -0.11843, -0.558766, -0.0349315),
  se      = c(0.123458, 0.193801, 0.164304, 0.125305,
              0.114112, 0.122906, 0.142118, 0.197277, 0.088232),
  n       = c(242L, 92L, 161L, 394L, 252L, 257L, 210L, 79L, 793L),
  log10p  = c(3.70451, 1.22279, 1.57969, 2.05372,
              1.27093, 1.07955, 0.392908, 2.33535, 0.159784),
  is_meta = FALSE,
  stringsAsFactors = FALSE
)
cohorts$pval <- 10^(-cohorts$log10p)

# ---------- meta-analysis summary (METAL, effect allele = A) -----------
meta <- data.frame(
  label   = "Meta-analysis\n(15 cohorts)",
  beta    = -0.2409,
  se      =  0.0430,
  n       = NA_integer_,
  log10p  = NA_real_,
  pval    = 2.195e-08,
  is_meta = TRUE,
  stringsAsFactors = FALSE
)

df <- bind_rows(cohorts, meta)

# Order: cohorts top-to-bottom by |beta|, then meta at bottom
cohort_order <- cohorts$label[order(abs(cohorts$beta), decreasing = TRUE)]
df$label <- factor(df$label,
                   levels = rev(c(cohort_order, "Meta-analysis\n(15 cohorts)")))

df <- df %>%
  mutate(
    ci_lo  = beta - 1.96 * se,
    ci_hi  = beta + 1.96 * se,
    plabel = ifelse(pval < 0.001,
                    formatC(pval, format = "e", digits = 2),
                    formatC(pval, format = "f", digits = 3)),
    nlabel = ifelse(is.na(n), "", paste0("N = ", formatC(n, format = "d", big.mark = ",")))
  )

# ---------- colours ----------------------------------------------------
clr_cohort <- "#2166AC"
clr_meta   <- "#B2182B"

# ---------- plot -------------------------------------------------------
x_left  <- min(df$ci_lo, na.rm = TRUE) - 0.10
x_right <- max(df$ci_hi, na.rm = TRUE) + 0.10

p <- ggplot(df, aes(y = label, x = beta)) +

  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +

  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi,
        colour    = is_meta,
        linewidth = is_meta),
    height = 0.25
  ) +

  geom_point(
    aes(colour = is_meta,
        shape  = is_meta,
        size   = is_meta)
  ) +

  # N annotation left
  geom_text(
    aes(x = x_left - 0.04, label = nlabel),
    hjust = 1, size = 4.5, colour = "grey30"
  ) +

  # P-value annotation right
  geom_text(
    aes(x = x_right + 0.04, label = paste0("P = ", plabel)),
    hjust = 0, size = 4.5, colour = "grey20"
  ) +

  scale_colour_manual(values = c("FALSE" = clr_cohort, "TRUE" = clr_meta),
                      guide = "none") +
  scale_linewidth_manual(values = c("FALSE" = 0.8, "TRUE" = 1.4),
                         guide = "none") +
  scale_shape_manual(values  = c("FALSE" = 16, "TRUE" = 18),
                     guide = "none") +
  scale_size_manual(values   = c("FALSE" = 4,  "TRUE" = 6),
                    guide = "none") +

  labs(
    title    = "Forest plot — chr7:73043561:G:A (SPDYE11 locus)",
    subtitle = "Oligodendrocyte bulk GWAS · effect allele: A · 9 of 15 cohorts genotyped",
    x        = "Beta (95% CI)",
    y        = NULL
  ) +

  scale_x_continuous(expand = expansion(add = c(0.45, 0.45))) +

  theme_classic(base_size = 18) +
  theme(
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 15, colour = "grey40"),
    axis.text.y        = element_text(size = 16),
    axis.text.x        = element_text(size = 14),
    panel.grid.major.x = element_line(colour = "grey90"),
    plot.margin        = margin(20, 140, 20, 20)
  )

# ---------- save -------------------------------------------------------
out <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts/plots/forest/Oligo_SPDYE11_top_SNP_forest.png"
ggsave(out, plot = p, width = 13, height = 7, dpi = 200, bg = "white")
message("Saved: ", out)
