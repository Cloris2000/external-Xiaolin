#!/usr/bin/env Rscript
# Forest plot for top SNP in Oligodendrocyte bulk 15-cohort meta-analysis
# Lead variant: chr7:12252119:G:GT (effect allele = GT, insertion)
# Only 3 of 15 cohorts carry this indel (ROSMAP, Mayo, MSBB — all WGS or WGS-imputed).
# No biallelic SNP proxy with consistent direction across these 3 cohorts was found
# in the TMEM106B region; the signal appears specific to the insertion.
# All betas are reported for the GT allele (insertion).
# Meta Effect from METAL is for G allele (Allele1); flipped here to match cohort direction.

library(ggplot2)
library(dplyr)

# ---------- data -------------------------------------------------------
# Cohort-level: from REGENIE step2 (effect allele = GT)
cohorts <- data.frame(
  label  = c("ROSMAP", "Mayo", "MSBB"),
  beta   = c(0.340655,  0.0416694, 0.0274115),
  se     = c(0.0451664, 0.0633103, 0.0734952),
  n      = c(795L,      256L,      250L),
  pval   = c(4.621e-14, 5.104e-01, 7.092e-01),
  is_meta = FALSE,
  stringsAsFactors = FALSE
)

# Meta-analysis (METAL Effect is for G allele; flip sign → GT allele direction)
meta <- data.frame(
  label   = "Meta-analysis\n(15 cohorts)",
  beta    = 0.1973,
  se      = 0.0329,
  n       = NA_integer_,
  pval    = 1.976e-09,
  is_meta = TRUE,
  stringsAsFactors = FALSE
)

df <- bind_rows(cohorts, meta)
df$label <- factor(df$label, levels = rev(c("ROSMAP", "Mayo", "MSBB",
                                             "Meta-analysis\n(15 cohorts)")))

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
clr_cohort <- "#2166AC"   # blue
clr_meta   <- "#B2182B"   # red

# ---------- plot -------------------------------------------------------
p <- ggplot(df, aes(y = label, x = beta)) +

  # zero line
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +

  # CI segments
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi,
        colour = is_meta, linewidth = is_meta),
    height = 0.25
  ) +

  # point / diamond via shape
  geom_point(
    aes(colour = is_meta,
        shape  = is_meta,
        size   = is_meta)
  ) +

  # sample-size annotation (left of CI)
  geom_text(
    aes(x = min(ci_lo) - 0.04, label = nlabel),
    hjust = 1, size = 5, colour = "grey30"
  ) +

  # p-value annotation (right of CI)
  geom_text(
    aes(x = max(ci_hi) + 0.04, label = paste0("P = ", plabel)),
    hjust = 0, size = 5, colour = "grey20"
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
    title    = "Forest plot — chr7:12252119:G:GT (TMEM106B locus)",
    subtitle = "Oligodendrocyte bulk GWAS · effect allele: GT (insertion) · 3 of 15 WGS cohorts genotyped · no biallelic proxy with consistent direction found",
    x        = "Beta (95% CI)",
    y        = NULL
  ) +

  # expand x to give room for annotations
  scale_x_continuous(
    expand = expansion(add = c(0.30, 0.38))
  ) +

  theme_classic(base_size = 18) +
  theme(
    plot.title    = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 15, colour = "grey40"),
    axis.text.y   = element_text(size = 16),
    axis.text.x   = element_text(size = 14),
    panel.grid.major.x = element_line(colour = "grey90"),
    plot.margin   = margin(20, 120, 20, 20)
  )

# ---------- save -------------------------------------------------------
out <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts/plots/forest/Oligo_TMEM106B_top_SNP_forest.png"
ggsave(out, plot = p, width = 12, height = 5, dpi = 200, bg = "white")
message("Saved: ", out)
