#!/usr/bin/env Rscript
# Forest plot for TMEM106B locus in VIP 15-cohort meta-analysis
# SNP: chr7:12284378:G:A (52 bp from lead indel chr7:12284430:T:TA)
# Effect allele: A (ALT). Present in 11/15 cohorts.
# Missing from GTEx_v10, AMP_AD_Mayo, AMP_AD_Rush (overlapping indel at 12284365),
# and ROSMAP_array (SNP array). NIMH HBCC cohorts are predominantly African-American;
# different LD structure at TMEM106B accounts for heterogeneity in those cohorts.
# METAL Allele1 = a (A allele), Effect = +0.1178; cohort betas also for A allele.

library(ggplot2)
library(dplyr)

# ---------- per-cohort data (REGENIE, effect allele = A) ---------------
# ancestry: "EUR" = European WGS, "AFR" = African-American array/imputed
cohorts <- data.frame(
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

# ---------- meta (METAL, 11 cohorts, Allele1 = A) ----------------------
meta <- data.frame(
  label    = "Meta-analysis",
  beta     =  0.1178,
  se       =  0.0220,
  n        = NA_integer_,
  pval     = 8.831e-08,
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
    # map ancestry to a display label for the legend
    anc_label = dplyr::case_when(
      ancestry == "EUR"  ~ "European",
      ancestry == "AFR"  ~ "Mixed ancestries",
      ancestry == "META" ~ "Meta-analysis"
    )
  )
df$anc_label <- factor(df$anc_label,
                       levels = c("European", "Mixed ancestries", "Meta-analysis"))

# ---------- colours ----------------------------------------------------
clr_eur  <- "#2166AC"
clr_afr  <- "#D95F02"
clr_meta <- "#B2182B"

# ---------- plot -------------------------------------------------------
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
    aes(x = min(ci_lo, na.rm = TRUE) - 0.04, label = nlabel),
    hjust = 1, size = 5.2, colour = "grey35"
  ) +

  geom_text(
    aes(x = max(ci_hi, na.rm = TRUE) + 0.04, label = paste0("P = ", plabel)),
    hjust = 0, size = 5.2, colour = "grey20"
  ) +

  scale_colour_manual(
    name   = "Ancestry",
    values = c("European" = clr_eur,
               "Mixed ancestries" = clr_afr,
               "Meta-analysis" = clr_meta),
    guide  = guide_legend(override.aes = list(shape = 16, size = 4,
                                              linewidth = 0.8))
  ) +
  scale_linewidth_manual(values = c("FALSE" = 0.8, "TRUE" = 1.5), guide = "none") +
  scale_shape_manual(values    = c("FALSE" = 16, "TRUE" = 18),    guide = "none") +
  scale_size_manual(values     = c("FALSE" = 3.5, "TRUE" = 5.5),  guide = "none") +

  labs(
    title = "chr7:12284378:G:A — TMEM106B locus (VIP)",
    x     = "Beta (95% CI), effect allele A",
    y     = NULL
  ) +

  scale_x_continuous(expand = expansion(add = c(0.42, 0.36))) +

  theme_classic(base_size = 20) +
  theme(
    plot.title         = element_text(face = "bold", size = 21, hjust = 0),
    axis.text.y        = element_text(size = 17, colour = "black"),
    axis.text.x        = element_text(size = 16, colour = "black"),
    axis.title.x       = element_text(size = 17),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
    legend.position    = "right",
    legend.background  = element_rect(fill = "white", colour = "grey80",
                                      linewidth = 0.4),
    legend.title       = element_text(size = 16, face = "bold"),
    legend.text        = element_text(size = 15),
    legend.key.size    = unit(1.1, "lines"),
    plot.margin        = margin(16, 150, 16, 16)
  )

# ---------- save -------------------------------------------------------
out <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts/plots/forest/VIP_TMEM106B_chr7_12284378_forest.png"
ggsave(out, plot = p, width = 15, height = 8, dpi = 200, bg = "white")
message("Saved: ", out)
