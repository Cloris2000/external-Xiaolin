#!/usr/bin/env Rscript
# Figure 5: snRNA-seq CTP GWAS concordance analysis
# Panel A: Direction concordance (same vs opposite) among loci found in sn — all 19 CTs
# Panel B: Bulk β vs sn β scatter over all bulk-suggestive loci found in sn;
#          dot size ∝ -log10(bulk P), coloured by direction, with Spearman r + n
#
# Outputs:
#   manuscript_figure/figure5_sn_concordance.png
#   manuscript_figure/figure5_sn_concordance.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

OUT_DIR   <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure"

# Which sn-vs-bulk concordance analysis to plot.
#   "hodge5u" = 5-cohort sn meta, de-duplicated (ROSMAP_Green + PsychAD_HBCC
#               + PsychAD_MSSM + Ruz_MSSM/Ruzicka + ROSMAP_Mathys_unique, i.e.
#               the 128 Mathys donors NOT shared with ROSMAP_Green)  [current]
#   "hodge5m" = 5-cohort sn meta with the FULL ROSMAP_Mathys_Hodge (323 donors;
#               194 shared with Green -> double-counted, superseded by hodge5u)
#   "hodge5"  = 5-cohort sn meta using the PsychAD_RADC_snRNA stand-in (89 donors)
#   "hodge3"  = original 3-cohort sn meta (ROSMAP_Green + PsychAD_HBCC + PsychAD_MSSM)
ANALYSIS_TAG <- "hodge5u"
CONC_FILE <- file.path(
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results",
  paste0("sn_bulk_meta_similarity_", ANALYSIS_TAG,
         "/top_hits/cell_type_concordance_summary.tsv")
)
HITS_FILE <- file.path(
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results",
  paste0("sn_bulk_meta_similarity_", ANALYSIS_TAG,
         "/top_hits/bulk_suggestive_hits_sn_direction.tsv")
)

# Panel A palette — cnsplots 'Nature' scheme (navy / red)
CLR_SAME <- "#3C5488"
CLR_OPP  <- "#E64B35"
CLR_BULK <- "#2166AC"
CLR_SN   <- "#D95F02"

# Panel B palette — cnsplots 'Cell' scheme (teal / amber), deliberately
# distinct from Panel A so the two panels don't share colours.
CLR_B_SAME <- "#2F7E8F"   # cnsplots Cell teal
CLR_B_OPP  <- "#E1A22E"   # cnsplots Cell amber

FIG_DPI  <- 200

# ─────────────────────────────────────────────────────────────────────────────
# Panel A — Same vs opposite only (no "not in sn" grey bars)
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Building Panel A ---\n")
conc_dt <- fread(CONC_FILE)
conc_dt[, pct_of_found := ifelse(n_found_sn > 0,
                                 round(n_concordant / n_found_sn * 100), NA_real_)]
setorder(conc_dt, -pct_of_found)
ct_order <- conc_dt$bulk_cell_type

bar_dt <- melt(
  conc_dt[, .(bulk_cell_type,
              `Same direction` = n_concordant,
              `Opposite`       = n_discordant)],
  id.vars       = "bulk_cell_type",
  variable.name = "category",
  value.name    = "n"
)
bar_dt[, bulk_cell_type := factor(bulk_cell_type, levels = rev(ct_order))]
bar_dt[, category := factor(category, levels = c("Opposite", "Same direction"))]

pct_lbl <- conc_dt[n_found_sn > 0,
                   .(bulk_cell_type,
                     lbl   = paste0(pct_of_found, "%"),
                     total = n_found_sn)]
pct_lbl[, bulk_cell_type := factor(bulk_cell_type, levels = rev(ct_order))]

pA <- ggplot(bar_dt, aes(x = n, y = bulk_cell_type, fill = category)) +
  geom_col(width = 0.7, position = position_stack()) +
  geom_text(
    data = pct_lbl,
    aes(x = total, y = bulk_cell_type, label = lbl),
    inherit.aes = FALSE,
    hjust = -0.15, size = 3.2, colour = "grey30"
  ) +
  scale_fill_manual(
    values = c("Same direction" = CLR_SAME, "Opposite" = CLR_OPP),
    name = NULL,
    guide = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of loci", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title           = element_blank(),
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.text          = element_text(size = 10),
    legend.key.size      = unit(0.8, "lines"),
    legend.spacing.x     = unit(12, "pt"),
    legend.margin        = margin(t = 4),
    axis.text.y          = element_text(size = 9, colour = "black"),
    axis.text.x          = element_text(size = 9),
    axis.title.x         = element_text(size = 10),
    plot.margin          = margin(6, 8, 6, 8)
  )

cat("Panel A built.\n")

# ─────────────────────────────────────────────────────────────────────────────
# Panel B — Bulk β vs sn β scatter across all bulk-suggestive loci found in sn.
#   x = bulk CTP-GWAS β, y = sn CTP-GWAS β, dot size ∝ -log10(bulk P),
#   coloured by effect-direction concordance. Annotated with overall Spearman r,
#   n loci, and % same-direction.
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Building Panel B (bulk vs sn beta scatter) ---\n")
hits_dt <- fread(HITS_FILE)

scat <- hits_dt[sn_found == TRUE]
scat[, `:=`(
  bulk_beta = as.numeric(bulk_effect),
  sn_beta   = as.numeric(sn_effect),
  bulk_p    = as.numeric(bulk_p)
)]
scat <- scat[is.finite(bulk_beta) & is.finite(sn_beta) & is.finite(bulk_p)]
scat[, neg_log10_bulk_p := -log10(bulk_p)]
scat[, direction := ifelse(sign(bulk_beta) == sign(sn_beta),
                           "Same direction", "Opposite")]
scat[, direction := factor(direction, levels = c("Same direction", "Opposite"))]

# Overall summary stats
n_loci   <- nrow(scat)
rho      <- suppressWarnings(cor(scat$bulk_beta, scat$sn_beta, method = "spearman"))
rho_p    <- suppressWarnings(
  cor.test(scat$bulk_beta, scat$sn_beta, method = "spearman")$p.value)
pct_same <- 100 * mean(scat$direction == "Same direction")

cat(sprintf("Panel B: n=%d loci, Spearman rho=%.3f (p=%.2g), same-direction=%.1f%%\n",
            n_loci, rho, rho_p, pct_same))

rho_p_lbl <- if (rho_p < 2.2e-16) {
  "P < 2.2e-16"
} else {
  paste0("P = ", formatC(rho_p, format = "g", digits = 2))
}
annot_lbl <- paste0(
  "Spearman r = ", sprintf("%.2f", rho), "\n",
  rho_p_lbl, "\n",
  "n = ", n_loci, " loci\n",
  "Same direction: ", sprintf("%.0f%%", pct_same)
)

# Symmetric axis limits so the y = x diagonal is meaningful
lim <- max(abs(c(scat$bulk_beta, scat$sn_beta)), na.rm = TRUE) * 1.05

pB <- ggplot(scat, aes(x = bulk_beta, y = sn_beta)) +
  # shade the two concordant (same-sign) quadrants
  annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim,
           fill = CLR_B_SAME, alpha = 0.06) +
  annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0,
           fill = CLR_B_SAME, alpha = 0.06) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey45", linewidth = 0.5) +
  geom_point(aes(size = neg_log10_bulk_p, fill = direction),
             shape = 21, colour = "white", stroke = 0.25, alpha = 0.85) +
  scale_fill_manual(
    values = c("Same direction" = CLR_B_SAME, "Opposite" = CLR_B_OPP),
    name = NULL,
    guide = guide_legend(order = 1, override.aes = list(size = 4))
  ) +
  scale_size_continuous(
    name = expression(-log[10] * "(" * P[bulk] * ")"),
    range = c(1.4, 7), guide = guide_legend(order = 2)
  ) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  annotate("text", x = -lim * 0.97, y = lim * 0.97,
           label = annot_lbl, hjust = 0, vjust = 1,
           size = 3.5, colour = "grey15", lineheight = 0.95) +
  labs(x = expression("Bulk CTP GWAS " * beta),
       y = expression("sn CTP GWAS " * beta)) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_blank(),
    axis.text          = element_text(size = 9, colour = "black"),
    axis.title         = element_text(size = 10),
    legend.position    = "right",
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 9),
    legend.key.size    = unit(0.9, "lines"),
    panel.grid.major   = element_line(colour = "grey94", linewidth = 0.3),
    plot.margin        = margin(6, 10, 6, 6)
  )

cat("Panel B built.\n")

# ─────────────────────────────────────────────────────────────────────────────
# Assemble + save
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Assembling Figure 5 ---\n")
fig5 <- plot_grid(
  pA, pB,
  ncol        = 1,
  rel_heights = c(1.0, 1.15),
  labels      = c("A", "B"),
  label_size  = 14,
  label_fontface = "bold"
)

out_png <- file.path(OUT_DIR, "figure5_sn_concordance.png")
out_pdf <- file.path(OUT_DIR, "figure5_sn_concordance.pdf")

cowplot::save_plot(out_png, fig5,
                   base_width = 10, base_height = 13,
                   dpi = FIG_DPI, bg = "white")
cowplot::save_plot(out_pdf, fig5,
                   base_width = 10, base_height = 13,
                   bg = "white")

cat("\n=== Figure 5 outputs ===\n")
cat("  PNG:", out_png, "\n")
cat("  PDF:", out_pdf, "\n")
cat("Done.\n")
