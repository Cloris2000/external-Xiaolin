#!/usr/bin/env Rscript
# Figure 5: snRNA-seq CTP GWAS concordance analysis
# Panel A: Direction concordance (same vs opposite) among loci found in sn — all 19 CTs
# Panel B: Bulk vs sn beta (95% CI) at the strongest bulk lead found in sn — all 19 CTs
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
CONC_FILE <- file.path(
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results",
  "sn_bulk_meta_similarity_hodge3/top_hits/cell_type_concordance_summary.tsv"
)
HITS_FILE <- file.path(
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results",
  "sn_bulk_meta_similarity_hodge3/top_hits/bulk_suggestive_hits_sn_direction.tsv"
)

CLR_SAME <- "#3C5488"
CLR_OPP  <- "#E64B35"
CLR_BULK <- "#2166AC"
CLR_SN   <- "#D95F02"
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
# Panel B — Bulk vs sn beta at top lead found in sn, all 19 cell types
# ─────────────────────────────────────────────────────────────────────────────
cat("--- Building Panel B (all cell types) ---\n")
hits_dt <- fread(HITS_FILE)

# Nearest-gene labels (topr); fallback to chr:pos
nearest_gene_labels <- function(chr, pos, p) {
  ch <- sub("^chr", "", as.character(chr))
  fallback <- paste0(ch, ":", format(as.integer(pos), scientific = FALSE, trim = TRUE))
  out <- tryCatch({
    ann <- topr::annotate_with_nearest_gene(
      data.frame(CHROM = ch, POS = as.integer(pos), P = as.numeric(p))
    )
    g <- if ("Gene_Symbol" %in% names(ann)) as.character(ann$Gene_Symbol) else NULL
    if (is.null(g) || length(g) != length(fallback)) fallback else g
  }, error = function(e) fallback)
  out[is.na(out) | out == ""] <- fallback[is.na(out) | out == ""]
  out
}

# Strongest bulk lead that is also present in sn, per cell type
focus_loci <- hits_dt[sn_found == TRUE][
  order(bulk_p), .SD[1], by = bulk_cell_type
]
if (!nrow(focus_loci)) stop("No sn-found loci in ", HITS_FILE)

# Keep panel A cell-type order
focus_loci <- focus_loci[match(ct_order, bulk_cell_type)]
focus_loci <- focus_loci[!is.na(bulk_cell_type)]

focus_loci[, `:=`(
  cell_type = bulk_cell_type,
  gene      = nearest_gene_labels(chrom, pos, bulk_p),
  bulk_beta = as.numeric(bulk_effect),
  bulk_se   = as.numeric(bulk_stderr),
  sn_beta   = as.numeric(sn_effect),
  sn_se     = as.numeric(sn_stderr)
)]
focus_loci[, bulk_lo := bulk_beta - 1.96 * bulk_se]
focus_loci[, bulk_hi := bulk_beta + 1.96 * bulk_se]
focus_loci[, sn_lo   := sn_beta   - 1.96 * sn_se]
focus_loci[, sn_hi   := sn_beta   + 1.96 * sn_se]

focus_loci[, row_label := paste0(cell_type, "  (", gene, ")")]
focus_loci[, row_label := factor(row_label, levels = rev(row_label))]

cat("Panel B loci (top bulk lead found in sn per CT):\n")
print(focus_loci[, .(cell_type, gene, marker, bulk_p, sn_p, concordant)])

long_B <- rbind(
  focus_loci[, .(row_label, cell_type, gwas = "Bulk CTP GWAS",
                 beta = bulk_beta, ci_lo = bulk_lo, ci_hi = bulk_hi)],
  focus_loci[, .(row_label, cell_type, gwas = "sn CTP GWAS",
                 beta = sn_beta, ci_lo = sn_lo, ci_hi = sn_hi)]
)
long_B[, gwas := factor(gwas, levels = c("Bulk CTP GWAS", "sn CTP GWAS"))]
long_B[, y_off := ifelse(gwas == "Bulk CTP GWAS", 0.18, -0.18)]
long_B[, y_num := as.numeric(row_label) + y_off]

pB <- ggplot(long_B, aes(x = beta, y = y_num, colour = gwas, shape = gwas)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.10, linewidth = 0.55) +
  geom_point(size = 2.8) +
  scale_colour_manual(
    values = c("Bulk CTP GWAS" = CLR_BULK, "sn CTP GWAS" = CLR_SN),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Bulk CTP GWAS" = 16, "sn CTP GWAS" = 17),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(
    breaks = as.numeric(focus_loci$row_label),
    labels = levels(focus_loci$row_label)
  ) +
  labs(x = expression(beta ~ "(95% CI)"), y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_blank(),
    axis.text.y        = element_text(size = 8.5, colour = "black"),
    axis.text.x        = element_text(size = 9),
    axis.title.x       = element_text(size = 10),
    legend.position    = "bottom",
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(1.0, "lines"),
    legend.margin      = margin(t = -2),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
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
  rel_heights = c(1.0, 1.35),
  labels      = c("A", "B"),
  label_size  = 14,
  label_fontface = "bold"
)

out_png <- file.path(OUT_DIR, "figure5_sn_concordance.png")
out_pdf <- file.path(OUT_DIR, "figure5_sn_concordance.pdf")

cowplot::save_plot(out_png, fig5,
                   base_width = 12, base_height = 16,
                   dpi = FIG_DPI, bg = "white")
cowplot::save_plot(out_pdf, fig5,
                   base_width = 12, base_height = 16,
                   bg = "white")

cat("\n=== Figure 5 outputs ===\n")
cat("  PNG:", out_png, "\n")
cat("  PDF:", out_pdf, "\n")
cat("Done.\n")
