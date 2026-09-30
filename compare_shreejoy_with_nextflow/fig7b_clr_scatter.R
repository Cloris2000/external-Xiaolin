# CLR scatterplots between matched snRNA-seq datasets.
#   fig7b: nine-type shortlist, Green vs Mathys
#   fig7c: Sst_25 across the six independent pairs
#
# Writes under compare_shreejoy_with_nextflow/. Does not edit the wp3 tree.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

OUT_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

dt <- fread(file.path(OUT_DIR, "data/fig7b_clr_scatter.tsv"))

SHORTLIST <- c(
  "Astro_6-SEAAD", "OPC_2_2-SEAAD", "OPC_2",
  "Sst_25", "Sst_3", "Micro-PVM_2_1-SEAAD",
  "Sst_23", "Micro-PVM_2", "Sst_2"
)
dt[, cell_type := factor(cell_type, levels = SHORTLIST)]
dt[, reactive := grepl("SEAAD$", cell_type)]
dt[, key := fcase(
  reactive, "Reactive (-SEAAD)",
  group == "neu", "Neuronal",
  default = "Non-neuronal"
)]
dt[, key := factor(key, levels = c("Neuronal", "Non-neuronal", "Reactive (-SEAAD)"))]

NEU <- "#3C5488"
GLIA <- "#E64B35"
COL <- c("Neuronal" = NEU, "Non-neuronal" = GLIA, "Reactive (-SEAAD)" = GLIA)
fontsize <- 8

theme_cns <- theme_classic(base_size = fontsize, base_family = "sans") +
  theme(
    text = element_text(size = fontsize, colour = "black", family = "sans"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.text = element_text(size = fontsize, colour = "black"),
    axis.title = element_text(size = fontsize, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = fontsize, colour = "black", face = "plain"),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

stats_for <- function(d) {
  d[, {
    ok <- is.finite(clr_a) & is.finite(clr_b)
    va <- if (sum(ok) >= 3) var(clr_a[ok]) else NA_real_
    icc <- if (is.finite(va) && va > 0) cov(clr_a[ok], clr_b[ok]) / va else NA_real_
    list(n = sum(ok), icc = icc)
  }, by = .(cell_type, pair)]
}

square_limits <- function(d) {
  d[, {
    r <- range(c(clr_a, clr_b), na.rm = TRUE)
    pad <- max(diff(r) * 0.06, 0.15)
    lo <- r[1] - pad
    hi <- r[2] + pad
    list(clr_a = c(lo, hi), clr_b = c(lo, hi))
  }, by = .(cell_type, pair)]
}

panel <- function(d, xlab, ylab) {
  st <- stats_for(d)
  st[, lab := sprintf("ICC = %.2f", icc)]
  lims <- square_limits(d)
  ggplot(d, aes(clr_a, clr_b)) +
    geom_blank(data = lims, aes(clr_a, clr_b), inherit.aes = FALSE) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55",
                linewidth = 0.3, linetype = "dashed") +
    geom_point(aes(colour = key), size = 0.9, alpha = 0.55, stroke = 0) +
    geom_text(data = st, aes(label = lab),
              x = -Inf, y = Inf, hjust = -0.08, vjust = 1.25,
              size = fontsize / ggplot2::.pt,
              family = "sans", colour = "black", inherit.aes = FALSE) +
    scale_colour_manual(values = COL) +
    labs(x = xlab, y = ylab) +
    theme_cns +
    theme(aspect.ratio = 1)
}

# ---- Green x Mathys, nine-type shortlist -----------------------------------
gm <- dt[pair == "Green x Mathys" & !is.na(cell_type)]
p_b <- panel(gm,
             "Centred log-ratio of cell-type proportion (Green, ln)",
             "Centred log-ratio of cell-type proportion (Mathys, ln)") +
  facet_wrap(~ cell_type, nrow = 3, scales = "free")

png_b <- file.path(FIG_DIR, "fig7b_clr_scatter_shortlist.png")
pdf_b <- file.path(FIG_DIR, "fig7b_clr_scatter_shortlist.pdf")
ggsave(png_b, p_b, width = 7.2, height = 7.0, dpi = 300, bg = "white")
ggsave(pdf_b, p_b, width = 7.2, height = 7.0, bg = "white")
message("wrote ", png_b)
message("wrote ", pdf_b)

# ---- Sst_25 across independent pairs ---------------------------------------
sst <- dt[cell_type == "Sst_25"]
sst[, pair := factor(pair, levels = unique(dt$pair))]
p_c <- panel(sst,
             "Centred log-ratio of cell-type proportion (first study, ln)",
             "Centred log-ratio of cell-type proportion (second study, ln)") +
  facet_wrap(~ pair, nrow = 2, scales = "free")

png_c <- file.path(FIG_DIR, "fig7c_clr_scatter_sst25_pairs.png")
pdf_c <- file.path(FIG_DIR, "fig7c_clr_scatter_sst25_pairs.pdf")
ggsave(png_c, p_c, width = 7.2, height = 5.0, dpi = 300, bg = "white")
ggsave(pdf_c, p_c, width = 7.2, height = 5.0, bg = "white")
message("wrote ", png_c)
message("wrote ", pdf_c)
