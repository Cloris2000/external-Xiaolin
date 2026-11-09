# Figure 8 supplement: per-supertype layer dispersion vs concordance, and the
# depth negative control. See fig8_dispersion_data.py.
#   (a) Xenium depth spread (sigma) vs median symmetric ICC across the 5 pairs
#   (b) Xenium white-matter share (p_wm) vs the same ICC
#   (c) sigma vs the supertype's own depth sensitivity: pooled slope of its
#       |residual CLR difference| on Wasserstein depth mismatch (pair as fixed
#       effect). Prediction: dispersed types flat. Excitatory types define x,
#       so they are hollow, reference only.
#   (d) p_wm vs depth sensitivity
# Class colours as Figure 1 (Allen-style). Spearman rho over inhibitory +
# glial supertypes (the non-circular set) in each panel.
# Inputs: data/fig8_supertype_dispersion.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig8_dispersion.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(ggrepel)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure8"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")

CLS_COL <- c("Excitatory" = "#D1491C", "Inhibitory" = "#2F4FA8", "Non-neuronal" = "#2E8B3C")
N_LABEL <- 4

d <- fread(file.path(HERE, "data/fig8_supertype_dispersion.tsv"))
d <- d[!is.na(sigma)]
d[, cls := factor(fifelse(cls == "Glia", "Non-neuronal", cls), levels = names(CLS_COL))]
d[, exc := cls == "Excitatory"]
cat(sprintf("%d supertypes with a reference depth; %d with ICC; %d with a depth slope\n",
            nrow(d), sum(!is.na(d$median_icc_sym)), sum(!is.na(d$depth_slope))))

stat_lab <- function(x, y, sub) {
  s <- spear(x, y)
  sprintf("%s: rho = %+.2f, p = %s, n = %d", sub, s$estimate, fmt_p(s$p.value), length(x))
}
extremes <- function(dt, y) {
  dt <- dt[!is.na(get(y))]
  dt[order(get(y))][c(1:N_LABEL, (.N - N_LABEL + 1):.N)]
}

panel <- function(xvar, yvar, xlab, ylab, hollow_exc) {
  dd <- d[!is.na(get(yvar))]
  test <- dd[exc == FALSE]
  sub <- stat_lab(test[[xvar]], test[[yvar]], "inhibitory + glial")
  if (!hollow_exc) {
    e <- dd[exc == TRUE]
    sub <- paste0(sub, "\n", stat_lab(e[[xvar]], e[[yvar]], "excitatory"))
  }
  p <- ggplot(dd, aes(.data[[xvar]], .data[[yvar]])) +
    geom_smooth(data = test, method = "lm", formula = y ~ x, colour = "grey30", fill = "grey85",
                linewidth = 0.4, linetype = "dashed") +
    geom_point(aes(colour = cls, fill = cls, shape = exc), size = 1.4, stroke = 0.4) +
    geom_text_repel(data = extremes(test, yvar), aes(label = supertype), size = fontsize / .pt * 0.7,
                    min.segment.length = 0, segment.size = 0.25, max.overlaps = Inf, seed = 1,
                    colour = "grey20") +
    scale_colour_manual(values = CLS_COL, name = NULL) +
    scale_fill_manual(values = CLS_COL, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = if (hollow_exc) 1 else 21), guide = "none") +
    labs(x = xlab, y = ylab, subtitle = sub) +
    theme_cns +
    theme(plot.subtitle = element_text(size = fontsize - 1, colour = "grey30"))
  p
}

ICC_LAB <- "Cross-dataset reliability of the supertype\n(median symmetric ICC, 5 pairs)"
SLOPE_LAB <- "Depth sensitivity of the supertype\n(slope of |Δ residual CLR| on depth mismatch)"
p_a <- panel("sigma", "median_icc_sym", "Layer dispersion (Xenium depth spread, sigma)", ICC_LAB, FALSE)
p_b <- panel("p_wm", "median_icc_sym", "White-matter share of the supertype (Xenium p_wm)", ICC_LAB, FALSE)
p_c <- panel("sigma", "depth_slope", "Layer dispersion (Xenium depth spread, sigma)", SLOPE_LAB, TRUE)
p_d <- panel("p_wm", "depth_slope", "White-matter share of the supertype (Xenium p_wm)", SLOPE_LAB, TRUE)

fig <- (p_a + p_b) / (p_c + p_d) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  caption = "Hollow points (c, d): excitatory supertypes define the depth axis, shown for reference only") &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = fontsize, face = "bold", colour = "black"),
        plot.caption = element_text(size = fontsize - 1, colour = "grey30", hjust = 0))

save_fig(fig, file.path(FIG_DIR, "figure8_dispersion_vs_concordance"), height_in = 7.2, width_mm = 183)
