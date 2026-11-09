# Draft schematic for Fig 9a (user 2026-11-09): layout guide for the BioRender version.
# Left: how one CelMod model is trained, applied and scored (as in panels c-d).
# Right: the three training strategies compared in panel b, over the six test cohorts.
# Sourced by fig9.R (needs ggplot2 + STRAT_COL); also writes a standalone PDF/SVG.
make_schematic <- function(col = STRAT_COL, pt = 5.5) {
  fs <- pt / .pt
  light <- function(c, a = 0.25) grDevices::adjustcolor(c, alpha.f = a)
  box <- function(x0, x1, y0, y1, fill = "white", colour = "grey45", lw = 0.3)
    annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, colour = colour, linewidth = lw)
  txt <- function(x, y, l, ...) annotate("text", x = x, y = y, label = l, size = fs, lineheight = 0.85, ...)
  arr <- function(x, y, xe, ye, colour = "grey30", lw = 0.3, ends = "last")
    annotate("segment", x = x, y = y, xend = xe, yend = ye, colour = colour, linewidth = lw,
             arrow = arrow(length = unit(0.8, "mm"), type = "closed", ends = ends))
  BULK_F <- "grey93"; SN_F <- "#FBE3E0"
  AGE_F <- "#EEE8F4"; PSY_F <- "#E4F1EC"

  # ---- left: train -> apply -> score -------------------------------------------
  L <- list(
    txt(20, 58, "One CelMod model", fontface = "bold"),
    box(1, 19, 45, 53, BULK_F), txt(10, 49, "Bulk RNA-seq\n(training cohort)"),
    box(21, 39, 45, 53, SN_F), txt(30, 49, "snRNA-seq\ncomposition\n(same donors)"),
    arr(10, 45, 16, 39.5), arr(30, 45, 24, 39.5),
    box(8, 32, 33, 39.5, "white", "black", 0.45), txt(20, 36.2, "CelMod\n(genes \u2192 cell types)"),
    box(1, 12.5, 18, 26, BULK_F), txt(6.75, 22, "Bulk\nRNA-seq\n(test cohort)"),
    arr(20, 33, 20, 26.5),
    arr(12.5, 22, 14.5, 22),
    box(14.5, 25.5, 18, 26, "white", "black", 0.45), txt(20, 22, "Predicted\ncomposition"),
    box(27.5, 40, 18, 26, SN_F), txt(33.75, 22, "snRNA-seq\nground truth\n(test donors)"),
    arr(25.5, 22, 27.5, 22, ends = "both"),
    txt(20, 11, "Correlation per Core Supertype\n(c: median; d: per type)\ntraining donors excluded;\nsame dataset \u2192 5-fold CV",
        colour = "grey20"),
    arr(30.75, 18, 30.75, 15, colour = "grey50")
  )

  # ---- right: three strategies --------------------------------------------------
  coh <- c("ROSMAP", "MSBB", "Diverse", "CMC", "HBCC", "NABEC")
  xs <- seq(49.5, 95.5, length.out = 6); hw <- 4.1
  cfill <- c(rep(AGE_F, 3), rep(PSY_F, 3))
  targets <- function(y0, y1) c(lapply(1:6, function(i) box(xs[i] - hw, xs[i] + hw, y0, y1, cfill[i])),
                                lapply(1:6, function(i) txt(xs[i], (y0 + y1) / 2, coh[i])))
  label <- function(y, s, l) txt(44, y, l, hjust = 0, colour = col[[s]], fontface = "bold")
  R <- c(
    list(annotate("segment", x = 42, xend = 42, y = 2, yend = 59, colour = "grey80", linewidth = 0.3),
         annotate("segment", x = xs[1] - hw, xend = xs[3] + hw, y = 56.3, yend = 56.3, colour = "grey50", linewidth = 0.3),
         annotate("segment", x = xs[4] - hw, xend = xs[6] + hw, y = 56.3, yend = 56.3, colour = "grey50", linewidth = 0.3),
         txt(xs[2], 58.2, "Aging / AD cohorts"), txt(xs[5], 58.2, "Adult / psychiatric cohorts"),
         # S1
         label(53, "S1", "S1 paired: train on each cohort's own snRNA-seq")),
    unlist(lapply(1:6, function(i) list(
      box(xs[i] - 3.2, xs[i] + 3.2, 46.5, 50, light(col[["S1"]], 0.35), col[["S1"]]),
      txt(xs[i], 48.25, "own sn"),
      arr(xs[i], 46.5, xs[i], 43.3, col[["S1"]]))), recursive = FALSE),
    targets(39.5, 43.3),
    # S2
    list(label(35, "S2", "S2 global: one model for every cohort"),
         box(56, 89, 28.5, 32.5, light(col[["S2"]], 0.3), col[["S2"]]),
         txt(72.5, 30.5, "best ROSMAP multi-dataset consensus\n(chosen on the other cohorts)")),
    lapply(1:6, function(i) arr(72.5 + (xs[i] - 72.5) * 0.35, 28.5, xs[i], 25.3, col[["S2"]])),
    targets(21.5, 25.3),
    # S3
    list(label(17, "S3", "S3 class-matched: train within diagnosis/age class"),
         box(xs[1] - hw, xs[3] + hw, 10.5, 14.5, light(col[["S3"]], 0.3), col[["S3"]]),
         txt(xs[2], 12.5, "best aging/AD-trained model"),
         box(xs[4] - hw, xs[6] + hw, 10.5, 14.5, light(col[["S3"]], 0.3), col[["S3"]]),
         txt(xs[5], 12.5, "best psychiatric-trained model")),
    lapply(1:6, function(i) arr(xs[i], 10.5, xs[i], 7.3, col[["S3"]])),
    targets(3.5, 7.3)
  )
  ggplot() + L + R +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 60), expand = FALSE, clip = "off") +
    theme_void()
}

if (sys.nframe() == 0) {   # standalone: Rscript fig9_schematic.R
  HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure9"
  source(file.path(HERE, "..", "common", "theme_cns.R"))
  STRAT_COL <- c(S1 = "#7B8C9E", S2 = "#2F7E8F", S3 = "#E1A22E", prod = "black")
  p <- make_schematic()
  out <- file.path(HERE, "figures", "figure9a_schematic_draft")
  ggsave(paste0(out, ".pdf"), p, width = 110 / 25.4, height = 66 / 25.4, bg = "white",
         device = cairo_pdf)   # cairo keeps the arrow glyph
  ggsave(paste0(out, ".png"), p, width = 110 / 25.4, height = 66 / 25.4, dpi = 300, bg = "white")
  ggsave(file.path(HERE, "figures", "figure9a_schematic_draft.svg"), p, width = 110 / 25.4,
         height = 66 / 25.4, bg = "white", device = svglite::svglite)
}
