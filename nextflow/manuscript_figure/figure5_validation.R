#!/usr/bin/env Rscript
# ============================================================
# Figure 5: Cross-modality validation of bulk cell-type GWAS signals
# ============================================================
# Layout:
#   A (top, full width) - Direction concordance bar chart across all 19 cell types
#   B (bottom-left)     - VIP TMEM106B cohort forest (15 cohorts)
#   C (bottom-middle)   - L5.6.IT.Car3 CACNA1C cohort forest (15 cohorts)
#   D (bottom-right)    - Microglia PRKN chr6 cohort forest (11 cohorts)
#
# All panels are pre-computed PNGs assembled with cowplot.
# Output: manuscript_figure/figure5_validation.{png,pdf,svg}
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(png)
  library(grid)
  library(scales)
})

# ── paths ─────────────────────────────────────────────────────────────────────
ROOT       <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
SN_DIR     <- file.path(ROOT, "results/sn_bulk_meta_similarity_design_matrix/top_hits")
FOREST_DIR <- file.path(ROOT, "results/meta_analysis_15cohorts/plots/forest")
OUT_DIR    <- file.path(ROOT, "manuscript_figure")

PNG_A <- file.path(SN_DIR,    "bulk_top_hits_direction_bar.png")
PNG_B <- file.path(FOREST_DIR, "VIP_TMEM106B_chr7_12284378_forest.png")
PNG_C <- file.path(FOREST_DIR, "L5.6.IT.Car3_CACNA1C_chr12_2324042_forest.png")
PNG_D <- file.path(FOREST_DIR, "Microglia_PRKN_chr6_164862615_forest.png")

LABEL_SIZE <- 16L

# ── helper: load PNG into ggdraw panel with bold label ────────────────────────
# Preserves natural aspect ratio (no forced fill); image is letterboxed in the
# cell rather than stretched.
load_png_panel <- function(path, label) {
  img  <- readPNG(path)
  grob <- rasterGrob(img, interpolate = TRUE)   # natural aspect ratio
  ggdraw() +
    draw_grob(grob) +
    draw_label(label, x = 0.01, y = 0.995, hjust = 0, vjust = 1,
               fontface = "bold", size = LABEL_SIZE) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# ── shared ancestry legend (replaces per-panel legends) ───────────────────────
make_shared_legend <- function() {
  clr_eur  <- "#2166AC"
  clr_afr  <- "#D95F02"
  clr_meta <- "#B2182B"
  leg_df <- data.frame(
    x   = 1:3,
    grp = factor(c("European", "Mixed ancestries", "Meta-analysis"),
                 levels = c("European", "Mixed ancestries", "Meta-analysis"))
  )
  # Return full ggplot (not extracted grob) so plot_grid renders it reliably
  ggplot(leg_df, aes(x = x, y = x, colour = grp, shape = grp)) +
    geom_point(size = 4, alpha = 0) +    # invisible points; only legend matters
    scale_colour_manual(
      name   = "Ancestry (panels B-D)",
      values = c("European"         = clr_eur,
                 "Mixed ancestries" = clr_afr,
                 "Meta-analysis"    = clr_meta)
    ) +
    scale_shape_manual(
      name   = "Ancestry (panels B-D)",
      values = c("European" = 16, "Mixed ancestries" = 16, "Meta-analysis" = 18)
    ) +
    guides(colour = guide_legend(nrow = 1,
                                 override.aes = list(size = 4, alpha = 1)),
           shape = "none") +
    theme_classic(base_size = 12) +
    theme(
      axis.line         = element_blank(),
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      axis.title        = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.title      = element_text(size = 11, face = "bold"),
      legend.text       = element_text(size = 10),
      legend.key.size   = unit(1.0, "lines"),
      legend.spacing.x  = unit(0.6, "lines")
    )
}

# ── load panels ───────────────────────────────────────────────────────────────
cat("Loading panels...\n")
pA   <- load_png_panel(PNG_A, "A")
pB   <- load_png_panel(PNG_B, "B")
pC   <- load_png_panel(PNG_C, "C")
pD   <- load_png_panel(PNG_D, "D")
pLeg <- make_shared_legend()   # full ggplot; theme_void shows only the legend
cat("  Panels A, B, C, D loaded.\n")

# ── assemble ──────────────────────────────────────────────────────────────────
cat("Assembling Figure 5...\n")

# Bottom row: B, C, D side by side; shared ancestry legend centred below
forest_row <- plot_grid(pB, pC, pD, ncol = 3, nrow = 1)
bottom_row <- plot_grid(
  forest_row,
  pLeg,              # gtable grob — plot_grid handles it directly
  ncol = 1, nrow = 2,
  rel_heights = c(1, 0.22)   # legend row ~22% of forest height (~1.8")
)

# Full figure: A on top (full width), bottom row below.
#   Panel A (2600×1100 px, 13"×5.5"):  18" wide → natural h = 18 × 0.423 = 7.62"
#   bottom_row (forests 8" + legend 22%): total ~9.76"
# rel_heights match natural heights exactly → zero letterbox gap.
fig5 <- plot_grid(
  pA,
  bottom_row,
  ncol = 1,
  nrow = 2,
  rel_heights = c(7.62, 9.76)
)

# ── save ──────────────────────────────────────────────────────────────────────
FIG_W <- 18
FIG_H <- 17   # 7.62 + 9.76 = 17.38"

out_png <- file.path(OUT_DIR, "figure5_validation.png")
out_pdf <- file.path(OUT_DIR, "figure5_validation.pdf")
out_svg <- file.path(OUT_DIR, "figure5_validation.svg")

png(out_png, width = FIG_W, height = FIG_H, units = "in",
    res = 200, type = "cairo")
print(fig5)
dev.off()

cairo_pdf(out_pdf, width = FIG_W, height = FIG_H)
print(fig5)
dev.off()

svg(out_svg, width = FIG_W, height = FIG_H, onefile = TRUE)
print(fig5)
dev.off()

cat("\nFigure 5 saved:\n  ", out_png, "\n  ", out_pdf, "\n  ", out_svg, "\n")
