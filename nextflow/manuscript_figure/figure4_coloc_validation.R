#!/usr/bin/env Rscript
# ============================================================
# Figure 4: Disease colocalization of cell-type proportion GWAS signals
# ============================================================
# Layout (2×2 grid):
#   A (top-left)    - Colocalization heatmap (19 cell types × 6 diseases)
#   B (top-right)   - VIP × MDD  TMEM106B  regional plot  (PP.H4 = 1.000)
#   C (bottom-left) - L5.6.IT.Car3 × BD  CACNA1C regional (PP.H4 = 0.877)
#   D (bottom-right)- Microglia × BD  chr6 regional       (PP.H4 = 0.616)
#
# All regional PNGs pre-computed by scripts/coloc/05_regional_coloc_plots.R.
# Output: manuscript_figure/figure4_coloc_validation.{png,pdf}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(png)
  library(grid)
})

# ── paths ────────────────────────────────────────────────────────────────────
ROOT       <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
COLOC_FILE <- file.path(ROOT, "results/coloc/full/coloc_all_results.tsv")
REG_DIR    <- file.path(ROOT, "results/coloc/full/plots/regional")
OUT_DIR    <- file.path(ROOT, "manuscript_figure")

REG_VIP    <- file.path(REG_DIR, "VIP_chr7_12284430_MDD_MDD2025_regional.png")
REG_CAR3   <- file.path(REG_DIR, "L5.6.IT.Car3_chr12_2324042_BD_bip2024_regional.png")
REG_MICRO  <- file.path(REG_DIR, "Microglia_chr6_164862615_BD_bip2024_regional.png")

LABEL_SIZE <- 16L   # uniform bold panel-label pt size for A, B, C, D

# ── helpers ──────────────────────────────────────────────────────────────────
load_png_panel <- function(path, label) {
  img  <- readPNG(path)
  grob <- rasterGrob(img, interpolate = TRUE)
  ggdraw() +
    draw_grob(grob) +
    draw_label(label, x = 0.01, y = 0.995, hjust = 0, vjust = 1,
               fontface = "bold", size = LABEL_SIZE)
}

# ============================================================
# SECTION 1: Panel A — Colocalization heatmap
# ============================================================
cat("Building Panel A: colocalization heatmap...\n")

coloc <- fread(COLOC_FILE)

# Clean disease labels
disease_map <- c(
  MDD_MDD2025        = "MDD",
  BD_bip2024         = "BD",
  SCZ_Trubetskoy2022 = "SCZ",
  AD_Bellenguez2022  = "AD",
  LBD_Chia2021       = "LBD",
  PD_Nalls2019       = "PD"
)
coloc[, disease_clean := disease_map[disease]]
coloc <- coloc[!is.na(disease_clean)]   # drop unlabelled diseases

# Cell-type class for Y-axis ordering
ct_class_map <- c(
  IT = "Excitatory", L4.IT = "Excitatory", L5.ET = "Excitatory",
  `L5.6.IT.Car3` = "Excitatory", `L5.6.NP` = "Excitatory",
  L6.CT = "Excitatory", L6b = "Excitatory",
  VIP = "Inhibitory", SST = "Inhibitory", PVALB = "Inhibitory",
  LAMP5 = "Inhibitory", PAX6 = "Inhibitory",
  Astrocyte = "Non-neuronal", Oligodendrocyte = "Non-neuronal",
  OPC = "Non-neuronal", Microglia = "Non-neuronal",
  Endothelial = "Non-neuronal", Pericyte = "Non-neuronal",
  VLMC = "Non-neuronal"
)
coloc[, ct_class := ct_class_map[cell_type]]
coloc[is.na(ct_class), ct_class := "Other"]

# Max PP.H4 per cell_type × disease (across all loci for that disease)
heat <- coloc[, .(max_h4 = max(PP.H4, na.rm = TRUE)),
              by = .(cell_type, disease_clean, ct_class)]

# Y-axis: within each class, order by total PP.H4 signal (descending)
ct_rank <- heat[, .(score = sum(max_h4)), by = .(cell_type, ct_class)][
  order(match(ct_class, c("Non-neuronal", "Inhibitory", "Excitatory")), -score)]
ct_levels <- ct_rank$cell_type   # non-neuronal first → shows at top of heatmap

heat[, cell_type    := factor(cell_type,     levels = rev(ct_levels))]
heat[, disease_clean := factor(disease_clean,
                                levels = c("MDD", "BD", "SCZ", "AD", "LBD", "PD"))]

# Focus cell types highlighted in red (all 3 including Microglia)
FOCUS_CTS <- c("VIP", "L5.6.IT.Car3", "Microglia")
FOCUS_COL <- "#c0392b"

# Build per-label colour/face vectors aligned to factor levels (bottom→top on Y)
y_colours <- ifelse(levels(heat$cell_type) %in% FOCUS_CTS, FOCUS_COL, "grey20")
y_faces   <- ifelse(levels(heat$cell_type) %in% FOCUS_CTS, "bold",    "plain")

# Number of cell types per class (for separator positions)
n_exc  <- sum(ct_rank$ct_class == "Excitatory")
n_inh  <- sum(ct_rank$ct_class == "Inhibitory")
# Separators: between Excitatory/Inhibitory and Inhibitory/Non-neuronal
sep_pos <- c(n_exc + 0.5, n_exc + n_inh + 0.5)

pA_gg <- ggplot(heat, aes(x = disease_clean, y = cell_type, fill = max_h4)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(data = heat[max_h4 > 0.5],
            aes(label = "\u2605"), colour = "white", size = 3.5, fontface = "bold") +
  scale_fill_gradientn(
    name    = "Max PP.H4",
    colours = c("grey97", "#FEE8C8", "#FDBB84", "#E34A33", "#8B0000"),
    values  = scales::rescale(c(0, 0.05, 0.3, 0.6, 1.0)),
    limits  = c(0, 1),
    breaks  = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  geom_hline(yintercept = sep_pos, linetype = "dashed",
             colour = "grey70", linewidth = 0.4) +
  annotate("text",
           x     = rep(6.65, 3),
           y     = c(n_exc / 2,
                     n_exc + n_inh / 2,
                     n_exc + n_inh + sum(ct_rank$ct_class == "Non-neuronal") / 2),
           label = c("Excitatory", "Inhibitory", "Non-neuronal"),
           hjust = 0, size = 2.8, colour = "grey45", fontface = "italic") +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +          # title handled via draw_label below
  theme_classic(base_size = 13) +
  theme(
    axis.text.x       = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.y       = element_text(size = 11, colour = y_colours, face = y_faces),
    axis.ticks        = element_blank(),
    axis.line         = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(1.2, "cm"),
    legend.title      = element_text(size = 10, face = "bold"),
    legend.text       = element_text(size = 9),
    panel.grid        = element_blank(),
    plot.margin       = margin(18, 55, 5, 5)  # top margin for draw_label
  )

# Wrap with uniform panel label "A" (same size/style as B/C/D)
pA <- ggdraw() +
  draw_plot(pA_gg) +
  draw_label("A", x = 0.01, y = 0.995, hjust = 0, vjust = 1,
             fontface = "bold", size = LABEL_SIZE)

cat("  Panel A done.\n")

# ============================================================
# SECTION 2: Panels B, C, D — Regional plots (load pre-computed PNGs)
# ============================================================
cat("Loading regional plot PNGs...\n")

pB <- load_png_panel(REG_VIP,   "B")
pC <- load_png_panel(REG_CAR3,  "C")
pD <- load_png_panel(REG_MICRO, "D")

cat("  Panels B, C, D loaded.\n")

# ============================================================
# SECTION 3: Assemble Figure 4
# ============================================================
cat("Assembling Figure 4...\n")

# 2 × 2 grid:  A (top-left)  |  B (top-right)
#              C (bot-left)  |  D (bot-right)
fig4 <- plot_grid(
  pA, pB,
  pC, pD,
  ncol = 2,
  nrow = 2,
  rel_heights = c(1.1, 1)
)

# ── save ─────────────────────────────────────────────────────────────────────
FIG_W <- 18   # inches: 9" per column
FIG_H <- 20   # inches: 10" per row

out_png <- file.path(OUT_DIR, "figure4_coloc_validation.png")
out_pdf <- file.path(OUT_DIR, "figure4_coloc_validation.pdf")
out_svg <- file.path(OUT_DIR, "figure4_coloc_validation.svg")

png(out_png, width = FIG_W, height = FIG_H, units = "in",
    res = 200, type = "cairo")
print(fig4)
dev.off()

cairo_pdf(out_pdf, width = FIG_W, height = FIG_H)
print(fig4)
dev.off()

svg(out_svg, width = FIG_W, height = FIG_H, onefile = TRUE)
print(fig4)
dev.off()

cat("\nFigure 4 saved:\n  ", out_png, "\n  ", out_pdf, "\n  ", out_svg, "\n")
