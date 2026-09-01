#!/usr/bin/env Rscript
# =============================================================================
# Figure 1 – Panel B and Panel C (revised, clean design)
# =============================================================================
# Panels B and C are designed to sit side-by-side beneath an existing Panel A
# workflow schematic. They are intentionally simple and data-driven.
#
# Panel B: Clean horizontal bar chart of the 15 bulk GWAS cohort datasets.
# Panel C: Left = snRNA-seq validation sample sizes; Right = 4-step flow.
#
# All sample sizes are read from verified pipeline outputs. No values invented.
# =============================================================================

# ─── 0. User switches ────────────────────────────────────────────────────────
INCLUDE_PANEL_LETTERS <- TRUE   # set FALSE to strip "B" / "C" before Inkscape

# ─── 1. Paths ─────────────────────────────────────────────────────────────────
ROOT   <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
MS_DIR <- file.path(ROOT, "manuscript_figure")
dir.create(MS_DIR, showWarnings = FALSE, recursive = TRUE)

# ─── 2. Packages ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(cowplot)
  library(scales)
  library(svglite)
})

# ─── 3. Shared style ──────────────────────────────────────────────────────────
BASE_PT   <- 8          # minimum font size (points) at final figure size
FONT_FAM  <- "Helvetica"
PANEL_W   <- 7.2        # inches — identical for B and C
PANEL_H   <- 4.6        # inches — identical for B and C

# Ancestry palette
COL_EUR   <- "#2166AC"   # EUR-majority  (blue)
COL_AFR   <- "#D95F02"   # AFR-enriched  (orange)
COL_MIX   <- "#1B9E77"   # Mixed/diverse (teal)

theme_fig <- function() {
  theme_minimal(base_size = BASE_PT, base_family = FONT_FAM) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.y  = element_blank(),
      panel.grid.major.x  = element_line(linewidth = 0.2, colour = "grey88"),
      axis.text           = element_text(size = BASE_PT,   colour = "grey15"),
      axis.title.x        = element_text(size = BASE_PT,   colour = "grey15"),
      axis.title.y        = element_blank(),
      axis.ticks.y        = element_blank(),
      axis.line.y         = element_blank(),
      plot.title          = element_text(size = BASE_PT + 1.5, face = "bold",
                                         hjust = 0, colour = "grey5",
                                         margin = margin(b = 1)),
      plot.subtitle       = element_text(size = BASE_PT - 0.5, hjust = 0,
                                         colour = "grey40",
                                         margin = margin(b = 4)),
      plot.margin         = margin(5, 4, 4, 4),
      legend.position     = "bottom",
      legend.justification = "left",
      legend.key.size     = unit(0.30, "cm"),
      legend.text         = element_text(size = BASE_PT - 0.5),
      legend.title        = element_text(size = BASE_PT - 0.5, face = "bold"),
      legend.margin       = margin(t = 2),
      legend.background   = element_rect(fill = NA, colour = NA)
    )
}

save_panel <- function(p, stem, w = PANEL_W, h = PANEL_H) {
  svglite(file.path(MS_DIR, paste0(stem, ".svg")), width = w, height = h)
  print(p); dev.off()
  cairo_pdf(file.path(MS_DIR, paste0(stem, ".pdf")), width = w, height = h)
  print(p); dev.off()
  png(file.path(MS_DIR, paste0(stem, ".png")),
      width = w * 600, height = h * 600, res = 600,
      type = "cairo", bg = "white")
  print(p); dev.off()
  message("  Saved: ", stem)
}

# =============================================================================
# PANEL B ── Bulk cohorts included in the CTP GWAS
# =============================================================================

# ─── B.1  Source data ─────────────────────────────────────────────────────────
# Analyzed N: regenie step2 raw_p N column (constant across all 19 traits).
# Brain region: tissue_filter in nextflow.config.combined.* files.
# Ancestry composition: ancestry keep-list counts + config comments (conservative).
#   EUR-majority: EUR keep-list N ≥ 85% of GWAS N, or only EUR keep-list present.
#   AFR-enriched: AFR keep-list > EUR keep-list.
#   Mixed/diverse: substantial AFR + AMR components without EUR dominance.

bulk <- data.table(
  # display_name, n_gwas, family, region, ancestry_cat
  # Ordered bottom-to-top (index 1 = bottom bar in horizontal chart)
  display_name  = c(
    "AMP-AD Mayo",    "AMP-AD Rush",
    "HBCC (Omni5M)",  "HBCC (h650)",   "HBCC (1M)",
    "BrainGVEX",
    "NABEC",
    "GTEx v10",
    "CMC PITT",  "CMC PENN",  "CMC MSSM",
    "MSBB",      "Mayo",
    "ROSMAP (array)", "ROSMAP (WGS)"
  ),
  n_gwas = c(
    227L, 120L,
    79L, 97L, 202L,
    394L,
    210L,
    283L,
    161L, 92L, 242L,
    233L, 243L,
    162L, 756L
  ),
  family = c(
    "AMP-AD Diverse", "AMP-AD Diverse",
    "HBCC",           "HBCC",          "HBCC",
    "BrainGVEX",
    "NABEC",
    "GTEx",
    "CMC", "CMC", "CMC",
    "AMP-AD",  "AMP-AD",
    "ROSMAP", "ROSMAP"
  ),
  region = c(
    "DLPFC", "DLPFC",
    "DLPFC", "DLPFC", "DLPFC",
    "DLPFC",
    "FCX",
    "FCX",
    "DLPFC", "DLPFC", "DLPFC",
    "STG",   "TCX",
    "DLPFC", "DLPFC"
  ),
  # Ancestry classification (conservative, based on keep-list files):
  #   HBCC h650: AFR keep-list = 53 of N=97 (~55%) → AFR-enriched
  #   HBCC 1M:   AFR=101, EUR=88 of N=202 → Mixed/diverse
  #   HBCC Omni5M: EUR=50 of N=79 (~63%) → EUR-majority
  #   AMP-AD Rush: AFR keep-list = 64 of N=120 (~53%) → AFR-enriched
  #   AMP-AD Mayo: AFR=52, AMR=178 of N=227, no EUR keep-list → Mixed/diverse
  #   All others: EUR keep-list ≥ 85% of GWAS N → EUR-majority
  ancestry = c(
    "Mixed/diverse",  "AFR-enriched",
    "EUR-majority",   "AFR-enriched",  "Mixed/diverse",
    "EUR-majority",
    "EUR-majority",
    "EUR-majority",
    "EUR-majority", "EUR-majority", "EUR-majority",
    "EUR-majority", "EUR-majority",
    "EUR-majority", "EUR-majority"
  )
)

total_n_bulk <- sum(bulk$n_gwas)

# Factor so bottom→top matches display_name order (index 1 = y=1 = bottom)
bulk[, display_name := factor(display_name, levels = display_name)]
bulk[, y_pos := as.integer(display_name)]

# Ancestry color map
anc_map   <- c("EUR-majority"  = COL_EUR,
               "AFR-enriched"  = COL_AFR,
               "Mixed/diverse" = COL_MIX)
anc_labels <- c("EUR-majority"  = "EUR-majority",
                "AFR-enriched"  = "AFR-enriched",
                "Mixed/diverse" = "Mixed/diverse")

# Family separator y-positions (gap between consecutive families)
fam_seq   <- bulk$family
breaks_y  <- which(fam_seq[-1] != fam_seq[-length(fam_seq)]) + 0.5

BAR_H <- 0.78   # bar height in y-axis units (thick, clearly visible)
x_max <- 870    # x-axis upper limit for bars (N labels + region placed beyond)

# ─── B.2  Build a single ggplot embedding the region annotation ───────────────
# Region annotation is placed INSIDE the bar chart at x positions beyond x_max,
# so y-alignment is guaranteed (no separate plot needed).
region_cols <- c(
  DLPFC = "#4A9C6D", FCX = "#3975B7", TCX = "#C46B3A", STG = "#8E63AE"
)

# x-positions for region tiles and N-labels
REG_X  <- x_max + 85      # center of region tile column
REG_W  <- 68               # tile width (x-axis units)
N_X    <- x_max + 13       # x for N labels
X_LIM  <- REG_X + REG_W / 2 + 12   # total x-axis limit

# Build region tile lookup for bulk
bulk[, reg_col := region_cols[region]]

pB <- ggplot(bulk) +
  # family gap lines
  geom_hline(yintercept = breaks_y, colour = "grey82",
             linewidth = 0.28, linetype = "solid") +
  # main bars (ancestry-colored)
  geom_col(aes(y = y_pos, x = n_gwas, fill = ancestry),
           width = BAR_H, colour = NA) +
  # N label just right of bar
  geom_text(aes(y = y_pos, x = N_X, label = format(n_gwas, big.mark = ",")),
            hjust = 0, size = BASE_PT / ggplot2::.pt, colour = "grey20") +
  # region tiles (hardcoded fill via fill = I(reg_col))
  geom_tile(aes(y = y_pos, x = REG_X, fill = I(reg_col)),
            width = REG_W, height = BAR_H, colour = "white", linewidth = 0.25) +
  geom_text(aes(y = y_pos, x = REG_X, label = region),
            size = (BASE_PT - 1) / ggplot2::.pt,
            hjust = 0.5, vjust = 0.5,
            colour = "white", fontface = "bold") +
  # "Region" column header
  annotate("text", x = REG_X, y = max(bulk$y_pos) + 0.75,
           label = "Region", size = (BASE_PT - 1) / ggplot2::.pt,
           hjust = 0.5, vjust = 0, colour = "grey30", fontface = "bold") +
  # ancestry fill scale (bars only)
  scale_fill_manual(values = anc_map, name = "Cohort composition",
                    labels = anc_labels,
                    guide  = guide_legend(nrow = 1,
                                          override.aes = list(size = 3.5))) +
  scale_y_continuous(
    breaks = bulk$y_pos,
    labels = levels(bulk$display_name),
    expand = expansion(add = c(0.5, 0.85))   # extra top padding for "Region" header
  ) +
  scale_x_continuous(
    limits = c(0, X_LIM),
    breaks = c(0, 200, 400, 600, 800),
    labels = label_comma(),
    expand = c(0, 0)
  ) +
  labs(
    x        = "Analyzed donors",
    y        = NULL,
    title    = if (INCLUDE_PANEL_LETTERS)
                 "B   Bulk cohorts included in the CTP GWAS"
               else "Bulk cohorts included in the CTP GWAS",
    subtitle = sprintf("15 cohort datasets  |  Total analyzed N = %s",
                       format(total_n_bulk, big.mark = ","))
  ) +
  theme_fig() +
  theme(
    legend.position = "bottom",
    # hide x-axis ticks/text beyond x_max (region area)
    axis.line.x = element_blank()
  )

# =============================================================================
# PANEL C ── Independent snRNA-seq validation
# =============================================================================
# Cohorts: hodge5u de-duplicated 5-cohort meta (ANALYSIS_TAG in figure5_sn_concordance.R)
# CTP donors: from data_input/*/cell_proportions.csv row counts.
# GWAS donors: from results/*/regenie_step2 N columns.

# ─── C.1  Source data ─────────────────────────────────────────────────────────
sn <- data.table(
  display_name  = c(
    "ROSMAP Green",     "PsychAD HBCC",
    "PsychAD MSSM",     "ROSMAP Mathys\n(unique)",
    "Ruzicka MSSM"
  ),
  ctp_n   = c(420L, 259L, 203L, 128L, 53L),
  gwas_n  = c(366L, 258L, 187L, 121L, 53L),
  # nuclei available only for Mathys (occupancy_report.tsv):
  # total = sum(n_cells) = 1,420,318
  nuclei_label = c("", "", "", "1.4M nuclei", "")
)

# Bottom→top ordering
sn[, display_name := factor(display_name, levels = display_name)]
sn[, y_pos := as.integer(display_name)]

SN_BAR_H   <- 0.70
bar_pale   <- "#9DC3E6"   # pale blue  – CTP donors
bar_dark   <- "#2E75B6"   # dark blue  – GWAS donors
bar_pale_c <- "#5B9BD5"   # pale bar outline

x_max_sn   <- 460   # x-axis limit (labels placed via clip = "off")

# ─── C.2  Left: sample-size bar chart ─────────────────────────────────────────
pC_bars <- ggplot(sn) +
  # pale outlined bar = CTP donors
  geom_col(aes(x = ctp_n, y = y_pos),
           fill = bar_pale, colour = bar_pale_c,
           linewidth = 0.35, width = SN_BAR_H) +
  # dark filled bar = GWAS donors (overlaid)
  geom_col(aes(x = gwas_n, y = y_pos),
           fill = bar_dark, colour = NA,
           width = SN_BAR_H) +
  # end label: "CTP / GWAS donors" — placed just beyond bar end
  geom_text(aes(x = ctp_n + 8, y = y_pos,
                label = paste0(ctp_n, " / ", gwas_n)),
            hjust = 0, size = BASE_PT / ggplot2::.pt,
            colour = "grey20") +
  # nuclei annotation (small, italic) where available
  geom_text(data = sn[nuclei_label != ""],
            aes(x = 8, y = y_pos - 0.32, label = nuclei_label),
            hjust = 0, size = (BASE_PT - 1.5) / ggplot2::.pt,
            colour = "grey50", fontface = "italic") +
  scale_y_continuous(
    breaks = sn$y_pos,
    labels = levels(sn$display_name),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_x_continuous(
    limits = c(0, x_max_sn),
    breaks = c(0, 100, 200, 300, 400),
    labels = label_comma(),
    expand = c(0, 0)
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, x_max_sn), clip = "off") +
  labs(
    x     = "Number of donors",
    y     = NULL,
    title = if (INCLUDE_PANEL_LETTERS)
              "C   Independent snRNA-seq validation"
            else "Independent snRNA-seq validation"
  ) +
  theme_fig() +
  theme(
    legend.position = "none",
    plot.margin     = margin(5, 55, 4, 4)   # right margin to show clipped labels
  )

# ─── C.3  Inline legend ───────────────────────────────────────────────────────
leg_data <- data.table(
  label = c("snRNA-seq CTP available", "Included in GWAS"),
  fill  = c(bar_pale, bar_dark),
  col   = c(bar_pale_c, bar_dark),
  y     = c(2, 1)
)
pC_leg <- ggplot(leg_data, aes(x = 1, y = y)) +
  geom_tile(aes(fill = I(fill), colour = I(col)),
            width = 0.30, height = 0.55, linewidth = 0.4) +
  geom_text(aes(label = label), x = 1.22,
            hjust = 0, size = (BASE_PT - 0.5) / ggplot2::.pt,
            colour = "grey20", vjust = 0.5) +
  xlim(0.7, 5) + ylim(0.3, 2.7) +
  theme_void() + theme(plot.margin = margin(1, 0, 4, 4))

# Stack bars + legend
pC_left <- plot_grid(pC_bars, pC_leg,
                     ncol = 1, rel_heights = c(1, 0.18),
                     align = "v", axis = "l")

# ─── C.4  Right: 4-step horizontal flow ───────────────────────────────────────
# 4 nodes at x = 1, 3.2, 5.4, 7.6 on a 0–9 axis.
# Node labels in bold; sub-labels below steps 2 and 4 in lighter color.
# No boxes. Arrows connect nodes at mid-height.

NODE_X  <- c(1, 3.2, 5.4, 7.6)
NODE_Y  <- 0.72   # y of main labels
SUB_Y   <- 0.32   # y of sub-labels
ARR_Y   <- NODE_Y + 0.02

steps <- data.frame(
  x      = NODE_X,
  label  = c("Annotated\nnuclei",
             "Direct\ndonor-level\nCTPs",
             "snRNA-seq\nGWAS",
             "Comparison\nwith bulk\nmeta-analysis"),
  sub    = c("",
             "19 harmonized\ncell types",
             "",
             "Effect-direction\n& genome-wide\nconcordance"),
  stringsAsFactors = FALSE
)
steps$has_sub <- nchar(steps$sub) > 0

# Arrow x-ranges (from just right of one label to just left of next)
ARROWX1 <- NODE_X[-4] + 0.65
ARROWX2 <- NODE_X[-1] - 0.65

pC_flow <- ggplot() +
  # bold step labels, shifted up slightly when sub-label present
  geom_text(data = steps,
            aes(x = x,
                y = ifelse(has_sub, NODE_Y + 0.06, NODE_Y),
                label = label),
            size = (BASE_PT + 0.2) / ggplot2::.pt,
            hjust = 0.5, vjust = 0.5,
            colour = "grey10", fontface = "bold",
            lineheight = 0.92) +
  # sub-labels in blue
  geom_text(data = steps[steps$has_sub, ],
            aes(x = x, y = SUB_Y, label = sub),
            size = (BASE_PT - 1.5) / ggplot2::.pt,
            hjust = 0.5, vjust = 0.5,
            colour = "#2E75B6", lineheight = 0.88) +
  # horizontal arrows
  annotate("segment",
           x    = ARROWX1, xend = ARROWX2,
           y    = ARR_Y,   yend = ARR_Y,
           arrow = arrow(length = unit(0.11, "cm"), type = "closed"),
           linewidth = 0.5, colour = "grey50") +
  # light dashed connector from main label to sub-label
  geom_segment(data = steps[steps$has_sub, ],
               aes(x = x, xend = x,
                   y = NODE_Y - 0.09, yend = SUB_Y + 0.14),
               linewidth = 0.25, colour = "grey70", linetype = "dashed") +
  xlim(0, 9) + ylim(0.12, 1.05) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(plot.margin = margin(16, 6, 8, 6))

# ─── C.5  Assemble Panel C ────────────────────────────────────────────────────
pC <- plot_grid(
  pC_left,
  pC_flow,
  nrow = 1, rel_widths = c(0.50, 0.50),
  align = "h", axis = "tb"
)

# =============================================================================
# SAVE INDIVIDUAL PANELS
# =============================================================================
message("Saving Panel B ...")
save_panel(pB, "figure1_panelB_revised")

message("Saving Panel C ...")
save_panel(pC, "figure1_panelC_revised")

# =============================================================================
# COMBINED PREVIEW
# =============================================================================
message("Saving combined preview ...")
pPreview <- plot_grid(pB, pC,
                      nrow = 1, rel_widths = c(1, 1))

png(file.path(MS_DIR, "figure1_panels_BC_revised_preview.png"),
    width = PANEL_W * 2 * 600, height = PANEL_H * 600,
    res = 600, type = "cairo", bg = "white")
print(pPreview); dev.off()
message("  Saved: figure1_panels_BC_revised_preview.png")
message("Done.")
