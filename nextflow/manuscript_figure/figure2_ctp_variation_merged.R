# =============================================================================
# Figure 2 (revised). Bulk deconvolution captures reproducible, heritable
# inter-individual variation in cell-type proportions (CTPs).
#
# Story beat 1 of the manuscript: deconvolution is accurate (A), CTPs vary
# across individuals (B) and are consistent across cohorts (C), and that
# variation is partly genetic (D) — motivating the GWAS in Figure 3.
#
# Panels:
#   A — Bulk-to-snRNA estimation accuracy per cell type (benchmark)
#       (mean Pearson r across ROSMAP/HBCC/MSBB validation cohorts ± SD)
#   B — Estimated cell-type composition across cohorts (dot-range:
#       donor IQR thick bar, cohort-median range thin whisker)
#   C — Cohort-level differences in cell-type estimates (z-score heatmap)
#   D — SNP heritability (LDSC h2 ± SE) per cell type from the
#       15-cohort meta-analysis
#
# Inputs (all pre-computed):
#   manuscript_figure/figure2_celltype_accuracy.tsv
#   manuscript_figure/combined_bulk_snrna_paired.tsv
#   manuscript_figure/bulk_celltype_proportions_long.tsv
#   results/meta_analysis_15cohorts/ldsc/summary/ldsc_h2_summary.tsv
#     (falls back to the 13-cohort summary with a warning if not yet present)
#
# Output:
#   figure2_ctp_variation_merged.png / .pdf / .svg
# =============================================================================

NF_DIR  <- "/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow"
OUT_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"

ACC_FILE    <- file.path(OUT_DIR, "figure2_celltype_accuracy.tsv")
PAIRED_FILE <- file.path(OUT_DIR, "combined_bulk_snrna_paired.tsv")
PROP_FILE   <- file.path(OUT_DIR, "bulk_celltype_proportions_long.tsv")

H2_FILE_15  <- file.path(NF_DIR, "results/meta_analysis_15cohorts/ldsc/summary/ldsc_h2_summary.tsv")
H2_FILE_13  <- file.path(NF_DIR, "results/meta_analysis_13cohorts/ldsc/summary/ldsc_h2_summary.tsv")

# Benchmark panel drops L4 IT (kept consistent with the published benchmark)
EXCLUDE_CELLTYPES_ACC <- c("l4_it")

FIG_WIDTH  <- 16
FIG_HEIGHT <- 11
FIG_DPI    <- 300

# =============================================================================
# Packages
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(scales)
  library(RColorBrewer)
})

make_clean_fn <- function(x) {
  sub("^_+|_+$", "",
      tolower(gsub("[^a-zA-Z0-9]+", "_",
                   gsub("([a-z])([A-Z])", "\\1_\\2", x))))
}

# =============================================================================
# Shared cell-type metadata (order, labels, classes, colors)
# =============================================================================
CELLTYPE_ORDER <- c(
  "oligodendrocyte", "opc", "astrocyte", "microglia",
  "endothelial", "pericyte", "vlmc",
  "it", "l4_it", "l5_6_np", "l5_et", "l6_ct", "l5_6_it_car3", "l6b",
  "pvalb", "sst", "vip", "lamp5", "pax6"
)

CT_DISPLAY <- c(
  oligodendrocyte = "Oligodendrocyte", opc = "OPC", astrocyte = "Astrocyte",
  microglia = "Microglia", endothelial = "Endothelial", pericyte = "Pericyte",
  vlmc = "VLMC",
  it = "IT", l4_it = "L4 IT", l5_6_np = "L5/6 NP", l5_et = "L5 ET",
  l6_ct = "L6 CT", l5_6_it_car3 = "L5/6 IT Car3", l6b = "L6b",
  pvalb = "PVALB", sst = "SST", vip = "VIP", lamp5 = "LAMP5", pax6 = "PAX6"
)

CT_CLASS <- c(
  oligodendrocyte = "Non-neuronal", opc = "Non-neuronal", astrocyte = "Non-neuronal",
  microglia = "Non-neuronal", endothelial = "Non-neuronal", pericyte = "Non-neuronal",
  vlmc = "Non-neuronal",
  it = "Excitatory", l4_it = "Excitatory", l5_6_np = "Excitatory",
  l5_et = "Excitatory", l6_ct = "Excitatory", l5_6_it_car3 = "Excitatory",
  l6b = "Excitatory",
  pvalb = "Inhibitory", sst = "Inhibitory", vip = "Inhibitory",
  lamp5 = "Inhibitory", pax6 = "Inhibitory"
)
CLASS_LEVELS <- c("Excitatory", "Inhibitory", "Non-neuronal")
CLASS_COLORS <- c("Excitatory"   = "#D55E00",
                  "Inhibitory"   = "#009E73",
                  "Non-neuronal" = "#7570B3")

COHORT_COLORS <- c(
  ROSMAP  = "#0072B2",   # blue
  HBCC    = "#E69F00",   # orange
  MSBB    = "#CC79A7",   # mauve/pink
  Mathys  = "#56B4E9",   # sky blue
  Ruzicka = "#009E73"    # bluish green
)

pal_g <- colorRampPalette(brewer.pal(9, "Blues")[4:9])(7)
pal_e <- colorRampPalette(brewer.pal(9, "Oranges")[4:9])(7)
pal_i <- colorRampPalette(brewer.pal(9, "Greens")[4:9])(5)
ct_colors <- setNames(c(pal_g, pal_e, pal_i), CELLTYPE_ORDER)

base_th <- theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(size = 12, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text       = element_text(size = 9),
    axis.title      = element_text(size = 10.5),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10, face = "bold")
  )

# =============================================================================
# Panel A — Estimation accuracy per cell type (benchmark)
# =============================================================================
cat("--- Panel A: bulk-to-snRNA accuracy ---\n")

acc_mean_df <- read_tsv(ACC_FILE, show_col_types = FALSE) %>%
  filter(!ct_clean %in% EXCLUDE_CELLTYPES_ACC)

paired_val <- read_tsv(PAIRED_FILE, show_col_types = FALSE) %>%
  mutate(ct_clean = make_clean_fn(cell_type)) %>%
  filter(!ct_clean %in% EXCLUDE_CELLTYPES_ACC,
         snrna_proportion < 0.99)

acc_df <- paired_val %>%
  group_by(validation_cohort, ct_clean) %>%
  summarise(
    pearson_r = suppressWarnings(
      cor(bulk_proportion, snrna_proportion,
          method = "pearson", use = "pairwise.complete.obs")),
    .groups = "drop"
  ) %>%
  filter(!is.na(pearson_r))

ord_tbl <- acc_mean_df %>%
  mutate(class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
  arrange(class_f, desc(mean_r))
ct_order_a <- ord_tbl$ct_clean

mean_a <- acc_mean_df %>%
  mutate(ct_ord  = factor(ct_clean, levels = ct_order_a),
         class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS))

coh_a <- acc_df %>%
  mutate(ct_ord  = factor(ct_clean, levels = ct_order_a),
         class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
  filter(!is.na(ct_ord))

panel_A <- ggplot(mean_a, aes(x = ct_ord, y = mean_r)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey55") +
  geom_errorbar(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r),
                width = 0.25, linewidth = 0.45, color = "grey60") +
  geom_point(data = coh_a,
             aes(x = ct_ord, y = pearson_r, color = validation_cohort),
             size = 1.7, alpha = 0.85, shape = 16, inherit.aes = FALSE) +
  geom_point(size = 3.0, shape = 18, color = "#08519c") +
  facet_grid(cols = vars(class_f), scales = "free_x", space = "free_x") +
  scale_color_manual(values = COHORT_COLORS, name = "Cohort",
                     breaks = names(COHORT_COLORS),
                     guide = guide_legend(override.aes = list(size = 3))) +
  scale_x_discrete(labels = function(x) {
    lbl <- CT_DISPLAY[x]; ifelse(is.na(lbl), x, lbl)
  }) +
  labs(x = NULL, y = "Pearson r (bulk vs snRNA-seq)") +
  base_th +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        axis.text.x        = element_text(size = 9, angle = 45, hjust = 1),
        strip.text         = element_text(size = 10.5, face = "bold"),
        strip.background   = element_rect(fill = "grey94", color = NA),
        panel.spacing.x    = unit(0.4, "lines"),
        legend.position    = "bottom",
        legend.margin      = margin(t = -4),
        panel.border       = element_rect(color = "grey80", fill = NA,
                                          linewidth = 0.4))

# =============================================================================
# Panels B & C — composition/variability + cohort consistency
# =============================================================================
cat("--- Panels B & C: composition + cohort heatmap ---\n")

prop_long <- read_tsv(PROP_FILE, show_col_types = FALSE)

present_types <- unique(prop_long$cell_type)
ordered_types <- c(CELLTYPE_ORDER[CELLTYPE_ORDER %in% present_types],
                   setdiff(present_types, CELLTYPE_ORDER))
prop_long <- prop_long %>%
  mutate(cell_type = factor(cell_type, levels = ordered_types))

COHORT_LEVELS <- unique(prop_long$cohort)

cohort_medians <- prop_long %>%
  group_by(cohort, cell_type) %>%
  summarise(cohort_median = median(proportion_percent, na.rm = TRUE),
            .groups = "drop")

celltype_donor_summary <- prop_long %>%
  group_by(cell_type) %>%
  summarise(
    overall_median = median(proportion_percent, na.rm = TRUE),
    q25_donor      = quantile(proportion_percent, 0.25, na.rm = TRUE),
    q75_donor      = quantile(proportion_percent, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cohort_range_summary <- cohort_medians %>%
  group_by(cell_type) %>%
  summarise(rng_min = min(cohort_median, na.rm = TRUE),
            rng_max = max(cohort_median, na.rm = TRUE),
            .groups = "drop")

panel_B_data <- left_join(celltype_donor_summary, cohort_range_summary,
                          by = "cell_type")

panel_B <- ggplot(panel_B_data, aes(y = cell_type, color = cell_type)) +
  geom_point(data = cohort_medians, aes(x = cohort_median),
             shape = 16, size = 1.0, alpha = 0.30) +
  geom_linerange(aes(xmin = rng_min, xmax = rng_max),
                 linewidth = 0.5, alpha = 0.55) +
  geom_linerange(aes(xmin = q25_donor, xmax = q75_donor),
                 linewidth = 2.0, alpha = 0.80) +
  geom_point(aes(x = overall_median, fill = cell_type),
             shape = 21, size = 2.6, color = "white", stroke = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.35) +
  scale_color_manual(values = ct_colors, guide = "none") +
  scale_fill_manual( values = ct_colors, guide = "none") +
  scale_y_discrete(limits = rev(ordered_types),
                   labels = function(x) {
                     lbl <- CT_DISPLAY[x]; ifelse(is.na(lbl), x, lbl)
                   }) +
  labs(x = "MGP estimated proportion score", y = NULL) +
  base_th +
  theme(panel.grid.major.y = element_blank(),
        plot.margin = margin(l = 14, t = 5.5, r = 5.5, b = 5.5))

# --- Panel C: z-score heatmap of cohort means ---
cohort_z <- prop_long %>%
  group_by(cohort, cell_type) %>%
  summarise(mean_est = mean(proportion_percent, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(
    ct_mean = mean(mean_est, na.rm = TRUE),
    ct_sd   = sd(mean_est, na.rm = TRUE),
    z_score = if_else(ct_sd > 0, (mean_est - ct_mean) / ct_sd, 0)
  ) %>%
  ungroup() %>%
  mutate(z_clip = pmax(pmin(z_score, 2.5), -2.5))

ct_z_wide <- cohort_z %>%
  select(cohort, cell_type, z_clip) %>%
  pivot_wider(names_from = cell_type, values_from = z_clip, values_fill = 0)
cohort_mat <- as.matrix(select(ct_z_wide, -cohort))
rownames(cohort_mat) <- ct_z_wide$cohort
hc <- hclust(dist(cohort_mat), method = "ward.D2")
cohort_hc_order <- rownames(cohort_mat)[hc$order]
cohort_z$cohort <- factor(cohort_z$cohort, levels = cohort_hc_order)

panel_C <- ggplot(cohort_z, aes(x = cell_type, y = cohort, fill = z_clip)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c",
    midpoint = 0, limits = c(-2.5, 2.5),
    name = "Relative\ncohort z",
    labels = label_number(accuracy = 0.1)
  ) +
  scale_x_discrete(labels = function(x) {
    lbl <- CT_DISPLAY[x]; ifelse(is.na(lbl), x, lbl)
  }) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(x = NULL, y = NULL) +
  base_th +
  theme(
    axis.text.x       = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 7.5),
    axis.line         = element_blank(),
    panel.grid        = element_blank(),
    legend.key.height = unit(0.7, "cm"),
    legend.key.width  = unit(0.3, "cm")
  )

# =============================================================================
# Panel D — SNP heritability (LDSC h2) per cell type
# =============================================================================
cat("--- Panel D: LDSC h2 ---\n")

if (file.exists(H2_FILE_15)) {
  h2_file <- H2_FILE_15
  h2_label <- "15-cohort meta-analysis"
} else {
  h2_file <- H2_FILE_13
  h2_label <- "13-cohort meta-analysis"
  warning("15-cohort LDSC h2 summary not found yet; using the 13-cohort ",
          "summary as a placeholder. Re-run after SLURM job completes.")
}
cat("  Using:", h2_file, "(", h2_label, ")\n")

h2_df <- read_tsv(h2_file, show_col_types = FALSE) %>%
  mutate(
    ct_clean = make_clean_fn(cell_type),
    h2       = as.numeric(observed_h2),
    h2_se    = as.numeric(observed_h2_se),
    class_f  = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)
  ) %>%
  filter(!is.na(h2), !is.na(class_f)) %>%
  mutate(sig_pos = (h2 - h2_se) > 0)

# Order cell types within class by h2 (desc), matching Panel A's block layout
h2_order <- h2_df %>% arrange(class_f, desc(h2)) %>% pull(ct_clean)
h2_df <- h2_df %>% mutate(ct_ord = factor(ct_clean, levels = h2_order))

panel_D <- ggplot(h2_df, aes(x = h2, y = ct_ord, color = class_f)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = h2 - h2_se, xmax = h2 + h2_se),
                 height = 0.25, linewidth = 0.5, alpha = 0.8) +
  geom_point(aes(shape = sig_pos), size = 2.6, fill = "white", stroke = 0.6) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = expression(h^2 - SE > 0),
                                `FALSE` = "n.s."),
                     name = NULL) +
  scale_color_manual(values = CLASS_COLORS, name = "Cell class") +
  scale_y_discrete(limits = rev(h2_order),
                   labels = function(x) {
                     lbl <- CT_DISPLAY[x]; ifelse(is.na(lbl), x, lbl)
                   }) +
  labs(x = expression("SNP heritability" ~ (h^2 %+-% SE)), y = NULL) +
  base_th +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"))

# =============================================================================
# Assemble
# =============================================================================
cat("--- Combining panels ---\n")

bottom_row <- plot_grid(
  panel_B, panel_C, panel_D,
  ncol       = 3,
  rel_widths = c(0.85, 1.25, 1.0),
  labels     = c("B", "C", "D"),
  label_size = 15,
  label_fontface = "bold",
  align = "h", axis = "tb"
)

fig2 <- plot_grid(
  panel_A, bottom_row,
  nrow        = 2,
  rel_heights = c(0.85, 1.0),
  labels      = c("A", ""),
  label_size  = 15,
  label_fontface = "bold"
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

out_png <- file.path(OUT_DIR, "figure2_ctp_variation_merged.png")
out_pdf <- file.path(OUT_DIR, "figure2_ctp_variation_merged.pdf")
out_svg <- file.path(OUT_DIR, "figure2_ctp_variation_merged.svg")

ggsave(out_png, fig2, width = FIG_WIDTH, height = FIG_HEIGHT,
       dpi = FIG_DPI, bg = "white")
cat("Saved PNG:", out_png, "\n")
ggsave(out_pdf, fig2, width = FIG_WIDTH, height = FIG_HEIGHT, bg = "white")
cat("Saved PDF:", out_pdf, "\n")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(out_svg, fig2, width = FIG_WIDTH, height = FIG_HEIGHT,
         device = svglite::svglite)
} else {
  svg(out_svg, width = FIG_WIDTH, height = FIG_HEIGHT); print(fig2); dev.off()
}
cat("Saved SVG:", out_svg, "\n")
cat("h2 panel source:", h2_label, "\n")
cat("Done.\n")
