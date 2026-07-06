#!/usr/bin/env Rscript
# Visualise heterogeneity breakdown across cohort groups.
#
# Reads the three TSVs produced by heterogeneity_breakdown.py and generates:
#   1. Forest plots  — one PDF per top hit, cohorts colour-coded by each dimension.
#   2. Summary panel — I² distribution violin plots per dimension/group.
#   3. Heatmap       — cell type × cohort-group median I² at suggestive hits.
#
# Usage:
#   Rscript plot_heterogeneity_breakdown.R \
#     --per-cohort-effects   <dir>/per_cohort_effects.tsv \
#     --group-i2-summary     <dir>/group_i2_summary.tsv \
#     --overall-group-i2     <dir>/overall_group_i2.tsv \
#     --output-dir           <dir>/figures \
#     --forest-dir           <dir>/forest_plots \
#     [--max-forest-plots 50]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
args_raw <- commandArgs(trailingOnly = TRUE)

parse_args <- function(argv) {
  out <- list(
    per_cohort_effects = NULL,
    group_i2_summary   = NULL,
    overall_group_i2   = NULL,
    output_dir         = "figures",
    forest_dir         = "forest_plots",
    max_forest_plots   = 50L
  )
  i <- 1L
  while (i <= length(argv)) {
    flag <- argv[i]
    val  <- if (i + 1L <= length(argv)) argv[i + 1L] else NA_character_
    if      (flag == "--per-cohort-effects") { out$per_cohort_effects <- val; i <- i + 2L }
    else if (flag == "--group-i2-summary")   { out$group_i2_summary   <- val; i <- i + 2L }
    else if (flag == "--overall-group-i2")   { out$overall_group_i2   <- val; i <- i + 2L }
    else if (flag == "--output-dir")         { out$output_dir         <- val; i <- i + 2L }
    else if (flag == "--forest-dir")         { out$forest_dir         <- val; i <- i + 2L }
    else if (flag == "--max-forest-plots")   { out$max_forest_plots   <- as.integer(val); i <- i + 2L }
    else { i <- i + 1L }
  }
  out
}

args <- parse_args(args_raw)

if (is.null(args$per_cohort_effects) ||
    is.null(args$group_i2_summary)   ||
    is.null(args$overall_group_i2)) {
  stop("Required arguments: --per-cohort-effects, --group-i2-summary, --overall-group-i2")
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$forest_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
message("Loading data...")
eff   <- read.delim(args$per_cohort_effects, stringsAsFactors = FALSE)
grp   <- read.delim(args$group_i2_summary,   stringsAsFactors = FALSE)
ovrl  <- read.delim(args$overall_group_i2,   stringsAsFactors = FALSE)

# Convert types
eff$beta           <- as.numeric(eff$beta)
eff$se             <- as.numeric(eff$se)
eff$meta_HetISq    <- suppressWarnings(as.numeric(eff$meta_HetISq))
grp$within_group_I2 <- suppressWarnings(as.numeric(grp$within_group_I2))

# ---------------------------------------------------------------------------
# Colour palettes per dimension
# ---------------------------------------------------------------------------
diag_colours <- c(
  AD_neurological       = "#E41A1C",
  psychiatric_mixed     = "#FF7F00",
  neurologically_normal = "#4DAF4A",
  mixed_or_unclassified = "#984EA3",
  META                  = "#222222"
)

anc_colours <- c(
  EUR_homogeneous = "#2166AC",
  AFR_enriched    = "#D95F02",
  mixed_ancestry  = "#7B3294",
  META            = "#222222"
)

atype_colours <- c(
  DLPFC_array_imputed = "#1B9E77",
  DLPFC_WGS           = "#D95F02",
  non_DLPFC_brain     = "#7570B3",
  control_brain       = "#E7298A",
  control_brain_WGS   = "#E7298A",
  other               = "#999999",
  META                = "#222222"
)

dim_palette <- list(
  diagnosis_group = diag_colours,
  ancestry_class  = anc_colours,
  analysis_type   = atype_colours
)

dim_labels <- c(
  diagnosis_group = "Diagnosis group",
  ancestry_class  = "Ancestry class",
  analysis_type   = "Analysis type"
)

# ---------------------------------------------------------------------------
# Helper: safe p-value label
# ---------------------------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3)
}

# ---------------------------------------------------------------------------
# 1. Forest plots — one PDF per top hit, faceted by dimension
# ---------------------------------------------------------------------------
message("Generating forest plots...")

# Identify unique (cell_type, marker) pairs, prioritise GW hits
hit_keys <- eff %>%
  distinct(cell_type, marker, is_gw, meta_HetISq) %>%
  arrange(desc(is_gw), desc(meta_HetISq))

n_plots <- min(nrow(hit_keys), args$max_forest_plots)
message("  Plotting ", n_plots, " top hits (max = ", args$max_forest_plots, ")")

# Precompute per-cohort per-hit p-values from beta/se
eff <- eff %>%
  mutate(
    ci_lo = beta - 1.96 * se,
    ci_hi = beta + 1.96 * se,
    pval  = 2 * pnorm(-abs(beta / se)),
    plabel = sapply(pval, fmt_p),
    nlabel = ifelse(is.na(n) | n == "", "", paste0("N=", formatC(as.integer(n), format = "d", big.mark = ",")))
  )

make_forest_plot <- function(df_hit, cell_type_label, marker_label,
                             dim_col, palette, dim_label_str, meta_i2) {
  df <- df_hit %>%
    mutate(group_col = .data[[dim_col]]) %>%
    arrange(group_col, cohort)

  # Add meta row (diamond)
  meta_beta <- weighted.mean(df$beta, 1 / df$se^2, na.rm = TRUE)
  meta_se   <- sqrt(1 / sum(1 / df$se^2, na.rm = TRUE))
  meta_row  <- tibble(
    cohort    = paste0("Meta (I\u00b2=", ifelse(is.na(meta_i2), "?", round(meta_i2, 1)), "%)"),
    beta      = meta_beta, se = meta_se,
    ci_lo     = meta_beta - 1.96 * meta_se,
    ci_hi     = meta_beta + 1.96 * meta_se,
    pval      = 2 * pnorm(-abs(meta_beta / meta_se)),
    plabel    = fmt_p(2 * pnorm(-abs(meta_beta / meta_se))),
    nlabel    = "",
    group_col = "META",
    is_meta   = TRUE
  )

  df <- df %>% mutate(is_meta = FALSE)
  plot_df <- bind_rows(df, meta_row)

  # Factor for y-axis order: individual cohorts (grouped) then meta at top
  cohort_order <- c(df$cohort, meta_row$cohort)
  plot_df$cohort <- factor(plot_df$cohort, levels = rev(cohort_order))

  # Colours
  all_groups <- unique(plot_df$group_col)
  missing_cols <- setdiff(all_groups, names(palette))
  extra_cols <- setNames(rep("#999999", length(missing_cols)), missing_cols)
  full_palette <- c(palette, extra_cols)
  plot_df$pt_colour <- full_palette[plot_df$group_col]
  plot_df$pt_colour[is.na(plot_df$pt_colour)] <- "#999999"

  x_min <- min(plot_df$ci_lo, na.rm = TRUE)
  x_max <- max(plot_df$ci_hi, na.rm = TRUE)

  # Separator line between last individual cohort and meta
  sep_y <- 1.5

  p <- ggplot(plot_df, aes(y = cohort, x = beta)) +
    geom_hline(yintercept = sep_y, linetype = "dotted", colour = "grey60", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_errorbarh(
      aes(xmin = ci_lo, xmax = ci_hi),
      colour = plot_df$pt_colour,
      linewidth = ifelse(plot_df$is_meta, 1.3, 0.8),
      height = 0.25
    ) +
    geom_point(
      aes(shape = is_meta, size = is_meta),
      colour = plot_df$pt_colour
    ) +
    geom_text(
      aes(x = x_min - 0.02 * (x_max - x_min), label = nlabel),
      hjust = 1, size = 3.5, colour = "grey30"
    ) +
    geom_text(
      aes(x = x_max + 0.02 * (x_max - x_min), label = paste0("P=", plabel)),
      hjust = 0, size = 3.5, colour = "grey20"
    ) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values  = c("FALSE" = 3.0, "TRUE" = 5.0), guide = "none") +
    scale_x_continuous(expand = expansion(add = c(0.35, 0.35))) +
    labs(
      title    = paste0(cell_type_label, "  \u2014  ", marker_label),
      subtitle = paste0("Colour: ", dim_label_str,
                        "   |   Overall I\u00b2 = ",
                        ifelse(is.na(meta_i2), "NA", paste0(round(meta_i2, 1), "%"))),
      x = "Beta (95% CI)",
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 10, colour = "grey40"),
      axis.text.y        = element_text(size = 10),
      panel.grid.major.x = element_line(colour = "grey92"),
      plot.margin        = margin(15, 100, 15, 15)
    )

  p
}

for (k in seq_len(n_plots)) {
  ct  <- hit_keys$cell_type[k]
  mrk <- hit_keys$marker[k]
  i2  <- hit_keys$meta_HetISq[k]

  df_hit <- eff %>% filter(cell_type == ct, marker == mrk)
  if (nrow(df_hit) < 2) next

  # Safe filename
  safe_mrk <- gsub(":", "_", mrk)
  safe_ct  <- gsub("[^A-Za-z0-9._-]", "_", ct)

  # One PDF with three pages (one per dimension)
  pdf_path <- file.path(args$forest_dir,
                        paste0(safe_ct, "_", safe_mrk, "_forest.pdf"))
  pdf(pdf_path, width = 11, height = max(5, 0.5 * nrow(df_hit) + 2.5))
  for (dim_col in names(dim_palette)) {
    p <- make_forest_plot(df_hit, ct, mrk,
                          dim_col, dim_palette[[dim_col]],
                          dim_labels[dim_col], i2)
    print(p)
  }
  dev.off()
}
message("  Forest plots written to: ", args$forest_dir)

# ---------------------------------------------------------------------------
# 2. Summary panel — I² violin plots per dimension group
# ---------------------------------------------------------------------------
message("Generating I\u00b2 summary violin plots...")

for (dim_col in names(dim_palette)) {
  dim_data <- ovrl %>% filter(dimension == dim_col)
  if (nrow(dim_data) == 0) next

  palette <- dim_palette[[dim_col]]
  all_groups <- unique(dim_data$group)
  missing_g  <- setdiff(all_groups, names(palette))
  extra_g    <- setNames(rep("#999999", length(missing_g)), missing_g)
  full_pal   <- c(palette, extra_g)

  # Order groups by median I²
  group_order <- dim_data %>%
    group_by(group) %>%
    summarise(med = median(median_HetISq, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(group)

  dim_data$group <- factor(dim_data$group, levels = group_order)

  p_violin <- ggplot(dim_data, aes(x = group, y = median_HetISq, fill = group)) +
    geom_boxplot(width = 0.6, outlier.size = 1.5, outlier.alpha = 0.5) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = full_pal, guide = "none") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(
      title    = paste0("Heterogeneity (I\u00b2) by ", dim_labels[dim_col]),
      subtitle = "Median I\u00b2 across cell types at suggestive hits (p < 1e-5) — dashed line at 50%",
      x        = dim_labels[dim_col],
      y        = "Median I\u00b2 (%)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold"),
      axis.text.x   = element_text(angle = 30, hjust = 1, size = 11)
    )

  out_path <- file.path(args$output_dir,
                        paste0("i2_by_", dim_col, ".pdf"))
  ggsave(out_path, plot = p_violin, width = 8, height = 5, bg = "white")
  message("  Saved: ", out_path)
}

# ---------------------------------------------------------------------------
# 3. Four-panel combined summary (all three dimensions + overall)
# ---------------------------------------------------------------------------
message("Generating combined four-panel summary...")

panel_list <- list()
for (dim_col in names(dim_palette)) {
  dim_data <- ovrl %>% filter(dimension == dim_col)
  if (nrow(dim_data) == 0) next

  palette <- dim_palette[[dim_col]]
  all_groups <- unique(dim_data$group)
  missing_g  <- setdiff(all_groups, names(palette))
  extra_g    <- setNames(rep("#999999", length(missing_g)), missing_g)
  full_pal   <- c(palette, extra_g)

  group_order <- dim_data %>%
    group_by(group) %>%
    summarise(med = median(median_HetISq, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(group)
  dim_data$group <- factor(dim_data$group, levels = group_order)

  panel_list[[dim_col]] <- ggplot(dim_data,
                                  aes(x = group, y = median_HetISq, fill = group)) +
    geom_boxplot(width = 0.55, outlier.size = 1.2, outlier.alpha = 0.4) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
    scale_fill_manual(values = full_pal, guide = "none") +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    labs(title = dim_labels[dim_col], x = NULL, y = "Median I\u00b2 (%)") +
    theme_classic(base_size = 11) +
    theme(
      plot.title  = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9)
    )
}

# Arrange with patchwork if available, else save separately
if (length(panel_list) >= 2) {
  tryCatch({
    library(patchwork)
    combined <- wrap_plots(panel_list, ncol = 3) +
      plot_annotation(
        title    = "Heterogeneity (I\u00b2) breakdown across biological and technical dimensions",
        subtitle = "Box plots show median I\u00b2 across all cell types at suggestive hits (p < 1e-5)\nDashed line = I\u00b2 50% threshold",
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, colour = "grey40")
        )
      )
    ggsave(file.path(args$output_dir, "i2_summary_panel.pdf"),
           plot = combined, width = 16, height = 5.5, bg = "white")
    message("  Saved combined panel: i2_summary_panel.pdf")
  }, error = function(e) {
    message("  patchwork not available; panels saved individually.")
  })
}

# ---------------------------------------------------------------------------
# 4. Heatmap — cell type × cohort-group, colour = median I² at suggestive hits
# ---------------------------------------------------------------------------
message("Generating grouped I\u00b2 heatmap...")

for (dim_col in names(dim_palette)) {
  heat_data <- ovrl %>%
    filter(dimension == dim_col) %>%
    select(cell_type, group, median_HetISq) %>%
    mutate(cell_type = factor(cell_type, levels = sort(unique(cell_type))))

  if (nrow(heat_data) == 0) next

  group_order_h <- heat_data %>%
    group_by(group) %>%
    summarise(m = mean(median_HetISq, na.rm = TRUE), .groups = "drop") %>%
    arrange(m) %>%
    pull(group)
  heat_data$group <- factor(heat_data$group, levels = group_order_h)

  p_heat <- ggplot(heat_data, aes(x = group, y = cell_type, fill = median_HetISq)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = round(median_HetISq, 0)), size = 3, colour = "grey20") +
    scale_fill_gradient2(
      low      = "#F7FBFF",
      mid      = "#6BAED6",
      high     = "#08306B",
      midpoint = 50,
      limits   = c(0, 100),
      name     = "Median\nI\u00b2 (%)"
    ) +
    labs(
      title    = paste0("Median I\u00b2 by cell type and ", dim_labels[dim_col]),
      subtitle = "At suggestive hits (p < 1e-5) in the 15-cohort meta-analysis",
      x        = dim_labels[dim_col],
      y        = "Cell type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      axis.text.x   = element_text(angle = 35, hjust = 1, size = 10),
      axis.text.y   = element_text(size = 10),
      panel.grid    = element_blank()
    )

  out_h <- file.path(args$output_dir,
                     paste0("i2_heatmap_by_", dim_col, ".pdf"))
  ggsave(out_h, plot = p_heat,
         width  = max(6, length(unique(heat_data$group)) * 1.8 + 2),
         height = max(5, length(unique(heat_data$cell_type)) * 0.45 + 2),
         bg = "white")
  message("  Saved: ", out_h)
}

# ---------------------------------------------------------------------------
# 5. Within-group I² per top hit (group_i2_summary) — dot plot
# ---------------------------------------------------------------------------
message("Generating within-group I\u00b2 dot plot for top hits...")

for (dim_col in names(dim_palette)) {
  hit_data <- grp %>%
    filter(dimension == dim_col, !is.na(within_group_I2)) %>%
    mutate(
      hit_label = paste0(cell_type, "\n", marker),
      is_high   = within_group_I2 >= 50
    )

  if (nrow(hit_data) == 0) next

  palette   <- dim_palette[[dim_col]]
  all_groups <- unique(hit_data$group)
  missing_g  <- setdiff(all_groups, names(palette))
  extra_g    <- setNames(rep("#999999", length(missing_g)), missing_g)
  full_pal   <- c(palette, extra_g)

  # Limit to top 30 hits for readability
  top_markers <- hit_data %>%
    filter(is_gw == 1) %>%
    distinct(cell_type, marker) %>%
    slice_head(n = 30)
  if (nrow(top_markers) == 0) {
    top_markers <- hit_data %>% distinct(cell_type, marker) %>% slice_head(n = 30)
  }
  hit_data_sub <- hit_data %>%
    semi_join(top_markers, by = c("cell_type", "marker"))

  if (nrow(hit_data_sub) == 0) next

  p_dot <- ggplot(hit_data_sub,
                  aes(x = within_group_I2, y = reorder(hit_label, within_group_I2),
                      colour = group, shape = factor(is_gw))) +
    geom_vline(xintercept = 50, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_manual(values = full_pal, name = dim_labels[dim_col]) +
    scale_shape_manual(values = c("0" = 1, "1" = 16),
                       labels = c("0" = "Suggestive", "1" = "Genome-wide"),
                       name = "Significance") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(
      title    = paste0("Within-group I\u00b2 at top hits — ", dim_labels[dim_col]),
      subtitle = "Filled circle = genome-wide significant (p < 5e-8); open = suggestive",
      x        = "Within-group I\u00b2 (%)",
      y        = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      axis.text.y   = element_text(size = 8),
      legend.position = "right"
    )

  out_d <- file.path(args$output_dir,
                     paste0("within_group_i2_tophits_", dim_col, ".pdf"))
  ggsave(out_d, plot = p_dot,
         width  = 10,
         height = max(5, nrow(top_markers) * 0.4 + 2.5),
         bg = "white")
  message("  Saved: ", out_d)
}

message("All plots complete.")
