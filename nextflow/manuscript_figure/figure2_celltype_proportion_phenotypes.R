# =============================================================================
# Figure 2. Cell-type proportion phenotypes across bulk RNA-seq cohorts
#
# Panels:
#   A — Estimated cell-type composition across cohorts (dot-range plot)
#   B — Donor-level variability of cell-type estimates (IQR lollipop)
#   C — Cohort-level differences in cell-type estimates (z-score heatmap)
#   D — Bulk-to-snRNA validation across available cohorts (dot plot)
#
# Outputs saved to OUT_DIR:
#   bulk_celltype_proportions_long.tsv  (intermediate)
#   figure2_celltype_proportion_phenotypes.png
#   figure2_celltype_proportion_phenotypes.pdf
# =============================================================================

# =============================================================================
# SECTION 1 — File paths (edit here)
# =============================================================================

RESULTS_DIR <- "/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow/results"
OUT_DIR     <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"

# 15 canonical cohorts used in meta-analysis
COHORT_LIST <- c(
  "ROSMAP", "ROSMAP_array",
  "MSBB", "Mayo",
  "CMC_MSSM", "CMC_PENN", "CMC_PITT",
  "GTEx_v10", "NABEC", "GVEX",
  "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
  "AMP_AD_Mayo", "AMP_AD_Rush"
)

# Per-cohort proportion file (relative to results/<COHORT>/)
PROP_FILENAME <- "cell_proportions.csv"

# Validation file — Format A or Format B (set NULL to skip Panel D).
# Format A columns: validation_cohort [optional], cell_type, metric, value, n
#          OR:      cell_type, [sn_cell_type], n, pearson_r, [spearman_r], [p_value]
# Format B columns: sample_id, validation_cohort, cell_type, bulk_proportion, snrna_proportion
VALID_FILE <- file.path(RESULTS_DIR,
  "ROSMAP/scatter_mgp_vs_snrnaseq/scatter_mgp_vs_snrnaseq_correlations.tsv")
# Label shown in Panel D if only one cohort's data is in the file
VALID_DEFAULT_COHORT <- "ROSMAP"

# Figure dimensions (inches) and resolution
FIG_WIDTH  <- 14
FIG_HEIGHT <- 10
FIG_DPI    <- 300

# =============================================================================
# SECTION 2 — Package loading
# =============================================================================

required_pkgs <- c("tidyverse", "ggplot2", "cowplot", "scales",
                   "viridis", "RColorBrewer", "readr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  Installing missing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# janitor is optional — base-R fallback provided below
HAS_JANITOR <- requireNamespace("janitor", quietly = TRUE)
if (HAS_JANITOR) {
  suppressPackageStartupMessages(library(janitor))
} else {
  message("Note: package 'janitor' not found; using base-R fallback for name cleaning.")
}

clean_names_fn <- function(df) {
  if (HAS_JANITOR) return(janitor::clean_names(df))
  nms <- tolower(gsub("[^a-zA-Z0-9]+", "_",
                       gsub("([a-z])([A-Z])", "\\1_\\2", names(df))))
  names(df) <- sub("^_+|_+$", "", nms)
  df
}
make_clean_names_fn <- function(x) {
  if (HAS_JANITOR) return(janitor::make_clean_names(x))
  sub("^_+|_+$", "",
      tolower(gsub("[^a-zA-Z0-9]+", "_",
                   gsub("([a-z])([A-Z])", "\\1_\\2", x))))
}

cat("Packages loaded.\n")

# =============================================================================
# SECTION 3 — Helper functions
# =============================================================================

read_tabular <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv")            readr::read_csv(path, show_col_types = FALSE)
  else if (ext %in% c("tsv", "txt")) readr::read_tsv(path, show_col_types = FALSE)
  else stop("Unsupported extension '", ext, "': ", path)
}

make_placeholder <- function(msg, title = "") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 4, hjust = 0.5, vjust = 0.5,
             color = "grey40", fontface = "italic") +
    labs(title = title) +
    theme_void() +
    theme(
      panel.border    = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      plot.title      = element_text(size = 11, face = "bold", margin = margin(b = 6)),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

# =============================================================================
# SECTION 4 — Load and clean bulk proportions
# =============================================================================

cat("\n--- Loading per-cohort cell-type proportions ---\n")

prop_list <- lapply(COHORT_LIST, function(cohort) {
  path <- file.path(RESULTS_DIR, cohort, PROP_FILENAME)
  if (!file.exists(path)) {
    warning("File not found, skipping cohort '", cohort, "': ", path)
    return(NULL)
  }
  df <- read_tabular(path)
  df <- clean_names_fn(df)
  sid_col <- grep("^specimen_?id$|^sample_?id$", names(df),
                  value = TRUE, ignore.case = TRUE)[1]
  if (is.na(sid_col)) {
    warning("No specimen/sample ID column in '", cohort, "'. Skipping.")
    return(NULL)
  }
  df <- dplyr::rename(df, sample_id = !!sid_col)
  df$cohort <- cohort
  df
})

prop_list <- Filter(Negate(is.null), prop_list)
if (length(prop_list) == 0) stop("No proportion files loaded. Check RESULTS_DIR and COHORT_LIST.")

prop_wide <- dplyr::bind_rows(prop_list)
cat("  Loaded", nrow(prop_wide), "donors across", length(prop_list), "cohorts\n")

cell_type_cols <- setdiff(names(prop_wide), c("sample_id", "cohort"))

prop_long <- prop_wide %>%
  tidyr::pivot_longer(cols = all_of(cell_type_cols),
                      names_to  = "cell_type",
                      values_to = "proportion") %>%
  dplyr::filter(!is.na(sample_id), !is.na(cohort),
                !is.na(cell_type), !is.na(proportion))

cat("  Long-format rows:", nrow(prop_long), "\n")

# Proportion scale detection
max_prop <- max(prop_long$proportion, na.rm = TRUE)
has_negatives <- any(prop_long$proportion < 0, na.rm = TRUE)

if (max_prop <= 1 && !has_negatives) {
  prop_long$proportion_percent <- prop_long$proportion * 100
  PROP_AXIS_LABEL <- "Bulk-derived cell-type proportion (%)"
  cat("  Fractions detected → multiplied by 100\n")
} else {
  prop_long$proportion_percent <- prop_long$proportion
  PROP_AXIS_LABEL <- "MGP estimated proportion score"
  cat("  Non-fraction scale (or negative values); used as-is\n")
  cat("  Value range:", round(min(prop_long$proportion, na.rm = TRUE), 2),
      "to", round(max_prop, 2), "\n")
}

# Winsorized proportion (for any direct distribution plots)
q01 <- quantile(prop_long$proportion_percent, 0.01, na.rm = TRUE)
q99 <- quantile(prop_long$proportion_percent, 0.99, na.rm = TRUE)
prop_long$proportion_plot <- pmin(pmax(prop_long$proportion_percent, q01), q99)

# Save intermediate long-format file
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_long_path <- file.path(OUT_DIR, "bulk_celltype_proportions_long.tsv")
readr::write_tsv(prop_long, out_long_path)
cat("  Saved:", out_long_path, "\n")

# =============================================================================
# SECTION 5 — Cell-type order and color palette
# =============================================================================

CELLTYPE_ORDER <- c(
  "oligodendrocyte", "opc", "astrocyte", "microglia",
  "endothelial", "pericyte", "vlmc",
  "it", "l4_it", "l5_6_np", "l5_et", "l6_ct", "l5_6_it_car3", "l6b",
  "pvalb", "sst", "vip", "lamp5", "pax6"
)

present_types  <- unique(prop_long$cell_type)
ordered_types  <- CELLTYPE_ORDER[CELLTYPE_ORDER %in% present_types]
extra_types    <- setdiff(present_types, ordered_types)
ordered_types  <- c(ordered_types, extra_types)

prop_long$cell_type <- factor(prop_long$cell_type, levels = ordered_types)
prop_long$cohort    <- factor(prop_long$cohort,    levels = COHORT_LIST)

cat("\n--- Cell types (", length(ordered_types), "):",
    paste(ordered_types, collapse = ", "), "\n")

# Display labels for axes (replace underscores, uppercase)
ct_display <- setNames(
  gsub("_", " ", toupper(ordered_types)),
  ordered_types
)

# Color palette: glia = Blues, excitatory = Oranges, inhibitory = Greens
n_glia <- 7; n_exc <- 7; n_inh <- 5
pal_glia <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[4:9])(n_glia)
pal_exc  <- colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges")[4:9])(n_exc)
pal_inh  <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[4:9])(n_inh)

all_colors <- setNames(
  c(pal_glia, pal_exc, pal_inh)[seq_len(length(ordered_types))],
  ordered_types
)

# Common theme elements reused across panels
base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(size = 11, face = "bold"),
    plot.background    = element_rect(fill = "white", color = NA),
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 8, face = "bold"),
    axis.text          = element_text(size = 8),
    axis.title         = element_text(size = 9)
  )

# =============================================================================
# SECTION 6 — Panel A: Cross-cohort dot-range plot
# =============================================================================

cat("\n--- Building Panel A ---\n")

# Step 1: cohort × cell_type summaries (used in Panel A and Panel C)
cohort_medians <- prop_long %>%
  dplyr::group_by(cohort, cell_type) %>%
  dplyr::summarise(cohort_median = median(proportion_percent, na.rm = TRUE),
                   .groups = "drop")

# cohort_means used by Panel C z-score heatmap
cohort_means <- prop_long %>%
  dplyr::group_by(cohort, cell_type) %>%
  dplyr::summarise(mean_est = mean(proportion_percent, na.rm = TRUE),
                   .groups = "drop")

# Step 2: overall per-donor summaries per cell_type
# Use donor-level q25/q75 for the thick interval (shows typical donor range)
celltype_donor_summary <- prop_long %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    overall_median = median(proportion_percent, na.rm = TRUE),
    q25_donor      = quantile(proportion_percent, 0.25, na.rm = TRUE),
    q75_donor      = quantile(proportion_percent, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = ordered_types))

# Step 3: cross-cohort range of cohort medians (thin interval)
cohort_range_summary <- cohort_medians %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    rng_min = min(cohort_median, na.rm = TRUE),
    rng_max = max(cohort_median, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = ordered_types))

panel_A_data <- dplyr::left_join(celltype_donor_summary, cohort_range_summary,
                                  by = "cell_type")

cohort_med_plot <- cohort_medians %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = ordered_types))

panel_A <- ggplot(panel_A_data, aes(y = cell_type, color = cell_type)) +
  # Individual cohort median dots (faint, show cohort variation)
  geom_point(data = cohort_med_plot,
             aes(x = cohort_median),
             shape = 16, size = 1.0, alpha = 0.30) +
  # Cross-cohort range of cohort medians (thin whisker)
  geom_linerange(aes(xmin = rng_min, xmax = rng_max),
                 linewidth = 0.5, alpha = 0.55) +
  # Donor-level IQR (thick bar)
  geom_linerange(aes(xmin = q25_donor, xmax = q75_donor),
                 linewidth = 2.0, alpha = 0.80) +
  # Overall median point
  geom_point(aes(x = overall_median, fill = cell_type),
             shape = 21, size = 2.8, color = "white", stroke = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.35) +
  scale_color_manual(values = all_colors, guide = "none") +
  scale_fill_manual( values = all_colors, guide = "none") +
  scale_y_discrete(limits = rev(ordered_types), labels = ct_display) +
  labs(
    title = "A. Estimated cell-type composition across cohorts",
    x     = PROP_AXIS_LABEL,
    y     = NULL
  ) +
  base_theme +
  theme(panel.grid.major.y = element_blank())

# =============================================================================
# SECTION 7 — Panel B: Donor-level variability (IQR lollipop)
# =============================================================================

cat("--- Building Panel B ---\n")

variability <- prop_long %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    median_est = median(proportion_percent, na.rm = TRUE),
    iqr        = IQR(proportion_percent, na.rm = TRUE),
    mad        = mad(proportion_percent, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = ordered_types))

panel_B <- ggplot(variability,
                  aes(y = cell_type, x = iqr,
                      color = cell_type, fill = cell_type)) +
  geom_segment(aes(xend = 0, yend = cell_type),
               linewidth = 0.5, color = "grey70") +
  geom_point(shape = 21, size = 3.0, color = "white", stroke = 0.4) +
  scale_color_manual(values = all_colors, guide = "none") +
  scale_fill_manual( values = all_colors, guide = "none") +
  scale_y_discrete(limits = rev(ordered_types), labels = ct_display) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "B. Donor-level variability of cell-type estimates",
    x     = "Donor-level variability (IQR)",
    y     = NULL
  ) +
  base_theme +
  theme(panel.grid.major.y = element_blank())

# =============================================================================
# SECTION 8 — Panel C: Z-score heatmap of cohort differences
# =============================================================================

cat("--- Building Panel C ---\n")

# Per cell_type: z-score cohort means relative to cross-cohort distribution
cohort_z <- cohort_means %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(
    ct_mean = mean(mean_est, na.rm = TRUE),
    ct_sd   = sd(mean_est,   na.rm = TRUE),
    z_score = dplyr::if_else(ct_sd > 0, (mean_est - ct_mean) / ct_sd, 0)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    z_clip    = pmax(pmin(z_score, 2.5), -2.5),
    cell_type = factor(cell_type, levels = ordered_types),
    cohort    = factor(cohort,    levels = COHORT_LIST)
  )

# Check: if all z-values are nearly flat (max abs < 0.5), auto-focus on most variable CTs
variable_cts <- cohort_z %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(range_z = max(z_clip, na.rm=TRUE) - min(z_clip, na.rm=TRUE),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(range_z))

max_range <- max(variable_cts$range_z, na.rm = TRUE)
if (max_range < 0.5) {
  top_cts <- variable_cts$cell_type[1:min(10, nrow(variable_cts))]
  cohort_z <- dplyr::filter(cohort_z, cell_type %in% top_cts)
  message("Panel C: all z-values near 0; showing top 10 most variable cell types.")
}

# Order cohorts by hierarchical clustering for visual grouping
ct_z_wide <- cohort_z %>%
  dplyr::select(cohort, cell_type, z_clip) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = z_clip, values_fill = 0)
cohort_mat <- as.matrix(dplyr::select(ct_z_wide, -cohort))
rownames(cohort_mat) <- ct_z_wide$cohort

hc <- hclust(dist(cohort_mat), method = "ward.D2")
cohort_hc_order <- rownames(cohort_mat)[hc$order]

cohort_z$cohort <- factor(cohort_z$cohort, levels = cohort_hc_order)

panel_C <- ggplot(cohort_z,
                  aes(x = cell_type, y = cohort, fill = z_clip)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low      = "#2c7bb6",
    mid      = "#ffffbf",
    high     = "#d7191c",
    midpoint = 0,
    limits   = c(-2.5, 2.5),
    name     = "Relative\ncohort z",
    labels   = scales::label_number(accuracy = 0.1)
  ) +
  scale_x_discrete(labels = ct_display) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(
    title = "C. Cohort-level differences in cell-type estimates",
    x     = NULL,
    y     = NULL
  ) +
  base_theme +
  theme(
    axis.text.x       = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 7.5),
    panel.grid        = element_blank(),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.3, "cm")
  )

# =============================================================================
# SECTION 9 — Panel D: Bulk-to-snRNA validation
# =============================================================================

cat("--- Building Panel D ---\n")

panel_D <- tryCatch({

  if (is.null(VALID_FILE) || !file.exists(VALID_FILE)) {
    stop("Validation file not provided or not found.")
  }

  val_raw <- read_tabular(VALID_FILE)
  val_raw <- clean_names_fn(val_raw)

  # ---- Detect Format B: has bulk_proportion and snrna_proportion ----
  if ("bulk_proportion" %in% names(val_raw) &&
      "snrna_proportion" %in% names(val_raw)) {

    # Format B — compute correlations per validation_cohort × cell_type
    vcohort_col <- if ("validation_cohort" %in% names(val_raw)) "validation_cohort" else NULL

    val_corr <- val_raw %>%
      { if (!is.null(vcohort_col)) dplyr::group_by(., dplyr::across(dplyr::all_of(c(vcohort_col, "cell_type"))))
        else dplyr::group_by(., cell_type) } %>%
      dplyr::summarise(
        pearson_r  = cor(bulk_proportion, snrna_proportion,
                         method = "pearson",  use = "pairwise.complete.obs"),
        spearman_r = cor(bulk_proportion, snrna_proportion,
                         method = "spearman", use = "pairwise.complete.obs"),
        n_donors   = dplyr::n(),
        .groups    = "drop"
      )

    if (is.null(vcohort_col)) val_corr$validation_cohort <- VALID_DEFAULT_COHORT
    corr_col    <- "spearman_r"
    axis_label  <- "Spearman \u03c1 (bulk vs snRNA-seq)"

  } else {
    # ---- Format A — precomputed correlations ----
    # Supported layouts:
    #   (a) validation_cohort, cell_type, metric, value, n
    #   (b) cell_type, [sn_cell_type], n, pearson_r, [spearman_r], [p_value]
    has_metric_col   <- "metric"    %in% names(val_raw)
    has_spearman_col <- "spearman_r" %in% names(val_raw)
    has_pearson_col  <- "pearson_r"  %in% names(val_raw)

    if (has_metric_col) {
      # Layout (a): long format
      val_corr <- val_raw %>%
        dplyr::filter(metric %in% c("spearman_r", "pearson_r")) %>%
        tidyr::pivot_wider(names_from = metric, values_from = value)
    } else {
      # Layout (b): wide format — keep pearson_r / spearman_r columns
      val_corr <- val_raw
    }

    # Add validation_cohort column if missing
    if (!"validation_cohort" %in% names(val_corr)) {
      val_corr$validation_cohort <- VALID_DEFAULT_COHORT
    }

    # Standardise n column
    if ("n" %in% names(val_corr)) val_corr$n_donors <- val_corr$n

    # Choose primary correlation metric
    if ("spearman_r" %in% names(val_corr)) {
      corr_col   <- "spearman_r"
      axis_label <- "Spearman \u03c1 (bulk vs snRNA-seq)"
    } else if ("pearson_r" %in% names(val_corr)) {
      corr_col   <- "pearson_r"
      axis_label <- "Pearson r (bulk vs snRNA-seq)"
    } else {
      stop("Format A validation file has no pearson_r or spearman_r column.")
    }
  }

  # ---- Harmonise cell_type names ----
  val_corr$cell_type_clean <- make_clean_names_fn(val_corr$cell_type)

  # ---- Compute p_value significance labels ----
  if ("p_value" %in% names(val_corr)) {
    val_corr$sig <- dplyr::case_when(
      is.na(val_corr$p_value) ~ "",
      val_corr$p_value < 0.001 ~ "***",
      val_corr$p_value < 0.01  ~ "**",
      val_corr$p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  } else {
    val_corr$sig <- ""
  }

  # ---- Order cell types by median correlation across cohorts ----
  ct_median_corr <- val_corr %>%
    dplyr::group_by(cell_type_clean) %>%
    dplyr::summarise(med_corr = median(.data[[corr_col]], na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(med_corr)

  val_corr$cell_type_clean <- factor(
    val_corr$cell_type_clean,
    levels = ct_median_corr$cell_type_clean
  )

  # Multiple validation cohorts get distinct colors
  vcohorts <- unique(val_corr$validation_cohort)
  if (length(vcohorts) == 1) {
    vcohort_colors <- setNames("#2c7bb6", vcohorts)
  } else {
    vcohort_colors <- setNames(
      RColorBrewer::brewer.pal(max(3, length(vcohorts)), "Set2")[seq_along(vcohorts)],
      vcohorts
    )
  }

  val_corr$validation_cohort <- factor(val_corr$validation_cohort, levels = vcohorts)

  x_lo <- min(c(-0.05, val_corr[[corr_col]]), na.rm = TRUE) - 0.02

  p <- ggplot(val_corr,
              aes(x = .data[[corr_col]], y = cell_type_clean,
                  color = validation_cohort)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_segment(aes(xend = 0, yend = cell_type_clean, group = validation_cohort),
                 linewidth = 0.4, color = "grey75",
                 position = position_dodge(width = 0.6)) +
    geom_point(size = 3.2,
               position = position_dodge(width = 0.6)) +
    geom_text(aes(label = sig),
              hjust = -0.35, vjust = 0.5, size = 3.0, color = "grey25",
              show.legend = FALSE,
              position = position_dodge(width = 0.6)) +
    scale_color_manual(values = vcohort_colors,
                       name = "Validation\ncohort") +
    scale_y_discrete(labels = function(x) gsub("_", " ", toupper(x))) +
    scale_x_continuous(limits = c(x_lo, NA),
                       labels = scales::label_number(accuracy = 0.01)) +
    labs(
      title = "D. Bulk-to-snRNA validation across available cohorts",
      x     = axis_label,
      y     = NULL
    ) +
    base_theme +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position    = "right",
      legend.key.size    = unit(0.45, "cm")
    )

  # Facet if more than 3 validation cohorts
  if (length(vcohorts) > 3) {
    p <- p + facet_wrap(~ validation_cohort, nrow = 1) +
      theme(legend.position = "none")
  }
  p

}, error = function(e) {
  message("Panel D: ", conditionMessage(e), "\n  -> Placeholder.")
  make_placeholder("Validation data not provided.", "D. Bulk-to-snRNA validation")
})

# =============================================================================
# SECTION 10 — Combine and save
# =============================================================================

cat("\n--- Combining panels ---\n")

top_row <- cowplot::plot_grid(
  panel_A, panel_B,
  ncol       = 2,
  rel_widths = c(1.15, 1),
  align      = "h",
  axis       = "tb"
)

bottom_row <- cowplot::plot_grid(
  panel_C, panel_D,
  ncol       = 2,
  rel_widths = c(1.3, 1),
  align      = "h",
  axis       = "tb"
)

fig2 <- cowplot::plot_grid(
  top_row, bottom_row,
  nrow        = 2,
  rel_heights = c(1, 1)
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

out_png <- file.path(OUT_DIR, "figure2_celltype_proportion_phenotypes.png")
out_pdf <- file.path(OUT_DIR, "figure2_celltype_proportion_phenotypes.pdf")

ggsave(out_png, fig2, width = FIG_WIDTH, height = FIG_HEIGHT,
       dpi = FIG_DPI, bg = "white")
cat("Saved PNG:", out_png, "\n")

ggsave(out_pdf, fig2, width = FIG_WIDTH, height = FIG_HEIGHT, bg = "white")
cat("Saved PDF:", out_pdf, "\n")

cat("\nDone.\n")
