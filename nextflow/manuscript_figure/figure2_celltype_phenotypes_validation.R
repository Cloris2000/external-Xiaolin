# =============================================================================
# Figure 2. Bulk-derived cell-type proportion phenotypes and snRNA-seq validation
#
# Panels:
#   A — Estimated cell-type composition across cohorts (dot-range plot)
#   B — Validation scatter plots for prioritized cell types (VIP, L5.6.IT.Car3,
#         Microglia). Requires matched donor-level validation data (Format A).
#         Shows placeholder until VALID_PAIRED_FILE is provided.
#   C — Cohort-level variation in prioritized cell types (dot + IQR per cohort)
#   D — All-cell-type bulk-to-snRNA validation summary (dot/lollipop)
#         Uses precomputed correlations file by default.
#
# Output files:
#   figure2_celltype_phenotypes_validation.png
#   figure2_celltype_phenotypes_validation.pdf
# =============================================================================

# =============================================================================
# SECTION 1 — File paths and parameters (edit here)
# =============================================================================

RESULTS_DIR <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results"
OUT_DIR     <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure"

# 15 canonical cohorts
COHORT_LIST <- c(
  "ROSMAP", "ROSMAP_array",
  "MSBB", "Mayo",
  "CMC_MSSM", "CMC_PENN", "CMC_PITT",
  "GTEx_v10", "NABEC", "GVEX",
  "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
  "AMP_AD_Mayo", "AMP_AD_Rush"
)
PROP_FILENAME <- "cell_proportions.csv"

# Cell types to highlight across all panels
FOCUS_CELLTYPES <- c("vip", "l5_6_it_car3", "microglia")

# --------------------------------------------------------------------------
# Panel B — matched donor-level validation file (Format A).
# Columns: sample_id, validation_cohort, cell_type, bulk_proportion,
#          snrna_proportion
# Set to NULL to show a placeholder for Panel B.
# --------------------------------------------------------------------------
VALID_PAIRED_FILE <- NULL

# --------------------------------------------------------------------------
# Panel D — precomputed correlation summary (Format B).
# Supported column layouts:
#   (a) cell_type, [sn_cell_type], n, pearson_r, [spearman_r], [p_value]
#       — single-cohort precomputed table (validation_cohort added below)
#   (b) validation_cohort, cell_type, metric, value, n
#       — long/tidy format with explicit cohort column
# Set to NULL to show a placeholder for Panel D.
# --------------------------------------------------------------------------
VALID_PRECOMP_FILE <- file.path(RESULTS_DIR,
  "ROSMAP/scatter_mgp_vs_snrnaseq/scatter_mgp_vs_snrnaseq_correlations.tsv")
VALID_PRECOMP_COHORT_LABEL <- "ROSMAP"   # used when file has no cohort column

# Figure dimensions
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
    message("  Installing: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

HAS_JANITOR <- requireNamespace("janitor", quietly = TRUE)
if (HAS_JANITOR) {
  suppressPackageStartupMessages(library(janitor))
} else {
  message("Note: 'janitor' not found; using base-R fallback.")
}

clean_names_fn <- function(df) {
  if (HAS_JANITOR) return(janitor::clean_names(df))
  nms <- sub("^_+|_+$", "",
             tolower(gsub("[^a-zA-Z0-9]+", "_",
                          gsub("([a-z])([A-Z])", "\\1_\\2", names(df)))))
  names(df) <- nms
  df
}
make_clean_fn <- function(x) {
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
  if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported extension '", ext, "': ", path)
  }
}

make_placeholder <- function(msg, title = "") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 4, hjust = 0.5, vjust = 0.5,
             color = "grey45", fontface = "italic") +
    labs(title = title) +
    theme_void() +
    theme(panel.border    = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
          plot.title      = element_text(size = 11, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

spearman_label <- function(x, y) {
  ct    <- suppressWarnings(cor.test(x, y, method = "spearman"))
  r     <- round(ct$estimate, 2)
  p_str <- if (ct$p.value < 0.001) "p < 0.001" else paste0("p = ", round(ct$p.value, 3))
  paste0("r[s] = ", r, "\n", p_str)
}

# =============================================================================
# SECTION 4 — Load bulk proportions
# =============================================================================

cat("\n--- Loading bulk cell-type proportions ---\n")

prop_list <- lapply(COHORT_LIST, function(cohort) {
  path <- file.path(RESULTS_DIR, cohort, PROP_FILENAME)
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  df  <- read_tabular(path)
  df  <- clean_names_fn(df)
  sid <- grep("^specimen_?id$|^sample_?id$", names(df),
              value = TRUE, ignore.case = TRUE)[1]
  if (is.na(sid)) { warning("No ID column in ", cohort); return(NULL) }
  df  <- dplyr::rename(df, sample_id = !!sid)
  df$cohort <- cohort
  df
})
prop_list <- Filter(Negate(is.null), prop_list)

prop_wide <- dplyr::bind_rows(prop_list)
cat("  Donors:", nrow(prop_wide), "across", length(prop_list), "cohorts\n")

cell_type_cols <- setdiff(names(prop_wide), c("sample_id", "cohort"))

prop_long <- prop_wide %>%
  tidyr::pivot_longer(cols = all_of(cell_type_cols),
                      names_to  = "cell_type",
                      values_to = "proportion") %>%
  dplyr::filter(!is.na(sample_id), !is.na(cohort),
                !is.na(cell_type), !is.na(proportion))

max_prop <- max(prop_long$proportion, na.rm = TRUE)
has_neg  <- any(prop_long$proportion < 0, na.rm = TRUE)

if (max_prop <= 1 && !has_neg) {
  prop_long$proportion_pct <- prop_long$proportion * 100
  PROP_XLAB <- "Bulk-derived cell-type proportion (%)"
  cat("  Fractions detected; multiplied by 100\n")
} else {
  prop_long$proportion_pct <- prop_long$proportion
  PROP_XLAB <- "Bulk-derived cell-type proportion estimate (MGP)"
  cat("  Non-fraction/negative scale; used as-is\n")
}
prop_long$cohort <- factor(prop_long$cohort, levels = COHORT_LIST)

# =============================================================================
# SECTION 5 — Cell-type ordering, labels, and colors
# =============================================================================

CELLTYPE_ORDER <- c(
  "oligodendrocyte", "opc", "astrocyte", "microglia",
  "endothelial", "pericyte", "vlmc",
  "it", "l4_it", "l5_6_np", "l5_et", "l6_ct", "l5_6_it_car3", "l6b",
  "pvalb", "sst", "vip", "lamp5", "pax6"
)
present_ct  <- unique(prop_long$cell_type)
ordered_ct  <- c(CELLTYPE_ORDER[CELLTYPE_ORDER %in% present_ct],
                 setdiff(present_ct, CELLTYPE_ORDER))
prop_long$cell_type <- factor(prop_long$cell_type, levels = ordered_ct)

CT_DISPLAY <- setNames(toupper(gsub("_", " ", ordered_ct)), ordered_ct)
CT_DISPLAY["l5_6_it_car3"] <- "L5/6 IT Car3"
CT_DISPLAY["l5_6_np"]      <- "L5/6 NP"
CT_DISPLAY["l4_it"]        <- "L4 IT"
CT_DISPLAY["l5_et"]        <- "L5 ET"
CT_DISPLAY["l6_ct"]        <- "L6 CT"
CT_DISPLAY["l6b"]          <- "L6b"
CT_DISPLAY["it"]           <- "IT"

n_ct   <- length(ordered_ct)
pal_g  <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[3:8])(7)
pal_e  <- colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges")[3:8])(7)
pal_i  <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[3:8])(5)
ct_colors <- setNames(c(pal_g, pal_e, pal_i)[seq_len(n_ct)], ordered_ct)

y_face_A <- ifelse(rev(ordered_ct) %in% FOCUS_CELLTYPES, "bold", "plain")

base_th <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(size = 11, face = "bold"),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.text       = element_text(size = 7.5),
    legend.title      = element_text(size = 8.5, face = "bold"),
    axis.text         = element_text(size = 8),
    axis.title        = element_text(size = 9)
  )

cat("Cell types:", length(ordered_ct), "\n")

# =============================================================================
# SECTION 6 — Panel A: Cross-cohort dot-range summary
# =============================================================================

cat("\n--- Building Panel A ---\n")

cohort_med <- prop_long %>%
  dplyr::group_by(cohort, cell_type) %>%
  dplyr::summarise(coh_med = median(proportion_pct, na.rm = TRUE), .groups = "drop")

ct_summary <- cohort_med %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    ct_median = median(coh_med, na.rm = TRUE),
    q25       = quantile(coh_med, 0.25, na.rm = TRUE),
    q75       = quantile(coh_med, 0.75, na.rm = TRUE),
    rng_min   = min(coh_med, na.rm = TRUE),
    rng_max   = max(coh_med, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    cell_type = factor(cell_type, levels = ordered_ct),
    is_focus  = cell_type %in% FOCUS_CELLTYPES,
    lwd_iqr   = ifelse(is_focus, 2.4, 1.5),
    lwd_rng   = ifelse(is_focus, 0.7, 0.45),
    pt_size   = ifelse(is_focus, 3.8, 2.6)
  )

cohort_med_fac <- dplyr::mutate(cohort_med,
  cell_type = factor(cell_type, levels = ordered_ct))

panel_A <- ggplot(ct_summary, aes(y = cell_type, color = cell_type)) +
  geom_point(data   = cohort_med_fac, aes(x = coh_med),
             shape  = 16, size = 1.0, alpha = 0.22) +
  geom_linerange(aes(xmin = rng_min, xmax = rng_max, linewidth = lwd_rng),
                 alpha = 0.50) +
  geom_linerange(aes(xmin = q25, xmax = q75, linewidth = lwd_iqr),
                 alpha = 0.82) +
  geom_point(aes(x = ct_median, fill = cell_type, size = pt_size),
             shape = 21, color = "white", stroke = 0.4) +
  geom_point(data = dplyr::filter(ct_summary, is_focus),
             aes(x = ct_median),
             shape = 21, size = 4.2, fill = NA, color = "black", stroke = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.35) +
  scale_linewidth_identity() +
  scale_size_identity() +
  scale_color_manual(values = ct_colors, guide = "none") +
  scale_fill_manual( values = ct_colors, guide = "none") +
  scale_y_discrete(limits = rev(ordered_ct), labels = CT_DISPLAY) +
  labs(title = "A. Estimated cell-type composition across cohorts",
       x = PROP_XLAB, y = NULL) +
  base_th +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.5, face = y_face_A))

# =============================================================================
# SECTION 7 — Panel B: Scatter plots for focus cell types
# =============================================================================

cat("--- Building Panel B ---\n")

panel_B <- tryCatch({

  if (is.null(VALID_PAIRED_FILE) || !file.exists(VALID_PAIRED_FILE)) {
    stop("VALID_PAIRED_FILE not set or not found.")
  }

  pval_df <- read_tabular(VALID_PAIRED_FILE)
  pval_df <- clean_names_fn(pval_df)
  req_cols <- c("sample_id", "cell_type", "bulk_proportion", "snrna_proportion")
  missing  <- setdiff(req_cols, names(pval_df))
  if (length(missing) > 0) {
    stop("VALID_PAIRED_FILE missing columns: ", paste(missing, collapse = ", "))
  }
  if (!"validation_cohort" %in% names(pval_df)) {
    pval_df$validation_cohort <- "Validation cohort"
  }

  focus_data <- pval_df %>%
    dplyr::mutate(ct_clean = make_clean_fn(cell_type)) %>%
    dplyr::filter(ct_clean %in% FOCUS_CELLTYPES) %>%
    dplyr::mutate(
      cell_type_label   = factor(CT_DISPLAY[ct_clean], levels = CT_DISPLAY[FOCUS_CELLTYPES]),
      validation_cohort = factor(validation_cohort)
    )

  if (nrow(focus_data) == 0) stop("No paired data for focus cell types.")

  vcohorts <- levels(focus_data$validation_cohort)
  vcol <- if (length(vcohorts) == 1) {
    setNames("#2c7bb6", vcohorts)
  } else {
    setNames(RColorBrewer::brewer.pal(max(3, length(vcohorts)), "Set2")[seq_along(vcohorts)],
             vcohorts)
  }

  corr_labels <- focus_data %>%
    dplyr::group_by(cell_type_label) %>%
    dplyr::summarise(
      lbl   = spearman_label(bulk_proportion, snrna_proportion),
      x_pos = quantile(snrna_proportion, 0.02, na.rm = TRUE),
      y_pos = quantile(bulk_proportion,  0.98, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(focus_data, aes(x = snrna_proportion, y = bulk_proportion,
                          color = validation_cohort)) +
    geom_point(size = 1.2, alpha = 0.45) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7,
                alpha = 0.15, formula = y ~ x) +
    geom_text(data = corr_labels,
              aes(x = x_pos, y = y_pos, label = lbl),
              inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 2.8, color = "grey25",
              fontface = "italic") +
    facet_wrap(~ cell_type_label, nrow = 1, scales = "free") +
    scale_color_manual(values = vcol, name = "Cohort") +
    labs(title = "B. Validation examples for prioritized cell types",
         x = "snRNA-seq-derived proportion",
         y = "Bulk-derived estimate (MGP)") +
    base_th +
    theme(strip.text       = element_text(size = 9, face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.position  = if (length(vcohorts) == 1) "none" else "right")

}, error = function(e) {
  message("Panel B: ", conditionMessage(e))
  make_placeholder(
    paste0("Provide matched donor-level validation data\n",
           "(set VALID_PAIRED_FILE).\n\n",
           "Required columns:\n",
           "sample_id, validation_cohort, cell_type,\n",
           "bulk_proportion, snrna_proportion"),
    "B. Validation examples for prioritized cell types"
  )
})

# =============================================================================
# SECTION 8 — Panel C: Cohort-level variation in focus cell types
# =============================================================================

cat("--- Building Panel C ---\n")

focus_long <- prop_long %>%
  dplyr::filter(cell_type %in% FOCUS_CELLTYPES) %>%
  dplyr::mutate(
    cell_type_label = factor(CT_DISPLAY[as.character(cell_type)],
                             levels = CT_DISPLAY[FOCUS_CELLTYPES])
  )

cohort_sum_c <- focus_long %>%
  dplyr::group_by(cohort, cell_type_label) %>%
  dplyr::summarise(
    med = median(proportion_pct, na.rm = TRUE),
    q25 = quantile(proportion_pct, 0.25, na.rm = TRUE),
    q75 = quantile(proportion_pct, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Order cohorts by mean-of-medians
cohort_order_c <- cohort_sum_c %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(overall = mean(med, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(overall) %>%
  dplyr::pull(cohort) %>%
  as.character()

cohort_sum_c$cohort <- factor(cohort_sum_c$cohort, levels = cohort_order_c)

focus_fills  <- ct_colors[FOCUS_CELLTYPES]
names(focus_fills) <- CT_DISPLAY[FOCUS_CELLTYPES]

panel_C <- ggplot(cohort_sum_c,
                  aes(y = cohort, x = med,
                      color = cell_type_label,
                      fill  = cell_type_label)) +
  geom_segment(aes(x = q25, xend = q75, yend = cohort),
               linewidth = 1.4, alpha = 0.70) +
  geom_point(size = 2.6, shape = 21, color = "white", stroke = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey55", linewidth = 0.35) +
  facet_wrap(~ cell_type_label, nrow = 1, scales = "free_x") +
  scale_color_manual(values = focus_fills, guide = "none") +
  scale_fill_manual( values = focus_fills, guide = "none") +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(title = "C. Cohort-level variation in prioritized cell types",
       x = PROP_XLAB, y = NULL) +
  base_th +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(size = 9, face = "bold"),
    strip.background   = element_rect(fill = "grey95", color = NA),
    axis.text.y        = element_text(size = 7.5)
  )

# =============================================================================
# SECTION 9 — Panel D: All-cell-type validation summary
# =============================================================================

cat("--- Building Panel D ---\n")

panel_D <- tryCatch({

  if (is.null(VALID_PRECOMP_FILE) || !file.exists(VALID_PRECOMP_FILE)) {
    stop("VALID_PRECOMP_FILE not set or not found.")
  }

  precomp <- read_tabular(VALID_PRECOMP_FILE)
  precomp <- clean_names_fn(precomp)

  # --- Detect format ---
  # Format (b): has 'metric' and 'value' columns → pivot wide
  if ("metric" %in% names(precomp) && "value" %in% names(precomp)) {
    precomp <- precomp %>%
      tidyr::pivot_wider(names_from = metric, values_from = value)
  }

  # Add validation_cohort if missing
  if (!"validation_cohort" %in% names(precomp)) {
    precomp$validation_cohort <- VALID_PRECOMP_COHORT_LABEL
  }

  # Standardise n
  if ("n" %in% names(precomp) && !"n_donors" %in% names(precomp)) {
    precomp$n_donors <- precomp$n
  } else if (!"n_donors" %in% names(precomp)) {
    precomp$n_donors <- NA_real_
  }

  # Choose primary correlation metric
  if ("pearson_r" %in% names(precomp) && !"spearman_r" %in% names(precomp)) {
    precomp <- dplyr::rename(precomp, corr_val = pearson_r)
    D_XLAB  <- "Pearson r (bulk vs snRNA-seq)"
  } else if ("spearman_r" %in% names(precomp)) {
    precomp <- dplyr::rename(precomp, corr_val = spearman_r)
    D_XLAB  <- "Spearman rho (bulk vs snRNA-seq)"
  } else {
    stop("Precomputed file has no pearson_r or spearman_r column.")
  }

  # Clean cell-type names and match to ordered_ct
  precomp <- precomp %>%
    dplyr::mutate(ct_clean = make_clean_fn(cell_type)) %>%
    dplyr::filter(!is.na(corr_val))

  # Significance stars from p_value if available
  precomp$sig <- ""
  if ("p_value" %in% names(precomp)) {
    precomp$sig <- dplyr::case_when(
      is.na(precomp$p_value)        ~ "",
      precomp$p_value < 0.001       ~ "***",
      precomp$p_value < 0.01        ~ "**",
      precomp$p_value < 0.05        ~ "*",
      TRUE                          ~ ""
    )
  }

  # Order cell types by median correlation
  ct_order_d <- precomp %>%
    dplyr::group_by(ct_clean) %>%
    dplyr::summarise(med_r = median(corr_val, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(med_r) %>%
    dplyr::pull(ct_clean)

  precomp <- precomp %>%
    dplyr::mutate(
      ct_clean          = factor(ct_clean, levels = ct_order_d),
      is_focus          = ct_clean %in% FOCUS_CELLTYPES,
      pt_shape          = ifelse(is_focus, 23L, 21L),
      validation_cohort = factor(validation_cohort)
    )

  vcoh_d <- levels(precomp$validation_cohort)
  vcol_d <- if (length(vcoh_d) == 1) {
    setNames("#2c7bb6", vcoh_d)
  } else {
    setNames(RColorBrewer::brewer.pal(max(3, length(vcoh_d)), "Set2")[seq_along(vcoh_d)],
             vcoh_d)
  }

  x_lo <- min(c(-0.05, precomp$corr_val), na.rm = TRUE) - 0.02

  # y-axis labels with bold for focus cell types
  d_y_labels <- function(x) {
    lbl <- CT_DISPLAY[x]
    ifelse(is.na(lbl), x, lbl)
  }
  d_y_face <- ifelse(ct_order_d %in% FOCUS_CELLTYPES, "bold", "plain")

  ggplot(precomp,
         aes(x = corr_val, y = ct_clean,
             color = validation_cohort, fill = validation_cohort)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_segment(aes(xend = 0, yend = ct_clean, group = validation_cohort),
                 linewidth = 0.45, color = "grey75",
                 position = position_dodge(width = 0.65)) +
    geom_point(aes(shape = pt_shape, size = n_donors), stroke = 0.5,
               alpha = 0.9,
               position = position_dodge(width = 0.65)) +
    geom_text(aes(label = sig),
              hjust = -0.35, vjust = 0.4, size = 3.0, color = "grey25",
              show.legend = FALSE,
              position = position_dodge(width = 0.65)) +
    scale_shape_identity() +
    scale_size_continuous(range = c(2.5, 5.5), name = "Donors (n)",
                          labels = scales::label_comma()) +
    scale_color_manual(values = vcol_d, name = "Validation\ncohort") +
    scale_fill_manual( values = vcol_d, guide = "none") +
    scale_y_discrete(labels = d_y_labels) +
    scale_x_continuous(limits = c(x_lo, NA),
                       labels = scales::label_number(accuracy = 0.01)) +
    labs(title = "D. Bulk-to-snRNA validation across available cohorts",
         x = D_XLAB, y = NULL) +
    base_th +
    theme(
      panel.grid.major.y = element_blank(),
      legend.key.size    = unit(0.45, "cm"),
      axis.text.y        = element_text(size = 7.5, face = d_y_face)
    )

}, error = function(e) {
  message("Panel D: ", conditionMessage(e))
  make_placeholder("Set VALID_PRECOMP_FILE to a precomputed\ncorrelation summary (TSV/CSV).",
                   "D. Bulk-to-snRNA validation across available cohorts")
})

# =============================================================================
# SECTION 10 — Combine and save
# =============================================================================

cat("\n--- Combining panels ---\n")

top_row <- cowplot::plot_grid(
  panel_A, panel_B,
  ncol = 2, rel_widths = c(1.0, 1.5),
  align = "h", axis = "tb"
)
bot_row <- cowplot::plot_grid(
  panel_C, panel_D,
  ncol = 2, rel_widths = c(1.2, 1.0),
  align = "h", axis = "tb"
)
fig2 <- cowplot::plot_grid(top_row, bot_row, nrow = 2, rel_heights = c(1.05, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

out_png <- file.path(OUT_DIR, "figure2_celltype_phenotypes_validation.png")
out_pdf <- file.path(OUT_DIR, "figure2_celltype_phenotypes_validation.pdf")

ggsave(out_png, fig2, width = FIG_WIDTH, height = FIG_HEIGHT,
       dpi = FIG_DPI, bg = "white")
cat("Saved PNG:", out_png, "\n")

ggsave(out_pdf, fig2, width = FIG_WIDTH, height = FIG_HEIGHT, bg = "white")
cat("Saved PDF:", out_pdf, "\n")

cat("\nDone.\n")
