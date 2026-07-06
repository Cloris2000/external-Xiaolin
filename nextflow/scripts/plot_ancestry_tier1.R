#!/usr/bin/env Rscript
# Visualise Ancestry Tier 1 screening results.
#
# Classifies each hit by EUR vs HBCC stratum concordance and produces:
#   1. Pattern summary bar chart (overall + by cell type)
#   2. EUR vs HBCC meta-beta scatter (quadrants)
#   3. EUR vs HBCC -log10(p) scatter
#   4. Top-hit beta heatmap (pooled / EUR / HBCC)
#   5. Forest plots for top priority loci (from per-cohort detail)
#
# Usage:
#   Rscript plot_ancestry_tier1.R \
#     --screened-hits  results/meta_sensitivity/ancestry_tier1/ancestry_comparison_screened_hits.tsv \
#     --locus-detail   results/meta_sensitivity/ancestry_tier1/ancestry_locus_detail.tsv \
#     --output-dir     results/meta_sensitivity/ancestry_tier1/figures \
#     --forest-dir     results/meta_sensitivity/ancestry_tier1/forest_plots \
#     [--max-forest-plots 30] [--max-labels 15]

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
    screened_hits    = NULL,
    locus_detail     = NULL,
    output_dir       = "figures",
    forest_dir       = "forest_plots",
    max_forest_plots = 30L,
    max_labels       = 15L
  )
  i <- 1L
  while (i <= length(argv)) {
    flag <- argv[i]
    val  <- if (i + 1L <= length(argv)) argv[i + 1L] else NA_character_
    if      (flag == "--screened-hits")    { out$screened_hits    <- val; i <- i + 2L }
    else if (flag == "--locus-detail")     { out$locus_detail     <- val; i <- i + 2L }
    else if (flag == "--output-dir")       { out$output_dir       <- val; i <- i + 2L }
    else if (flag == "--forest-dir")       { out$forest_dir       <- val; i <- i + 2L }
    else if (flag == "--max-forest-plots") { out$max_forest_plots <- as.integer(val); i <- i + 2L }
    else if (flag == "--max-labels")       { out$max_labels       <- as.integer(val); i <- i + 2L }
    else { i <- i + 1L }
  }
  out
}

args <- parse_args(args_raw)
if (is.null(args$screened_hits)) stop("Required: --screened-hits")

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$forest_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
num_col <- function(df, ...) {
  for (nm in list(...)) {
    if (nm %in% names(df)) return(as.numeric(df[[nm]]))
  }
  rep(NA_real_, nrow(df))
}

int_col <- function(df, ...) {
  as.integer(num_col(df, ...) %||% 0)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

classify_pattern <- function(df) {
  eur_sig  <- (int_col(df, "eur_is_gw", "EUR_homogeneous_is_gw") == 1L) |
              (int_col(df, "eur_is_suggestive", "EUR_homogeneous_is_suggestive") == 1L)
  hbcc_sig <- (int_col(df, "hbcc_is_gw", "HBCC_AFR_enriched_is_gw") == 1L) |
              (int_col(df, "hbcc_is_suggestive", "HBCC_AFR_enriched_is_suggestive") == 1L)

  eur_b  <- num_col(df, "eur_meta_beta", "eur_beta", "EUR_homogeneous_beta")
  hbcc_b <- num_col(df, "hbcc_meta_beta", "hbcc_beta", "HBCC_AFR_enriched_beta")

  sign_match <- function(a, b) {
    !is.na(a) & !is.na(b) & a != 0 & b != 0 & sign(a) == sign(b)
  }

  pattern <- dplyr::case_when(
    eur_sig & hbcc_sig & sign_match(eur_b, hbcc_b) ~ "both_sig_same_direction",
    eur_sig & hbcc_sig & !sign_match(eur_b, hbcc_b) ~ "both_sig_opposite_direction",
    eur_sig & !hbcc_sig ~ "eur_only",
    hbcc_sig & !eur_sig ~ "hbcc_only",
    TRUE ~ "neither_stratum_sig"
  )

  eur_mean  <- num_col(df, "eur_mean_beta")
  hbcc_mean <- num_col(df, "hbcc_mean_beta")

  df %>%
    mutate(
      eur_sig  = eur_sig,
      hbcc_sig = hbcc_sig,
      eur_beta_plot  = ifelse(is.na(eur_b), eur_mean, eur_b),
      hbcc_beta_plot = ifelse(is.na(hbcc_b), hbcc_mean, hbcc_b),
      eur_p_plot  = num_col(df, "eur_meta_p", "eur_p", "EUR_homogeneous_p"),
      hbcc_p_plot = num_col(df, "hbcc_meta_p", "hbcc_p", "HBCC_AFR_enriched_p"),
      pooled_p    = as.numeric(pooled_p),
      pooled_beta = as.numeric(pooled_beta),
      pooled_i2   = as.numeric(pooled_i2),
      masking_score = as.integer(masking_score),
      is_priority   = as.integer(is_priority),
      pooled_is_gw  = as.integer(pooled_is_gw),
      sig_level = ifelse(pooled_is_gw == 1L, "Genome-wide", "Suggestive"),
      ancestry_pattern = factor(
        pattern,
        levels = c(
          "both_sig_same_direction",
          "both_sig_opposite_direction",
          "eur_only",
          "hbcc_only",
          "neither_stratum_sig"
        )
      ),
      locus_label = paste0(cell_type, ": ", sub(":.*", "", marker))
    )
}

pattern_colours <- c(
  both_sig_same_direction      = "#1a9850",
  both_sig_opposite_direction    = "#d73027",
  eur_only                       = "#4575b4",
  hbcc_only                      = "#f46d43",
  neither_stratum_sig            = "#bdbdbd"
)

pattern_labels <- c(
  both_sig_same_direction      = "Both strata sig, same direction",
  both_sig_opposite_direction  = "Both strata sig, opposite direction",
  eur_only                     = "EUR only",
  hbcc_only                    = "HBCC only",
  neither_stratum_sig          = "Neither stratum sig"
)

save_plot <- function(p, path, w = 10, h = 7) {
  ggsave(path, plot = p, width = w, height = h, dpi = 200, bg = "white")
  message("Saved: ", path)
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
message("Loading screened hits...")
hits <- read.delim(args$screened_hits, stringsAsFactors = FALSE)
hits <- classify_pattern(hits)
message("  ", nrow(hits), " hits; patterns: ",
        paste(names(table(hits$ancestry_pattern)), table(hits$ancestry_pattern), sep = "=", collapse = ", "))

detail <- NULL
if (!is.null(args$locus_detail) && file.exists(args$locus_detail)) {
  detail <- read.delim(args$locus_detail, stringsAsFactors = FALSE)
  message("  ", nrow(detail), " per-cohort detail rows")
}

# ---------------------------------------------------------------------------
# 1. Overall pattern summary
# ---------------------------------------------------------------------------
pat_counts <- hits %>%
  count(ancestry_pattern, sig_level, name = "n")

p_summary <- ggplot(pat_counts, aes(x = sig_level, y = n, fill = ancestry_pattern)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(n >= 5, n, "")),
            position = position_stack(vjust = 0.5), size = 3.5, colour = "white") +
  scale_fill_manual(values = pattern_colours, labels = pattern_labels, name = "Ancestry pattern") +
  labs(
    title = "Ancestry Tier 1 — EUR vs HBCC stratum concordance",
    subtitle = "Screened hits (GW + suggestive with I\u00b2 \u2265 50% or coloc priority)",
    x = NULL, y = "Number of hits"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9))

save_plot(p_summary, file.path(args$output_dir, "ancestry_pattern_summary.png"), 11, 6)

# ---------------------------------------------------------------------------
# 2. By cell type (GW hits only for readability)
# ---------------------------------------------------------------------------
gw_hits <- hits %>% filter(pooled_is_gw == 1L)
if (nrow(gw_hits) > 0) {
  ct_counts <- gw_hits %>%
    count(cell_type, ancestry_pattern, name = "n")

  p_ct <- ggplot(ct_counts, aes(x = reorder(cell_type, n, sum), y = n, fill = ancestry_pattern)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_fill_manual(values = pattern_colours, labels = pattern_labels, name = "Ancestry pattern") +
    labs(
      title = "Genome-wide hits — ancestry pattern by cell type",
      x = NULL, y = "Count"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))

  save_plot(p_ct, file.path(args$output_dir, "ancestry_pattern_by_celltype_gw.png"), 10, 8)
}

# ---------------------------------------------------------------------------
# 3. EUR vs HBCC meta-beta scatter
# ---------------------------------------------------------------------------
scatter_df <- hits %>%
  filter(!is.na(eur_beta_plot) | !is.na(hbcc_beta_plot)) %>%
  mutate(
    neg_log10_eur  = -log10(pmax(eur_p_plot,  .Machine$double.xmin)),
    neg_log10_hbcc = -log10(pmax(hbcc_p_plot, .Machine$double.xmin)),
    plot_size = ifelse(is_priority == 1L, 3.5, 1.8),
    label_me = masking_score >= 6L & is_priority == 1L
  )

label_df <- scatter_df %>%
  filter(label_me) %>%
  arrange(desc(masking_score), pooled_p) %>%
  slice_head(n = args$max_labels)

p_beta <- ggplot(scatter_df, aes(x = eur_beta_plot, y = hbcc_beta_plot)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(aes(colour = ancestry_pattern, size = plot_size), alpha = 0.75) +
  geom_text(data = label_df, aes(label = cell_type),
            size = 2.8, hjust = -0.1, vjust = 0.5, check_overlap = TRUE) +
  scale_colour_manual(values = pattern_colours, labels = pattern_labels, name = "Pattern") +
  scale_size_identity() +
  labs(
    title = "EUR vs HBCC stratum effect sizes",
    subtitle = "Points on diagonal = consistent direction; off-diagonal quadrants = ancestry-specific or discordant",
    x = "EUR-homogeneous meta beta",
    y = "HBCC AFR-enriched meta beta"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

save_plot(p_beta, file.path(args$output_dir, "eur_hbcc_beta_scatter.png"), 11, 8)

# ---------------------------------------------------------------------------
# 4. -log10(p) scatter
# ---------------------------------------------------------------------------
p_pval <- ggplot(
  scatter_df %>% filter(eur_sig | hbcc_sig),
  aes(x = neg_log10_eur, y = neg_log10_hbcc)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = -log10(5e-8), linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = ancestry_pattern, size = plot_size), alpha = 0.75) +
  geom_text(data = label_df, aes(label = cell_type),
            size = 2.8, hjust = -0.1, vjust = 0.5, check_overlap = TRUE) +
  scale_colour_manual(values = pattern_colours, labels = pattern_labels, name = "Pattern") +
  scale_size_identity() +
  labs(
    title = "EUR vs HBCC stratum significance",
    subtitle = "Dashed lines = genome-wide threshold (p < 5\u00d710\u207b\u2078)",
    x = expression(-log[10](p)[EUR]),
    y = expression(-log[10](p)[HBCC])
  ) +
  theme_classic(base_size = 12)

save_plot(p_pval, file.path(args$output_dir, "eur_hbcc_pval_scatter.png"), 11, 8)

# ---------------------------------------------------------------------------
# 5. Top-hit beta heatmap
# ---------------------------------------------------------------------------
top_hits <- hits %>%
  arrange(desc(masking_score), pooled_p) %>%
  slice_head(n = 40) %>%
  mutate(
    hit_id = paste0(cell_type, "\n", sub("^chr([0-9]+):([0-9]+).*", "chr\\1:\\2", marker)),
    eur_beta_plot  = ifelse(is.na(eur_beta_plot), 0, eur_beta_plot),
    hbcc_beta_plot = ifelse(is.na(hbcc_beta_plot), 0, hbcc_beta_plot)
  ) %>%
  select(hit_id, ancestry_pattern, pooled_beta, eur_beta_plot, hbcc_beta_plot) %>%
  pivot_longer(c(pooled_beta, eur_beta_plot, hbcc_beta_plot),
               names_to = "stratum", values_to = "beta") %>%
  mutate(
    stratum = factor(stratum,
      levels = c("pooled_beta", "eur_beta_plot", "hbcc_beta_plot"),
      labels = c("Pooled", "EUR meta", "HBCC meta")
    ),
    beta_scaled = ifelse(abs(beta) < 1e-6, 0, sign(beta) * log1p(abs(beta) * 100))
  )

p_heat <- ggplot(top_hits, aes(x = stratum, y = hit_id, fill = beta_scaled)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    name = "Signed\nlog beta"
  ) +
  labs(
    title = "Top 40 ancestry-screened hits — effect direction by stratum",
    subtitle = "Green pattern = both strata same direction; blue = EUR only",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

save_plot(p_heat, file.path(args$output_dir, "top_hits_beta_heatmap.png"), 8, 12)

# ---------------------------------------------------------------------------
# 5b. Dumbbell plot — top priority loci (pooled / EUR / HBCC betas)
# ---------------------------------------------------------------------------
top_dumb <- hits %>%
  filter(is_priority == 1L | masking_score >= 8L) %>%
  arrange(desc(masking_score), pooled_p) %>%
  distinct(cell_type, marker, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    hit_label = paste0(cell_type, " | ", sub("^chr([0-9]+):([0-9]+).*", "chr\\1:\\2", marker)),
    hit_label = factor(hit_label, levels = rev(unique(hit_label)))
  ) %>%
  pivot_longer(
    c(pooled_beta, eur_beta_plot, hbcc_beta_plot),
    names_to = "stratum", values_to = "beta"
  ) %>%
  mutate(stratum = factor(stratum,
    levels = c("pooled_beta", "eur_beta_plot", "hbcc_beta_plot"),
    labels = c("Pooled", "EUR", "HBCC")
  ))

p_dumb <- ggplot(top_dumb, aes(x = beta, y = hit_label, colour = stratum)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(aes(group = hit_label), colour = "grey75", linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = c(Pooled = "#333333", EUR = "#2166AC", HBCC = "#D95F02"), name = "Stratum") +
  facet_wrap(~ ancestry_pattern, scales = "free_y", ncol = 1) +
  labs(
    title = "Top priority loci — beta across strata",
    subtitle = "Lines connect pooled, EUR, and HBCC effects; gaps = stratum not significant",
    x = "Effect size (beta)", y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

save_plot(p_dumb, file.path(args$output_dir, "top_priority_dumbbell.png"), 10, 10)

# ---------------------------------------------------------------------------
# 6. Consistent vs ancestry-specific tables (TSV for report)
# ---------------------------------------------------------------------------
consistent <- hits %>%
  filter(ancestry_pattern == "both_sig_same_direction") %>%
  arrange(pooled_p) %>%
  select(cell_type, marker, pooled_p, pooled_beta, pooled_i2,
         eur_beta_plot, hbcc_beta_plot, masking_score, is_priority)

ancestry_specific <- hits %>%
  filter(ancestry_pattern %in% c("eur_only", "hbcc_only", "both_sig_opposite_direction")) %>%
  arrange(desc(masking_score), pooled_p) %>%
  select(cell_type, marker, ancestry_pattern, pooled_p, pooled_beta, pooled_i2,
         eur_beta_plot, hbcc_beta_plot, masking_score, is_priority)

write.table(consistent, file.path(args$output_dir, "hits_consistent_across_ancestries.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ancestry_specific, file.path(args$output_dir, "hits_ancestry_specific_or_discordant.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote summary TSVs: hits_consistent_across_ancestries.tsv (n=", nrow(consistent),
        "), hits_ancestry_specific_or_discordant.tsv (n=", nrow(ancestry_specific), ")")

# ---------------------------------------------------------------------------
# 7. Forest plots from per-cohort detail
# ---------------------------------------------------------------------------
if (!is.null(detail) && nrow(detail) > 0) {
  cohort_order <- c(
    "ROSMAP", "ROSMAP_array", "Mayo", "MSBB", "NABEC",
    "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GTEx_v10", "GVEX",
    "AMP_AD_Rush", "AMP_AD_Mayo",
    "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650"
  )

  ancestry_colours <- c(EUR = "#2166AC", AFR_enriched = "#D95F02", mixed = "#7570b3")

  forest_targets <- hits %>%
    filter(masking_score >= 6L | is_priority == 1L) %>%
    arrange(desc(masking_score), pooled_p) %>%
    distinct(cell_type, marker, .keep_all = TRUE) %>%
    slice_head(n = args$max_forest_plots)

  for (i in seq_len(nrow(forest_targets))) {
    ct  <- forest_targets$cell_type[i]
    mkr <- forest_targets$marker[i]
    pat <- as.character(forest_targets$ancestry_pattern[i])

    sub <- detail %>%
      filter(cell_type == ct, marker == mkr) %>%
      mutate(
        beta = as.numeric(beta),
        se   = as.numeric(se),
        ci_lo = beta - 1.96 * se,
        ci_hi = beta + 1.96 * se,
        cohort = factor(cohort, levels = rev(cohort_order[cohort_order %in% cohort]))
      )

    if (nrow(sub) < 2L) next

    safe_name <- gsub("[^A-Za-z0-9._-]", "_", paste0(ct, "_", sub(":.*", "", mkr)))
    subtitle <- paste0(
      pat, " | pooled p=", formatC(forest_targets$pooled_p[i], format = "e", digits = 2),
      " | I\u00b2=", round(forest_targets$pooled_i2[i], 1), "%"
    )

    p_forest <- ggplot(sub, aes(x = beta, y = cohort, colour = ancestry_group)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.8) +
      geom_point(size = 3) +
      scale_colour_manual(values = ancestry_colours, name = "Ancestry") +
      labs(
        title = paste0(ct, " - ", mkr),
        subtitle = subtitle,
        x = "Per-cohort beta (95% CI)", y = NULL
      ) +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom")

    out_pdf <- file.path(args$forest_dir, paste0(safe_name, "_forest.pdf"))
    ggsave(out_pdf, plot = p_forest, width = 9, height = max(4, nrow(sub) * 0.35 + 2), device = "pdf")
  }
  message("Forest plots written to: ", args$forest_dir)
}

message("Done.")
