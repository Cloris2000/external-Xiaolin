#!/usr/bin/env Rscript
# Ancestry forest plots for pooled 15-cohort lead SNPs.
#
# Input: results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv
#   (from scripts/extract_ancestry_lead_effects.py)
#
# For each lead: rows = EUR meta / AFR meta / Mayo AMR / Pooled 15-cohort meta.
# Style aligned with scripts/plot_forest_vip_tmem106b.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(optparse)
})

option_list <- list(
  make_option(c("--effects"), type = "character",
              default = "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"),
  make_option(c("--output_dir"), type = "character",
              default = "results/meta_sensitivity/ancestry_lead_effects/forests"),
  make_option(c("--cell_type"), type = "character", default = NULL,
              help = "Optional single cell type"),
  make_option(c("--max_leads_per_ct"), type = "integer", default = 6,
              help = "Max lead_rank to plot per cell type (default 6)"),
  make_option(c("--combined"), type = "logical", default = TRUE,
              help = "Also write one multi-panel PNG per cell type")
)
opt <- parse_args(OptionParser(option_list = option_list))

effects_path <- opt$effects
out_dir <- opt$output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Alleles T/F must stay character (readr would coerce "T" -> TRUE)
df <- read_tsv(
  effects_path,
  show_col_types = FALSE,
  col_types = cols(
    effect_allele = col_character(),
    other_allele = col_character(),
    stratum = col_character(),
    stratum_label = col_character(),
    marker = col_character(),
    cell_type = col_character(),
    harmonize_status = col_character()
  )
)
if (!is.null(opt$cell_type)) {
  df <- df %>% filter(cell_type == opt$cell_type)
}
df <- df %>% filter(lead_rank <= opt$max_leads_per_ct)

if (nrow(df) == 0) stop("No rows to plot in: ", effects_path)

stratum_levels <- c(
  "EUR meta (White)",
  "AFR meta (African American)",
  "Mayo AMR (Latino; single cohort)",
  "Pooled 15-cohort meta"
)
# Fall back to stratum codes if labels differ
if (!all(df$stratum_label %in% stratum_levels)) {
  stratum_levels <- c("EUR", "AFR", "AMR", "Pooled")
  df <- df %>% mutate(stratum_label = stratum)
}

clr <- c(
  "EUR meta (White)" = "#2166AC",
  "AFR meta (African American)" = "#D95F02",
  "Mayo AMR (Latino; single cohort)" = "#1B9E77",
  "Pooled 15-cohort meta" = "#B2182B",
  "EUR" = "#2166AC",
  "AFR" = "#D95F02",
  "AMR" = "#1B9E77",
  "Pooled" = "#B2182B"
)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                formatC(p, format = "f", digits = 3)))
}

plot_one_lead <- function(d) {
  d <- d %>%
    mutate(
      stratum_label = factor(stratum_label, levels = rev(stratum_levels)),
      is_meta = stratum %in% c("Pooled"),
      has_est = !is.na(beta) & !is.na(se),
      plabel = paste0("P = ", fmt_p(p)),
      nlabel = ifelse(has_est & !is.na(n),
                      paste0("N≈", formatC(as.integer(n), big.mark = ",")),
                      ""),
      # placeholder at 0 for missing so the row still appears
      beta_plot = ifelse(has_est, beta, 0),
      ci_lo_plot = ifelse(has_est, ci_lo, NA_real_),
      ci_hi_plot = ifelse(has_est, ci_hi, NA_real_)
    )

  ea <- unique(d$effect_allele)[1]
  marker <- unique(d$marker)[1]
  ct <- unique(d$cell_type)[1]
  title <- sprintf("%s — %s (effect allele %s)", ct, marker, ea)

  xmin <- min(d$ci_lo, na.rm = TRUE)
  xmax <- max(d$ci_hi, na.rm = TRUE)
  if (!is.finite(xmin) || !is.finite(xmax)) {
    xmin <- -0.2; xmax <- 0.2
  }
  xpad <- max(0.05, 0.15 * (xmax - xmin))

  ggplot(d, aes(y = stratum_label, x = beta_plot, colour = stratum_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_errorbarh(
      data = d %>% filter(has_est),
      aes(xmin = ci_lo, xmax = ci_hi, linewidth = is_meta),
      height = 0.25
    ) +
    geom_point(
      data = d %>% filter(has_est),
      aes(shape = is_meta, size = is_meta)
    ) +
    geom_point(
      data = d %>% filter(!has_est),
      aes(x = 0),
      shape = 4, size = 2.5, colour = "grey55"
    ) +
    geom_text(
      aes(x = xmin - xpad * 0.15, label = nlabel),
      hjust = 1, size = 2.8, colour = "grey35"
    ) +
    geom_text(
      data = d %>% filter(has_est),
      aes(x = xmax + xpad * 0.15, label = plabel),
      hjust = 0, size = 2.8, colour = "grey20"
    ) +
    geom_text(
      data = d %>% filter(!has_est),
      aes(x = xmax + xpad * 0.15, label = "not in sumstats"),
      hjust = 0, size = 2.6, colour = "grey50", fontface = "italic"
    ) +
    scale_colour_manual(values = clr, guide = "none") +
    scale_linewidth_manual(values = c("FALSE" = 0.7, "TRUE" = 1.2), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2.4, "TRUE" = 3.6), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.35, 0.45))) +
    labs(title = title, x = "Beta (95% CI)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      axis.text.y = element_text(size = 9, colour = "black"),
      axis.text.x = element_text(size = 8, colour = "black"),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      plot.margin = margin(6, 90, 6, 6)
    )
}

# Per-lead PNGs
leads <- df %>% distinct(cell_type, lead_rank, marker)
message("Plotting ", nrow(leads), " lead forests -> ", out_dir)

for (i in seq_len(nrow(leads))) {
  ct <- leads$cell_type[i]
  rk <- leads$lead_rank[i]
  mk <- leads$marker[i]
  d <- df %>% filter(cell_type == ct, lead_rank == rk)
  p <- plot_one_lead(d)
  safe_mk <- gsub("[:/]", "_", mk)
  out <- file.path(out_dir, sprintf("%s_lead%02d_%s_forest.png", ct, rk, safe_mk))
  ggsave(out, plot = p, width = 7.2, height = 3.2, dpi = 200, bg = "white")
}

# Combined small-multiples per cell type
if (isTRUE(opt$combined)) {
  for (ct in unique(df$cell_type)) {
    d_ct <- df %>% filter(cell_type == ct)
    n_leads <- n_distinct(d_ct$lead_rank)
    if (n_leads == 0) next

    d_ct <- d_ct %>%
      mutate(
        facet_lab = paste0("lead ", lead_rank, "\n", marker),
        stratum_label = factor(stratum_label, levels = rev(stratum_levels)),
        is_meta = stratum %in% c("Pooled"),
        has_est = !is.na(beta) & !is.na(se)
      )

    p_all <- ggplot(d_ct %>% filter(has_est),
                    aes(y = stratum_label, x = beta, colour = stratum_label)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, linewidth = is_meta), height = 0.22) +
      geom_point(aes(shape = is_meta, size = is_meta)) +
      facet_wrap(~ facet_lab, scales = "free_x", ncol = min(3L, n_leads)) +
      scale_colour_manual(values = clr, guide = "none") +
      scale_linewidth_manual(values = c("FALSE" = 0.6, "TRUE" = 1.1), guide = "none") +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
      scale_size_manual(values = c("FALSE" = 2.0, "TRUE" = 3.2), guide = "none") +
      labs(
        title = paste0(ct, ": ancestry effects at pooled lead SNPs"),
        x = "Beta (95% CI), harmonized to pooled effect allele",
        y = NULL
      ) +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        strip.text = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.5),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.25)
      )

    ncol <- min(3L, n_leads)
    nrow <- ceiling(n_leads / ncol)
    out <- file.path(out_dir, sprintf("%s_ancestry_lead_forests.png", ct))
    ggsave(out, plot = p_all,
           width = 3.6 * ncol + 0.5,
           height = 2.4 * nrow + 0.6,
           dpi = 200, bg = "white", limitsize = FALSE)
    message("  combined: ", out)
  }
}

message("Done.")
