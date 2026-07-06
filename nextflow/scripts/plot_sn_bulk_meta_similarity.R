#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--outdir", type = "character",
              default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity")
)
opt <- parse_args(OptionParser(option_list = option_list))

plot_dir <- file.path(opt$outdir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

theme_meta <- theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  fread(path)
}

snp <- read_if_exists(file.path(opt$outdir, "snp_similarity_summary.tsv"))
if (!is.null(snp) && nrow(snp) > 0) {
  p_z <- ggplot(snp, aes(x = reorder(sn_cell_type, z_cor_pearson), y = z_cor_pearson)) +
    geom_col(fill = "#2171b5") +
    coord_flip() +
    labs(title = "Genome-wide Z-score correlation", x = NULL, y = "Pearson r")
  ggsave(file.path(plot_dir, "z_correlation_by_cell_type.png"), p_z,
         width = 7, height = 5, dpi = 300)

  p_sign <- ggplot(snp, aes(x = reorder(sn_cell_type, sign_concordance_all), y = sign_concordance_all)) +
    geom_col(fill = "#238b45") +
    coord_flip() +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    labs(title = "Genome-wide effect direction concordance", x = NULL, y = "Fraction same sign")
  ggsave(file.path(plot_dir, "direction_concordance_by_cell_type.png"), p_sign,
         width = 7, height = 5, dpi = 300)
}

rg <- read_if_exists(file.path(opt$outdir, "ldsc", "rg_results", "ldsc_rg_summary.tsv"))
if (!is.null(rg) && nrow(rg) > 0) {
  rg[is.na(reliability_flags) | reliability_flags == "", reliability_flags := "(none)"]
  rg[, `:=`(
    rg_num = suppressWarnings(as.numeric(rg)),
    rg_se_num = suppressWarnings(as.numeric(rg_se)),
    status_group = ifelse(status == "ok" & is.finite(suppressWarnings(as.numeric(rg))),
                          "numeric rg", status)
  )]
  rg[, has_numeric_rg := is.finite(rg_num) & is.finite(rg_se_num)]
  rg[, y_fail_marker := 0]
  all_flags <- sort(unique(rg$reliability_flags))
  rg[, reliability_flags := factor(reliability_flags, levels = all_flags)]

  p_rg <- ggplot(rg, aes(x = sn_cell_type, colour = reliability_flags)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_pointrange(
      data = rg[has_numeric_rg == TRUE],
      aes(y = rg_num, ymin = rg_num - rg_se_num, ymax = rg_num + rg_se_num),
      na.rm = TRUE
    ) +
    geom_point(
      data = rg[has_numeric_rg == FALSE],
      aes(y = y_fail_marker),
      shape = 4, size = 2.8, stroke = 1.0, na.rm = TRUE
    ) +
    labs(title = "LDSC genetic correlation: sn vs bulk",
         subtitle = "Error bars are +/- 1 SE; X marker at 0 indicates NA/failed rg",
         x = NULL, y = "rg", colour = "Reliability flags") +
    scale_colour_discrete(drop = FALSE) +
    theme_meta
  ggsave(file.path(plot_dir, "ldsc_rg_by_cell_type.png"), p_rg,
         width = 11, height = 5, dpi = 300)

  p_status <- ggplot(rg, aes(x = sn_cell_type, fill = status_group)) +
    geom_bar(width = 0.75) +
    labs(title = "LDSC rg status by cell type", x = NULL, y = "Count", fill = "Status") +
    theme_meta
  ggsave(file.path(plot_dir, "ldsc_rg_status.png"), p_status,
         width = 9, height = 4.5, dpi = 300)
}

rep <- read_if_exists(file.path(opt$outdir, "replication_enrichment_summary.tsv"))
if (!is.null(rep) && nrow(rep) > 0) {
  rep_plot <- rep[selection %in% c("p<1e-05", "top1000")]
  if (nrow(rep_plot) > 0) {
    p_rep <- ggplot(rep_plot,
                    aes(x = sn_cell_type, y = target_nominal_fraction, fill = direction)) +
      geom_col(position = position_dodge(0.75), width = 0.7) +
      facet_wrap(~selection, ncol = 1) +
      labs(title = "Replication-style nominal enrichment", x = NULL,
           y = "Fraction target P < 0.05", fill = NULL) +
      theme_meta
    ggsave(file.path(plot_dir, "replication_nominal_enrichment.png"), p_rep,
           width = 11, height = 7, dpi = 300)
  }
}

locus <- read_if_exists(file.path(opt$outdir, "lead_locus_concordance_summary.tsv"))
if (!is.null(locus) && nrow(locus) > 0) {
  p_locus <- ggplot(locus, aes(x = sn_cell_type, y = n_region_pass, fill = direction)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    facet_wrap(~threshold, ncol = 1, scales = "free_y") +
    labs(title = "Lead loci with target signal in the same region", x = NULL,
         y = "Lead loci with target region P below threshold", fill = NULL) +
    theme_meta
  ggsave(file.path(plot_dir, "lead_locus_region_concordance.png"), p_locus,
         width = 11, height = 7, dpi = 300)
}

message("Plots written to: ", plot_dir)
