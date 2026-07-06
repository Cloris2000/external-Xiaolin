#!/usr/bin/env Rscript
# Aggregate coloc results across all cell types and diseases.
# Generates:
#   results/coloc/coloc_all_results.tsv        -- combined results table
#   results/coloc/plots/coloc_heatmap.png      -- PP.H4 heatmap (cell × disease)
#   results/coloc/plots/coloc_dotplot.png      -- dot plot of top hits
#   results/coloc/coloc_significant.tsv        -- PP.H4 >= 0.5 hits
#
# Usage:
#   Rscript 04_summarize_coloc.R \
#     --results_dir results/coloc/coloc_results \
#     --loci_dir    results/coloc/loci \
#     --output_dir  results/coloc

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--results_dir", type = "character", default = "results/coloc/coloc_results"),
  make_option("--loci_dir",    type = "character", default = "results/coloc/loci"),
  make_option("--output_dir",  type = "character", default = "results/coloc"),
  make_option("--pp_h4_sig",   type = "double",    default = 0.5,
              help = "PP.H4 threshold for 'significant' colocalization [default: 0.5]"),
  make_option("--method",      type = "character", default = "coloc.abf",
              help = "Method to use for heatmap: coloc.abf or coloc.susie [default: coloc.abf]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(file.path(opt$output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load all per-cell-type results
# ---------------------------------------------------------------------------
result_files <- list.files(opt$results_dir, pattern = "_coloc_results\\.tsv$",
                            full.names = TRUE)
if (length(result_files) == 0) {
  stop(paste("No *_coloc_results.tsv files found in", opt$results_dir))
}
cat(sprintf("Loading results from %d cell types...\n", length(result_files)))

all_res <- rbindlist(lapply(result_files, fread), fill = TRUE) %>%
  as.data.frame() %>%
  mutate(PP.H4 = as.numeric(PP.H4))

cat(sprintf("Total coloc pairs: %d\n", nrow(all_res)))

# Write combined table
combined_out <- file.path(opt$output_dir, "coloc_all_results.tsv")
fwrite(all_res, combined_out, sep = "\t", quote = FALSE)
cat(sprintf("Combined results written to %s\n", combined_out))

# Write significant hits
sig_res <- all_res %>% filter(PP.H4 >= opt$pp_h4_sig) %>% arrange(desc(PP.H4))
sig_out <- file.path(opt$output_dir, "coloc_significant.tsv")
fwrite(sig_res, sig_out, sep = "\t", quote = FALSE)
cat(sprintf("Significant hits (PP.H4 >= %.1f): %d written to %s\n",
            opt$pp_h4_sig, nrow(sig_res), sig_out))

# ---------------------------------------------------------------------------
# Heatmap: maximum PP.H4 per (cell_type × disease), coloc.abf only
# ---------------------------------------------------------------------------
heat_data <- all_res %>%
  filter(coloc_method == opt$method) %>%
  group_by(cell_type, disease) %>%
  summarise(max_PP_H4 = max(PP.H4, na.rm = TRUE), .groups = "drop")

# Cell type ordering (broad class grouping)
ct_order <- c(
  "Astrocyte", "Microglia", "Oligodendrocyte", "OPC",
  "Endothelial", "Pericyte", "VLMC",
  "IT", "L4.IT", "L5.ET", "L5.6.IT.Car3", "L5.6.NP", "L6.CT", "L6b",
  "LAMP5", "PAX6", "PVALB", "SST", "VIP"
)
ct_present <- intersect(ct_order, unique(heat_data$cell_type))
ct_missing <- setdiff(unique(heat_data$cell_type), ct_order)
ct_levels  <- c(ct_present, sort(ct_missing))

# Disease ordering (rough neuro grouping)
dis_order <- c(
  "AD_Kunkle2019", "LBD_Chia2021", "FTD_vanderLee2019",
  "PD_Nalls2019", "ALS_vanRheenen2021",
  "SCZ_Trubetskoy2022", "BD_Mullins2021", "MDD_Howard2019"
)
dis_present <- intersect(dis_order, unique(heat_data$disease))
dis_missing <- setdiff(unique(heat_data$disease), dis_order)
dis_levels  <- c(dis_present, sort(dis_missing))

heat_data <- heat_data %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(ct_levels)),
    disease   = factor(disease,   levels = dis_levels),
    # Cap at 1 for display
    pp_h4_plot = pmin(max_PP_H4, 1.0)
  )

p_heat <- ggplot(heat_data, aes(x = disease, y = cell_type, fill = pp_h4_plot)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pp_h4_plot >= 0.5,
                               sprintf("%.2f", pp_h4_plot), "")),
            size = 3, color = "black") +
  scale_fill_gradientn(
    name   = "Max PP.H4",
    colors = c("#f7f7f7", "#fee8c8", "#fdbb84", "#e34a33", "#8b0000"),
    values = scales::rescale(c(0, 0.2, 0.5, 0.8, 1.0)),
    limits = c(0, 1), na.value = "grey90"
  ) +
  labs(
    title    = "Colocalization: Cell-Type GWAS × Disease GWAS",
    subtitle = sprintf("Max PP.H4 per (cell type, disease) pair — %s", opt$method),
    x        = "Disease GWAS",
    y        = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    panel.grid       = element_blank()
  )

heat_out <- file.path(opt$output_dir, "plots", "coloc_heatmap.png")
ggsave(heat_out, p_heat, width = max(8, length(dis_levels) * 1.1),
       height = max(6, length(ct_levels) * 0.45), dpi = 150)
cat(sprintf("Heatmap written to %s\n", heat_out))

# ---------------------------------------------------------------------------
# Dot plot: top hits (PP.H4 >= 0.5)
# ---------------------------------------------------------------------------
if (nrow(sig_res) > 0) {
  dot_data <- sig_res %>%
    filter(coloc_method == opt$method) %>%
    mutate(
      locus_label = sub(paste0("^", cell_type, "_"), "", locus_id),
      label       = paste0(cell_type, " — ", locus_label)
    ) %>%
    arrange(desc(PP.H4)) %>%
    slice_head(n = 60)

  p_dot <- ggplot(dot_data,
                  aes(x = disease, y = reorder(label, PP.H4),
                      size = PP.H4, color = PP.H4)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(2, 8), name = "PP.H4") +
    scale_color_gradientn(
      colors = c("#fee8c8", "#e34a33", "#8b0000"),
      limits = c(0.5, 1), name = "PP.H4"
    ) +
    labs(
      title = "Top Colocalization Hits (PP.H4 ≥ 0.5)",
      x = "Disease", y = "Cell type — locus"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x    = element_text(angle = 35, hjust = 1),
      legend.position = "right"
    )

  dot_out <- file.path(opt$output_dir, "plots", "coloc_dotplot.png")
  ggsave(dot_out, p_dot,
         width  = max(8, length(unique(dot_data$disease)) * 1.2),
         height = max(6, nrow(dot_data) * 0.28 + 2), dpi = 150)
  cat(sprintf("Dot plot written to %s\n", dot_out))
}

# ---------------------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------------------
cat("\n=== Colocalization Summary (PP.H4 >= 0.5) ===\n")
cat(sprintf("%-30s  %s\n", "Locus", "PP.H4"))
if (nrow(sig_res) > 0) {
  top <- sig_res %>%
    filter(coloc_method == opt$method) %>%
    arrange(desc(PP.H4)) %>%
    slice_head(n = 20)
  for (i in seq_len(nrow(top))) {
    cat(sprintf("  %-42s  %-30s  %.3f\n",
                top$locus_id[i], top$disease[i], top$PP.H4[i]))
  }
} else {
  cat("  No loci with PP.H4 >= 0.5\n")
}

cat("\nDone.\n")
