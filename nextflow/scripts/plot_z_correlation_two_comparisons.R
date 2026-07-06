#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--base_dir", type = "character",
              default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity",
              help = "Directory for sn vs bulk comparison"),
  make_option("--design_dir", type = "character",
              default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity_design_matrix",
              help = "Directory for sn vs bulk (design matrix) comparison"),
  make_option("--out_png", type = "character",
              default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity_design_matrix/plots/z_correlation_combined_two_bulk_models.png",
              help = "Output PNG path")
)
opt <- parse_args(OptionParser(option_list = option_list))

base_tsv <- file.path(opt$base_dir, "snp_similarity_summary.tsv")
design_tsv <- file.path(opt$design_dir, "snp_similarity_summary.tsv")

if (!file.exists(base_tsv)) stop("Missing file: ", base_tsv)
if (!file.exists(design_tsv)) stop("Missing file: ", design_tsv)

base <- fread(base_tsv)[, .(
  sn_cell_type,
  z_cor_pearson = as.numeric(z_cor_pearson),
  model = "bulk meta"
)]
design <- fread(design_tsv)[, .(
  sn_cell_type,
  z_cor_pearson = as.numeric(z_cor_pearson),
  model = "bulk meta (design matrix)"
)]

both <- rbindlist(list(base, design), fill = TRUE)
setorder(base, -z_cor_pearson)
ct_levels <- base$sn_cell_type
both[, sn_cell_type := factor(sn_cell_type, levels = ct_levels)]

dir.create(dirname(opt$out_png), recursive = TRUE, showWarnings = FALSE)

p <- ggplot(both, aes(x = sn_cell_type, y = z_cor_pearson, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = c("bulk meta" = "#2171b5", "bulk meta (design matrix)" = "#6baed6")) +
  labs(
    title = "Z-score correlation by cell type: two bulk meta models",
    subtitle = "snRNAseq meta vs bulk meta (standard) and bulk meta (design matrix)",
    x = NULL,
    y = "Pearson correlation of SNP Z-scores",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(opt$out_png, p, width = 11, height = 5.5, dpi = 300)
message("Wrote: ", opt$out_png)
