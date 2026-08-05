#!/usr/bin/env Rscript
# Regenerate the snRNA concordance bar chart (Panel A of Figure 5)
# with colors that do NOT clash with the B/C/D ancestry palette
# (blue=#2166AC / orange=#D95F02 / red=#B2182B).
#
# New palette (ColorBrewer Dark2 / colorblind-safe):
#   same direction    → teal-green  #1B9E77
#   opposite direction → purple     #7570B3
#   not found in sn   → light grey  #D3D3D3

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

TSV <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity_design_matrix/top_hits/cell_type_concordance_summary.tsv"
OUT <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity_design_matrix/top_hits/bulk_top_hits_direction_bar.png"

# Nature NPG-inspired palette (colorblind-safe, widely used in Nature journals)
clr_same <- "#3C5488"   # deep navy-blue  (same direction)
clr_opp  <- "#E64B35"   # warm orange-red (opposite direction)
clr_miss <- "#C8C8C8"   # light grey      (not found in sn)

summary_dt <- fread(TSV)

# Melt to long format for stacked bar
bar_dt <- melt(
  summary_dt[, .(bulk_cell_type,
                 `same direction`      = n_concordant,
                 `opposite direction`  = n_discordant,
                 `not found in sn`     = n_not_in_sn)],
  id.vars       = "bulk_cell_type",
  variable.name = "category",
  value.name    = "count"
)
bar_dt[, category := factor(category,
                             levels = c("same direction",
                                        "opposite direction",
                                        "not found in sn"))]
bar_dt[, bulk_cell_type := factor(
  bulk_cell_type,
  levels = summary_dt[order(-pct_concordant)]$bulk_cell_type
)]

# Label: % concordant inside bar, total n above bar
# Short label: just the percentage — avoids overflow/overlap inside the bar
summary_dt[, pct_label := ifelse(n_concordant >= 1,
                                  paste0(round(n_concordant / n_found_sn * 100, 0), "%"),
                                  "")]
summary_dt[, label_y   := n_concordant / 2]
summary_dt[, bulk_cell_type := factor(
  bulk_cell_type,
  levels = summary_dt[order(-pct_concordant)]$bulk_cell_type
)]

p <- ggplot(bar_dt, aes(x = bulk_cell_type, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    data    = summary_dt[pct_label != "" & n_concordant >= 2],
    mapping = aes(x = bulk_cell_type, y = label_y,
                  label = pct_label, fill = NULL),
    colour = "white", fontface = "bold", size = 3.8, vjust = 0.5
  ) +
  scale_fill_manual(
    values = c("same direction"     = clr_same,
               "opposite direction" = clr_opp,
               "not found in sn"    = clr_miss),
    breaks = c("same direction", "opposite direction", "not found in sn")
  ) +
  labs(
    title    = "Direction concordance of bulk independent suggestive loci in snRNA-seq GWAS",
    subtitle = "Bulk hits p < 1e-05, clumped (\u00b1500 kb) | x-axis: cell type, % concordant (concordant / found in sn)",
    x    = NULL,
    y    = "# independent loci",
    fill = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y      = element_text(size = 12),
    axis.title.y     = element_text(size = 13),
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 11, colour = "grey40"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(1.0, "lines"),
    panel.grid.minor = element_blank()
  )

ggsave(OUT, p, width = 13, height = 5.5, dpi = 200, bg = "white")
message("Saved: ", OUT)
