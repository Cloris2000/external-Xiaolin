#!/usr/bin/env Rscript
# Reads all per-cell-type concordance TSVs and generates a stacked bar plot
# showing # of independent loci per cell type classified as:
#   same_direction / opposite_direction / not_found
# Bars are split by GW-significant (P < 5e-8) vs suggestive (P < 1e-5).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

CONC_DIR <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts/concordance"
OUT_DIR  <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts/plots"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------- load all TSVs ----------------------------------------------
tsv_files <- list.files(CONC_DIR, pattern = "_concordance\\.tsv$", full.names = TRUE)
if (length(tsv_files) == 0) stop("No concordance TSVs found in: ", CONC_DIR)
message("Loading ", length(tsv_files), " cell types")

dat <- rbindlist(lapply(tsv_files, fread), use.names = TRUE, fill = TRUE)

# ---------- summarise ---------------------------------------------------
# Create a combined group: sig level + category
dat[, group := paste0(bulk_sig, "__", category)]

# Count per cell type + group
counts <- dat[, .N, by = .(celltype, bulk_sig, category)]

# ---------- ordering / factor levels -----------------------------------
ct_order <- c(
  "Oligodendrocyte", "OPC", "Astrocyte", "Microglia", "Endothelial",
  "Pericyte", "VLMC",
  "IT", "L4.IT", "L5.6.IT.Car3", "L5.ET", "L5.6.NP", "L6.CT", "L6b",
  "LAMP5", "PAX6", "PVALB", "SST", "VIP"
)
# only keep present cell types, in order
ct_order <- ct_order[ct_order %in% unique(counts$celltype)]
counts$celltype <- factor(counts$celltype, levels = ct_order)

# category × sig level ordering & colours
counts[, fill_group := factor(
  paste0(category, "__", bulk_sig),
  levels = c(
    "same_direction__GW_significant",
    "same_direction__suggestive",
    "opposite_direction__GW_significant",
    "opposite_direction__suggestive",
    "not_found__GW_significant",
    "not_found__suggestive"
  )
)]

fill_colours <- c(
  "same_direction__GW_significant"    = "#1a6b3c",
  "same_direction__suggestive"        = "#74c476",
  "opposite_direction__GW_significant"= "#a50f15",
  "opposite_direction__suggestive"    = "#fb6a4a",
  "not_found__GW_significant"         = "#4d4d4d",
  "not_found__suggestive"             = "#bdbdbd"
)

fill_labels <- c(
  "same_direction__GW_significant"    = "Same direction (GW-sig)",
  "same_direction__suggestive"        = "Same direction (suggestive)",
  "opposite_direction__GW_significant"= "Opposite direction (GW-sig)",
  "opposite_direction__suggestive"    = "Opposite direction (suggestive)",
  "not_found__GW_significant"         = "No sn proxy within 500kb (GW-sig)",
  "not_found__suggestive"             = "No sn proxy within 500kb (suggestive)"
)

# ---------- plot -------------------------------------------------------
p <- ggplot(counts, aes(x = celltype, y = N, fill = fill_group)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = fill_colours, labels = fill_labels,
                    name = "Concordance with sn GWAS\n(direction only; no sn P threshold)") +
  labs(
    title    = "Bulk vs snRNA-seq GWAS concordance per cell type",
    subtitle = "Independent loci (P < 1e-5) from 15-cohort bulk meta · direction compared to 3-cohort sn meta (±500 kb proxy)",
    x        = NULL,
    y        = "Number of independent loci"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    legend.position  = "right",
    legend.text      = element_text(size = 11),
    legend.title     = element_text(size = 12, face = "bold"),
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(size = 12, colour = "grey40"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(reverse = FALSE))

ggsave(file.path(OUT_DIR, "bulk_sn_concordance_barplot.png"),
       plot = p, width = 14, height = 6, dpi = 200, bg = "white")
message("Saved: ", file.path(OUT_DIR, "bulk_sn_concordance_barplot.png"))

# ---------- also save summary table ------------------------------------
fwrite(counts[order(celltype, fill_group)],
       file.path(CONC_DIR, "concordance_summary.tsv"), sep = "\t")
message("Summary table: ", file.path(CONC_DIR, "concordance_summary.tsv"))
