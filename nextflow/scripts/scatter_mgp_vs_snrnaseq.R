#!/usr/bin/env Rscript
# scatter_mgp_vs_snrnaseq.R
#
# Generates scatter plots comparing bulk-derived MGP cell type estimates (y-axis)
# against snRNA-seq direct-count proportions (x-axis) for ROSMAP DLPFC.
#
# Outputs:
#   scatter_mgp_vs_snrnaseq_all_celltypes.png  — multi-panel grid (all mapped cell types)
#   scatter_mgp_vs_snrnaseq_<celltype>.png     — individual plots (one per cell type)
#   scatter_mgp_vs_snrnaseq_correlations.tsv   — Pearson r, p-value, n per cell type

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# ── Paths ────────────────────────────────────────────────────────────────────
mgp_file  <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/cell_proportions.csv"
sn_file   <- "/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv"
meta_file <- paste0(
  "/external/rprshnas01/external_data/rosmap/gene_expression/",
  "RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Metadata/",
  "RNAseq_Harmonization_ROSMAP_combined_metadata.csv"
)
# Provenance files: map synapse file ID → ROSMAP specimenID
prov_dir <- paste0(
  "/external/rprshnas01/external_data/rosmap/gene_expression/",
  "RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/"
)
prov_files <- c(
  file.path(prov_dir, "Rosmap_Batch1_Stranded/ROSMAP_batch1_provenance.csv"),
  file.path(prov_dir, "Rosmap_Batch2_Stranded/ROSMAP_batch2_provenance.csv"),
  file.path(prov_dir, "Rosmap_Batch3_Stranded/ROSMAP_batch3_provenance.csv"),
  file.path(prov_dir, "Rosmap_Batch4_Stranded/ROSMAP_batch4_provenance.csv")
)
outdir <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/scatter_mgp_vs_snrnaseq"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Cell-type mapping: pipeline MGP name → snRNA-seq column name ─────────────
# snRNA column spaces/slashes are replaced with '.' before matching
ct_map <- c(
  "Astrocyte"       = "Astro",
  "Endothelial"     = "Endo",
  "Microglia"       = "Micro",
  "Oligodendrocyte" = "Oligo",
  "OPC"             = "OPC",
  "PVALB"           = "Pvalb",
  "VIP"             = "Vip",
  "LAMP5"           = "Lamp5",
  "SST"             = "Sst",
  "PAX6"            = "Pax6",
  "L4.IT"           = "L4.IT",
  "L5.ET"           = "L5.ET",
  "L5.6.NP"         = "L5.6.NP",
  "L6.CT"           = "L6.CT",
  "L6b"             = "L6b",
  "L5.6.IT.Car3"    = "L6.IT.Car3"   # closest match in snRNA
  # IT, Pericyte, VLMC: no clear 1-to-1 equivalent in this snRNA-seq dataset
)

# ── 1. Load bulk MGP estimates ────────────────────────────────────────────────
cat("Loading bulk MGP estimates...\n")
mgp <- read_csv(mgp_file, show_col_types = FALSE)
cat("  Samples:", nrow(mgp), "| Cell types:", ncol(mgp) - 1, "\n")

# ── 2. Load snRNA-seq proportions, harmonise names ───────────────────────────
cat("Loading snRNA-seq proportions...\n")
sn_raw <- read_csv(sn_file, show_col_types = FALSE)
# Standardise column names: spaces and slashes → dots
colnames(sn_raw) <- gsub("[ /]", ".", colnames(sn_raw))
cat("  snRNA-seq donors:", nrow(sn_raw), "\n")

# ── 3. Build synapse ID → projid lookup via provenance + metadata ─────────────
# cell_proportions.csv uses synapse file IDs (e.g. "syn4212535") as specimenID.
# Provenance files map: synapse id → ROSMAP specimenID → projid (via metadata).
cat("Building synapse ID → projid lookup...\n")
prov <- bind_rows(lapply(prov_files, read_csv, show_col_types = FALSE))
# prov columns: id (synapse ID), specimenID (ROSMAP specimen)
cat("  Provenance records:", nrow(prov), "\n")

meta <- read_csv(meta_file, show_col_types = FALSE)
# meta columns include specimenID (ROSMAP) and projid
spec_to_projid <- meta %>%
  filter(!is.na(projid)) %>%
  select(specimenID, projid) %>%
  distinct(specimenID, .keep_all = TRUE) %>%
  mutate(projid = as.character(projid))

# synapse ID → projid
synid_to_projid <- prov %>%
  select(id, specimenID) %>%
  distinct(id, .keep_all = TRUE) %>%
  inner_join(spec_to_projid, by = "specimenID") %>%
  select(synapseID = id, projid)
cat("  Synapse IDs mapped to projid:", nrow(synid_to_projid), "\n")

# Add projid to snRNA-seq data
sn_matched <- sn_raw %>%
  mutate(ID = as.character(ID)) %>%
  inner_join(synid_to_projid, by = c("ID" = "projid"))
# sn_matched now has synapseID which matches mgp$specimenID
cat("  snRNA-seq donors with matched synapseID:", nrow(sn_matched), "\n")

# ── 4. Merge MGP and snRNA-seq on synapse ID ──────────────────────────────────
# Prefix snRNA-seq cell-type columns with "sn_" to avoid name clashes with MGP
sn_ct_cols <- setdiff(colnames(sn_matched), c("...1", "ID", "synapseID"))
colnames(sn_matched)[colnames(sn_matched) %in% sn_ct_cols] <-
  paste0("sn_", colnames(sn_matched)[colnames(sn_matched) %in% sn_ct_cols])

# Update ct_map values to use "sn_" prefix
ct_map <- setNames(paste0("sn_", ct_map), names(ct_map))

merged <- inner_join(mgp, sn_matched, by = c("specimenID" = "synapseID"))
cat("  Samples in final merged dataset:", nrow(merged), "\n")

# ── 5. Build per-cell-type scatter plots ─────────────────────────────────────
cat("Building scatter plots...\n")

# Shared theme
large_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 13, face = "bold"),
    plot.caption = element_text(size = 10, hjust = 0),
    panel.grid.minor = element_blank()
  )

cor_rows  <- list()
plot_list <- list()

for (mgp_ct in names(ct_map)) {
  sn_ct <- ct_map[[mgp_ct]]

  if (!mgp_ct %in% colnames(merged)) {
    cat("  SKIP", mgp_ct, "— not found in MGP file\n"); next
  }
  if (!sn_ct %in% colnames(merged)) {
    cat("  SKIP", mgp_ct, "→", sn_ct, "— snRNA column missing\n"); next
  }

  x <- merged[[sn_ct]]   # snRNA direct proportion (x-axis)
  y <- merged[[mgp_ct]]  # bulk MGP estimate       (y-axis)
  ok <- complete.cases(x, y)

  if (sum(ok) < 5) {
    cat("  SKIP", mgp_ct, "— too few complete cases (", sum(ok), ")\n"); next
  }

  ct_test <- cor.test(x[ok], y[ok], method = "pearson")
  r_val   <- round(ct_test$estimate, 3)
  p_val   <- ct_test$p.value
  n_val   <- sum(ok)

  # Format p-value label
  p_label <- if (p_val < 0.001) {
    formatC(p_val, format = "e", digits = 2)
  } else {
    round(p_val, 3)
  }
  annot_label <- paste0("r = ", r_val, "\np = ", p_label, "\nn = ", n_val)

  # Pretty cell-type title
  ct_title <- gsub("\\.", " ", mgp_ct)
  ct_title <- gsub("L5 6", "L5/6", ct_title)

  plot_df <- data.frame(x = x[ok], y = y[ok])

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, size = 2, colour = "#2171b5") +
    geom_smooth(method = "lm", se = TRUE, colour = "#ef6548",
                fill = "#ef6548", alpha = 0.15, linewidth = 1.2) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.08, vjust = 1.4,
             size  = 4,
             fontface = "bold",
             label = annot_label) +
    labs(
      title = ct_title,
      x     = "snRNA-seq proportion (direct count)",
      y     = "Bulk MGP estimate (z-score)"
    ) +
    large_theme

  # Save individual plot
  out_file <- file.path(outdir,
    paste0("scatter_mgp_vs_snrnaseq_", gsub("\\.", "_", mgp_ct), ".png"))
  ggsave(out_file, p, width = 7, height = 6, dpi = 200)

  plot_list[[mgp_ct]] <- p
  cor_rows[[mgp_ct]]  <- data.frame(
    cell_type    = mgp_ct,
    sn_cell_type = sn_ct,
    n            = n_val,
    pearson_r    = r_val,
    p_value      = p_val
  )
  cat(sprintf("  %-20s  r = %6.3f  p = %.3g  n = %d\n",
              mgp_ct, r_val, p_val, n_val))
}

# ── 6. Save correlation table ─────────────────────────────────────────────────
cor_df <- do.call(rbind, cor_rows)
write_tsv(cor_df,
  file.path(outdir, "scatter_mgp_vs_snrnaseq_correlations.tsv"))
cat("\nCorrelation table saved.\n")

# ── 7. Multi-panel grid ───────────────────────────────────────────────────────
cat("Saving multi-panel grid...\n")
n_plots <- length(plot_list)
n_cols  <- 4
n_rows  <- ceiling(n_plots / n_cols)

grid_plot <- plot_grid(
  plotlist = plot_list,
  ncol     = n_cols,
  nrow     = n_rows,
  align    = "hv"
)

ggsave(
  file.path(outdir, "scatter_mgp_vs_snrnaseq_all_celltypes.png"),
  grid_plot,
  width  = n_cols * 7,
  height = n_rows * 6,
  dpi    = 200,
  limitsize = FALSE
)

cat("\n========================================\n")
cat("All outputs written to:", outdir, "\n")
cat("  scatter_mgp_vs_snrnaseq_all_celltypes.png\n")
cat("  scatter_mgp_vs_snrnaseq_<celltype>.png (", n_plots, "cell types)\n")
cat("  scatter_mgp_vs_snrnaseq_correlations.tsv\n")
cat("========================================\n")
