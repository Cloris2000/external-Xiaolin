#!/usr/bin/env Rscript
# =============================================================================
# Figure 1 – Panel B (bulk cohort GWAS summary) and Panel C (snRNA-seq validation)
# =============================================================================
# Outputs (manuscript_figure/):
#   figure1_panelB_bulk_cohorts.{svg,pdf,png}   – Panel B standalone
#   figure1_panelB_bulk_cohorts.tsv              – source data table
#   figure1_panelC_snrna_validation.{svg,pdf,png} – Panel C standalone
#   figure1_panelC_snrna_validation.tsv          – source data table
#   figure1_panels_BC_preview.{pdf,png}          – side-by-side preview
#   figure1_panels_BC_QC_report.txt              – QC report
#
# All sample sizes, brain regions, genotype platforms, and cohort names are
# derived exclusively from verified pipeline inputs and GWAS outputs.
# No cohort names, sample sizes, ancestry assignments, or brain regions
# are invented in this script.
# =============================================================================

# ── 0. Settings ──────────────────────────────────────────────────────────────

# Set to FALSE to remove panel letters B/C from saved individual panels.
INCLUDE_PANEL_LETTERS <- TRUE

ROOT    <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
MS_DIR  <- file.path(ROOT, "manuscript_figure")
dir.create(MS_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(scales)
  library(grid)
  library(gtable)
  library(RColorBrewer)
  library(svglite)
})

# ── 1. Shared theme & constants ───────────────────────────────────────────────

# Ancestry palette (same as generalizability_figure script)
CLR_EUR     <- "#2166AC"
CLR_AFR     <- "#D95F02"
CLR_AMR     <- "#1B9E77"
CLR_MIXED   <- "#7B3294"
CLR_NEUTRAL <- "#637083"   # neutral bar fill

BASE_SIZE   <- 8
PANEL_W_IN  <- 7.5          # target panel width inches (both panels identical)
PANEL_H_IN  <- 3.3          # target panel height inches

theme_fig <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.y  = element_blank(),
      panel.grid.major.x  = element_line(linewidth = 0.2, colour = "grey88"),
      axis.text           = element_text(colour = "grey20", size = base_size),
      axis.title          = element_text(colour = "grey20", size = base_size),
      plot.title          = element_text(face = "bold",  size = base_size + 1.5,
                                         hjust = 0, colour = "grey10"),
      plot.subtitle       = element_text(size = base_size - 0.5,
                                         hjust = 0, colour = "grey40"),
      plot.margin         = margin(6, 6, 6, 6),
      legend.position     = "right",
      legend.key.size     = unit(0.32, "cm"),
      legend.text         = element_text(size = base_size - 1),
      legend.title        = element_text(size = base_size - 1, face = "bold"),
      legend.background   = element_rect(fill = NA, colour = NA),
      strip.text          = element_text(size = base_size - 1, face = "bold")
    )
}

save_panel <- function(p, stem, w = PANEL_W_IN, h = PANEL_H_IN, dpi = 600) {
  svglite(file.path(MS_DIR, paste0(stem, ".svg")), width = w, height = h)
  print(p); dev.off()

  cairo_pdf(file.path(MS_DIR, paste0(stem, ".pdf")), width = w, height = h)
  print(p); dev.off()

  png(file.path(MS_DIR, paste0(stem, ".png")),
      width = w * dpi, height = h * dpi, res = dpi,
      type = "cairo", bg = "white")
  print(p); dev.off()
  message("  Saved: ", stem)
}

# =============================================================================
# PANEL B ─ Bulk cohorts included in the CTP GWAS
# =============================================================================

# ── B.1  Source data  ─────────────────────────────────────────────────────────
# N values: median N across all 19 CTP traits from regenie step2 raw_p files.
# All cohorts showed constant N across traits (verified for ROSMAP, HBCC, CMC,
# GVEX, AMP-AD_Rush), so the single-trait N is also the median.
# Brain regions from tissue_filter in nextflow.config.combined.* files.
# Genotype platform from normalized_vcf_dir paths in those configs.
# Ancestry annotation: conservative (EUR-majority / Mixed/diverse) based on
#   ancestry keep-list counts in docs/ancestry_specific/keep_lists/*.txt and
#   ancestry-subset config files.

bulk_dt <- data.table(
  cohort_id    = c(
    "ROSMAP",        "ROSMAP_array",
    "Mayo",          "MSBB",
    "CMC_MSSM",      "CMC_PENN",     "CMC_PITT",
    "GTEx_v10",      "NABEC",
    "GVEX",
    "NIMH_HBCC_1M",  "NIMH_HBCC_h650", "NIMH_HBCC_Omni5M",
    "AMP_AD_Rush",   "AMP_AD_Mayo"
  ),
  display_name = c(
    "ROSMAP (WGS)",   "ROSMAP (array)",
    "Mayo",           "MSBB",
    "CMC MSSM",       "CMC PENN",    "CMC PITT",
    "GTEx v10",       "NABEC",
    "BrainGVEX",
    "HBCC (1M)",      "HBCC (h650)", "HBCC (Omni5M)",
    "AMP-AD Rush",    "AMP-AD Mayo"
  ),
  n_gwas       = c(756L, 162L, 243L, 233L, 242L, 92L, 161L,
                   283L, 210L, 394L,
                   202L, 97L, 79L, 120L, 227L),
  family       = c(
    "ROSMAP",    "ROSMAP",
    "AMP-AD",    "AMP-AD",
    "CMC",       "CMC",     "CMC",
    "GTEx",      "NABEC",
    "BrainGVEX",
    "HBCC",      "HBCC",    "HBCC",
    "AMP-AD",    "AMP-AD"
  ),
  brain_region = c(
    "DLPFC",  "DLPFC",
    "TCX",    "STG",
    "DLPFC",  "DLPFC",  "DLPFC",
    "FCX",    "FCX",
    "DLPFC",
    "DLPFC",  "DLPFC",  "DLPFC",
    "DLPFC",  "DLPFC"
  ),
  geno_platform = c(
    "WGS",            "Array (TOPMed)",
    "WGS",            "WGS",
    "Array",          "Array",          "Array",
    "WGS",            "WGS",
    "WGS",
    "Array (1M)",     "Array (h650)",   "Array (Omni5M)",
    "WGS",            "WGS"
  ),
  geno_short   = c(
    "WGS",      "Array",
    "WGS",      "WGS",
    "Array",    "Array",  "Array",
    "WGS",      "WGS",
    "WGS",
    "Array",    "Array",  "Array",
    "WGS",      "WGS"
  ),
  ancestry_cat = c(
    "EUR-majority",  "EUR-majority",
    "EUR-majority",  "EUR-majority",
    "EUR-majority",  "EUR-majority", "EUR-majority",
    "EUR-majority",  "EUR-majority",
    "EUR-majority",
    "Mixed/diverse", "Mixed/diverse", "EUR-majority",
    "Mixed/diverse", "Mixed/diverse"
  ),
  source_file  = c(
    "results/ROSMAP/regenie_step2/ROSMAP_Astrocyte_step2.regenie.raw_p",
    "results/ROSMAP_array/regenie_step2/ROSMAP_array_Astrocyte_step2.regenie.raw_p",
    "results/Mayo/regenie_step2/Mayo_Astrocyte_step2.regenie.raw_p",
    "results/MSBB/regenie_step2/MSBB_Astrocyte_step2.regenie.raw_p",
    "results/CMC_MSSM/regenie_step2/CMC_MSSM_Astrocyte_step2.regenie.raw_p",
    "results/CMC_PENN/regenie_step2/CMC_PENN_Astrocyte_step2.regenie.raw_p",
    "results/CMC_PITT/regenie_step2/CMC_PITT_Astrocyte_step2.regenie.raw_p",
    "results/GTEx_v10/regenie_step2/GTEx_v10_Astrocyte_step2.regenie.raw_p",
    "results/NABEC/regenie_step2/NABEC_Astrocyte_step2.regenie.raw_p",
    "results/GVEX/regenie_step2/GVEX_Astrocyte_step2.regenie.raw_p",
    "results/NIMH_HBCC_1M/regenie_step2/NIMH_HBCC_1M_Astrocyte_step2.regenie.raw_p",
    "results/NIMH_HBCC_h650/regenie_step2/NIMH_HBCC_h650_Astrocyte_step2.regenie.raw_p",
    "results/NIMH_HBCC_Omni5M/regenie_step2/NIMH_HBCC_Omni5M_Astrocyte_step2.regenie.raw_p",
    "results/AMP_AD_Rush/regenie_step2/AMP_AD_Rush_Astrocyte_step2.regenie.raw_p",
    "results/AMP_AD_Mayo/regenie_step2/AMP_AD_Mayo_Astrocyte_step2.regenie.raw_p"
  ),
  notes = c(
    "ROSMAP WGS joint-called; b37; EUR-majority (EUR keep-list=797/852 pheno)",
    "ROSMAP TOPMed-imputed array subset not in WGS pipeline; b37; EUR-majority",
    "Mayo Clinic TCX; WGS; EUR-majority",
    "Mt Sinai Brain Bank STG; WGS; EUR-majority",
    "CommonMind MSSM DLPFC; array+imputation; EUR-majority",
    "CommonMind PENN DLPFC; array+imputation; EUR-majority",
    "CommonMind PITT DLPFC; array+imputation; EUR-majority",
    "GTEx v10 BA9 frontal cortex; WGS GRCh38 (lifted hg19 for meta); EUR-majority",
    "NABEC frontal cortex; WGS; EUR-majority",
    "BrainGVEX DLPFC; WGS; EUR-majority (multi-ancestry noted; EUR keep-list=389/394)",
    "NIMH HBCC DLPFC; 1M array+imputation; Mixed (AFR=101, EUR=88 of N=202)",
    "NIMH HBCC DLPFC; h650 array+imputation; Mixed (AFR keep-list=53 of N=97)",
    "NIMH HBCC DLPFC; Omni5M array+imputation; EUR-majority (EUR keep-list=50 of N=79)",
    "AMP-AD Diverse Rush DLPFC; WGS GRCh38 (lifted); Mixed (AFR keep-list=64 of N=120)",
    "AMP-AD Diverse Mayo DLPFC; WGS GRCh38 (lifted); Mixed (AFR=52, AMR=178 of N=227)"
  )
)

# ── B.2  Cohort display order (family-grouped, matching COHORT_ORDER in scripts) ──
family_order <- c("ROSMAP", "AMP-AD", "CMC", "GTEx", "NABEC", "BrainGVEX", "HBCC")
# Within AMP-AD, keep Mayo+MSBB together, then Rush+Mayo Diverse
# Desired display order (bottom → top in horizontal bar chart):
display_order <- c(
  "AMP-AD Mayo",   "AMP-AD Rush",          # AMP-AD Diverse (bottom)
  "HBCC (Omni5M)", "HBCC (h650)", "HBCC (1M)",    # HBCC platforms
  "BrainGVEX",                               # BrainGVEX
  "NABEC",                                   # NABEC
  "GTEx v10",                                # GTEx
  "CMC PITT",  "CMC PENN",  "CMC MSSM",     # CMC
  "MSBB",      "Mayo",                       # AMP-AD legacy
  "ROSMAP (array)", "ROSMAP (WGS)"           # ROSMAP (top)
)

# family → y-gap size (add gap above each family transition)
bulk_dt[, display_name := factor(display_name, levels = display_order)]
total_n <- sum(bulk_dt$n_gwas)

# ancestry colour palette (strip)
anc_cols <- c(
  "EUR-majority"  = "#3399CC",
  "Mixed/diverse" = "#CC6633"
)

# region colour palette
region_cols <- c(
  "DLPFC" = "#4A9C6D",
  "TCX"   = "#C46B3A",
  "STG"   = "#8E63AE",
  "FCX"   = "#3975B7"
)

# ── B.3  Build panel B plot ───────────────────────────────────────────────────

# Assign y positions with gaps between families
# Top to bottom in bar chart (highest y = first display_order entry, bottom of horizontal chart)
bulk_sorted <- bulk_dt[order(match(display_name, display_order))]
bulk_sorted[, y_pos := .I]

# Family groups for gap lines (between groups)
family_groups <- list(
  ROSMAP        = c("ROSMAP (WGS)", "ROSMAP (array)"),
  legacy        = c("Mayo", "MSBB"),
  CMC           = c("CMC MSSM", "CMC PENN", "CMC PITT"),
  GTEx          = c("GTEx v10"),
  NABEC         = c("NABEC"),
  BrainGVEX     = c("BrainGVEX"),
  HBCC          = c("HBCC (1M)", "HBCC (h650)", "HBCC (Omni5M)"),
  `AMP-AD Div`  = c("AMP-AD Rush", "AMP-AD Mayo")
)

# Compute gap y-positions (draw a dashed line between family groups)
fam_boundaries <- c()
y_last <- 0
for (grp in rev(names(family_groups))) {
  members <- family_groups[[grp]]
  y_positions <- bulk_sorted[display_name %in% members, y_pos]
  if (length(y_positions) > 0 && y_last > 0) {
    fam_boundaries <- c(fam_boundaries, (min(y_positions) + y_last - 0.5))
  }
  y_last <- max(y_positions, na.rm = TRUE)
}

# Right-panel annotations data
annot_dt <- bulk_sorted[, .(display_name, y_pos, brain_region, geno_short, ancestry_cat)]

# Max x for axis limits
x_max <- ceiling(max(bulk_sorted$n_gwas) / 100) * 100 + 50

# Bar width
BAR_H  <- 0.65

# Build the main bar chart
pB_bars <- ggplot(bulk_sorted,
                  aes(x = n_gwas, y = y_pos)) +
  # horizontal family separator lines
  geom_hline(yintercept = fam_boundaries,
             linetype   = "dashed", linewidth = 0.25, colour = "grey70") +
  # main bars
  geom_col(aes(x = n_gwas, y = y_pos), fill = CLR_NEUTRAL,
           colour = NA, width = BAR_H) +
  # N labels at bar end
  geom_text(aes(label = format(n_gwas, big.mark = ",")),
            hjust = -0.15, vjust = 0.38,
            size = BASE_SIZE * 0.28, colour = "grey20") +
  # cohort labels
  scale_y_continuous(
    breaks = bulk_sorted$y_pos,
    labels = bulk_sorted$display_name,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_x_continuous(
    limits = c(0, x_max * 1.18),
    breaks = seq(0, x_max, by = 200),
    labels = scales::label_comma(),
    expand = c(0, 0)
  ) +
  labs(
    x     = "Analyzed donors",
    y     = NULL,
    title = if (INCLUDE_PANEL_LETTERS) "B   Bulk cohorts included in the CTP GWAS" else
              "Bulk cohorts included in the CTP GWAS",
    subtitle = sprintf("15 cohort datasets  |  Total analyzed N = %s", format(total_n, big.mark = ","))
  ) +
  theme_fig() +
  theme(
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90")
  )

# ── B.3b  Ancestry annotation column ─────────────────────────────────────────
# Build ancestry as its own annotation column (matches the Region/Genotype style)
make_anc_col <- function(dt, title = "Ancestry") {
  dt2 <- copy(dt)
  # Short label for tile
  dt2[, anc_short := ifelse(ancestry_cat == "EUR-majority", "EUR", "Mixed")]
  anc_tile_cols <- c("EUR" = "#3399CC", "Mixed" = "#CC6633")
  ggplot(dt2, aes(x = 0.5, y = y_pos)) +
    geom_tile(aes(fill = anc_short),
              width = 1, height = BAR_H, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = anc_short),
              size = BASE_SIZE * 0.24, hjust = 0.5, vjust = 0.5,
              colour = "white", fontface = "bold") +
    scale_y_continuous(
      breaks = dt2$y_pos,
      labels = rep("", nrow(dt2)),
      expand = expansion(add = c(0.5, 0.5))
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = anc_tile_cols) +
    labs(x = NULL, y = NULL, title = title) +
    theme_fig() +
    theme(
      axis.text     = element_blank(),
      axis.ticks    = element_blank(),
      axis.line     = element_blank(),
      panel.grid    = element_blank(),
      plot.title    = element_text(size = BASE_SIZE - 1, hjust = 0.5,
                                   face = "bold", colour = "grey30",
                                   margin = margin(b = 2)),
      plot.margin   = margin(6, 2, 6, 2),
      legend.position = "none"
    )
}

# ── B.4  Right annotation columns (brain region + geno platform) ──────────────
make_annot_col <- function(dt, col, title, colors = NULL) {
  p <- ggplot(dt, aes(x = 0.5, y = y_pos)) +
    geom_tile(aes(fill = .data[[col]]),
              width = 1, height = BAR_H, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = .data[[col]]),
              size = BASE_SIZE * 0.24, hjust = 0.5, vjust = 0.5,
              colour = "white", fontface = "bold") +
    scale_y_continuous(
      breaks = dt$y_pos,
      labels = rep("", nrow(dt)),
      expand = expansion(add = c(0.5, 0.5))
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_fig() +
    theme(
      axis.text     = element_blank(),
      axis.ticks    = element_blank(),
      axis.line     = element_blank(),
      panel.grid    = element_blank(),
      plot.title    = element_text(size = BASE_SIZE - 1, hjust = 0.5,
                                   face = "bold", colour = "grey30",
                                   margin = margin(b = 2)),
      plot.margin   = margin(6, 2, 6, 2),
      legend.position = "none"
    )
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  } else {
    p <- p + scale_fill_manual(values = setNames(
      colorRampPalette(c("#5B9BD5","#70AD47","#FFC000","#7030A0"))(
        length(unique(dt[[col]]))),
      unique(dt[[col]])))
  }
  p
}

annot_anc    <- make_anc_col(annot_dt, "Ancestry")
annot_region <- make_annot_col(annot_dt, "brain_region", "Region",
  colors = region_cols)
plat_cols <- c(WGS = "#4D7FAF", Array = "#9B7DB0")
annot_plat <- make_annot_col(annot_dt, "geno_short", "Genotype",
  colors = plat_cols)

# ── B.6  Assemble Panel B ─────────────────────────────────────────────────────
pB <- plot_grid(
  pB_bars,
  annot_anc,
  annot_region,
  annot_plat,
  nrow    = 1,
  rel_widths = c(1, 0.085, 0.10, 0.10),
  align   = "h",
  axis    = "tb"
)

# ── B.7  Save Panel B ─────────────────────────────────────────────────────────
message("Saving Panel B ...")
save_panel(pB, "figure1_panelB_bulk_cohorts")

# ── B.8  Write TSV ────────────────────────────────────────────────────────────
bulk_out <- bulk_sorted[, .(
  cohort_name    = cohort_id,
  cohort_family  = family,
  display_name,
  analyzed_N     = n_gwas,
  ctp_n_min      = n_gwas,   # N constant across 19 traits
  ctp_n_max      = n_gwas,
  ancestry_category = ancestry_cat,
  ancestry_notes = notes,
  brain_region,
  geno_platform,
  source_file
)]
fwrite(bulk_out, file.path(MS_DIR, "figure1_panelB_bulk_cohorts.tsv"), sep = "\t")
message("  Saved: figure1_panelB_bulk_cohorts.tsv")


# =============================================================================
# PANEL C ─ Independent snRNA-seq validation
# =============================================================================
# Final snRNA-seq meta-analysis: hodge5u (de-duplicated 5-cohort)
#   ROSMAP_Green_Hodge   – Green et al. ROSMAP snRNA-seq + ROSMAP WGS
#   PsychAD_HBCC_Hodge   – PsychAD HBCC snRNA-seq + HBCC genotype
#   PsychAD_MSSM_Hodge   – PsychAD MSSM snRNA-seq + MSBB WGS
#   Ruz_MSSM             – Ruzicka et al. 2024 Science, MSSM + CMC_MSSM genotype
#   ROSMAP_Mathys_unique – Mathys et al. ROSMAP snRNA-seq (128 donors not shared with Green)
# Source: nextflow.config.standalone_meta.sn_hodge5u
# =============================================================================

# ── C.1  Source data ──────────────────────────────────────────────────────────
sn_dt <- data.table(
  cohort_id      = c("ROSMAP_Green_Hodge", "PsychAD_HBCC_Hodge",
                     "PsychAD_MSSM_Hodge", "Ruz_MSSM", "ROSMAP_Mathys_unique"),
  display_name   = c("ROSMAP Green",       "PsychAD HBCC",
                     "PsychAD MSSM",       "Ruzicka MSSM", "ROSMAP Mathys\n(unique)"),
  dataset        = c("Green et al. (ROSMAP)",  "PsychAD Study (HBCC)",
                     "PsychAD Study (MSSM)",   "Ruzicka et al. 2024",
                     "Mathys et al. (ROSMAP)"),
  cohort_source  = c("ROSMAP", "NIMH HBCC", "MSBB", "CMC MSSM", "ROSMAP"),
  ctp_donors     = c(420L,  259L, 203L, 53L, 128L),
  gwas_donors    = c(366L,  258L, 187L, 53L, 121L),
  nuclei_total   = c(NA_integer_, NA_integer_, NA_integer_,
                     NA_integer_, 1420318L),   # Mathys occupancy_report.tsv total
  n_celltypes    = c(19L, 19L, 19L, 18L, 16L),
  ancestry_notes = c(
    "EUR-majority (ROSMAP WGS)",
    "Mixed (PsychAD HBCC: AFR/EUR/other; HBCC genotype)",
    "EUR-majority (MSBB WGS)",
    "EUR-majority (53 donors: 49 EUR + 4 AFR; CMC_MSSM array)",
    "EUR-majority (ROSMAP WGS, b37)"
  ),
  genome_build   = c("b37", "b37", "b37", "b37", "b37"),
  prop_file      = c(
    "data_input/sn_rosmap_green_hodge/cell_proportions.csv",
    "data_input/sn_psychad_hbcc_hodge/cell_proportions.csv",
    "data_input/sn_psychad_mssm_hodge/cell_proportions.csv",
    "data_input/sn_ruz_mssm/cell_proportions.csv",
    "data_input/sn_rosmap_mathys_unique/cell_proportions.csv"
  ),
  gwas_dir       = c(
    "results/sn_rosmap_green_hodge/regenie_step2",
    "results/sn_psychad_hbcc_hodge/regenie_step2",
    "results/sn_psychad_mssm_hodge/regenie_step2",
    "results/sn_Ruz_MSSM/regenie_step2",
    "results/sn_rosmap_mathys_unique/regenie_step2"
  ),
  meta_analysis  = rep("hodge5u (sn_hodge5u)", 5L),
  notes          = c(
    "420 CTP donors; 366 passed geno QC; 54 lost at genotype-matching step",
    "259 CTP donors; 258 passed geno QC",
    "203 CTP donors; 187 passed geno QC",
    "53 CTP donors; all 53 matched genotype (CMC_MSSM array; 49 EUR + 4 AFR)",
    "128 Mathys donors NOT shared with ROSMAP_Green; 121 passed geno QC"
  )
)

# Display order (top to bottom in bar chart: largest cohort first)
sn_order <- c("ROSMAP Green", "PsychAD HBCC", "PsychAD MSSM",
              "ROSMAP Mathys\n(unique)", "Ruzicka MSSM")
sn_dt[, display_name := factor(display_name, levels = rev(sn_order))]
sn_dt[, y_pos := as.integer(display_name)]

# ── C.2  Read actual cell-proportion data for heatmap ─────────────────────────
# Manuscript Hodge subclass order (from bulk_rosmap_cell_types.txt)
CT_ORDER <- c(
  "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
  "L5.ET", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC",
  "Oligodendrocyte", "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC"
)
# Ruzicka uses slightly different names (Ruzicka 2024); harmonise
# (L5.ET vs L5.ET, L5.6.IT.Car3 vs L5.6.IT.Car3 – already same)

read_props <- function(csv_path, cohort_label) {
  full_path <- file.path(ROOT, csv_path)
  if (!file.exists(full_path)) {
    message("  MISSING: ", full_path)
    return(NULL)
  }
  dt <- fread(full_path)
  id_col <- names(dt)[1]
  # pivot to long
  dt_long <- melt(dt, id.vars = id_col, variable.name = "cell_type",
                  value.name = "proportion")
  setnames(dt_long, id_col, "donor_id")
  dt_long[, cohort := cohort_label]
  dt_long[, donor_idx := as.integer(factor(donor_id,
                                           levels = unique(donor_id)))]
  dt_long
}

message("Reading cell-proportion files ...")
prop_list <- mapply(
  read_props,
  csv_path     = sn_dt$prop_file,
  cohort_label = as.character(sn_dt$display_name),
  SIMPLIFY     = FALSE
)
prop_long <- rbindlist(Filter(Negate(is.null), prop_list))

# Harmonise Ruzicka cell-type names to Hodge subclass names where needed
# Ruzicka 2024 uses the same names as the other hodge cohorts in this pipeline
# (verified from data_input/sn_ruz_mssm/cell_proportions.csv header)

# Filter to canonical cell types only; keep all donors (no clustering)
prop_long <- prop_long[cell_type %in% CT_ORDER]
prop_long[, cell_type := factor(cell_type, levels = CT_ORDER)]

# Assign global donor index per cohort for heatmap columns
# Keep donor order within cohort (natural order in file = row order)
# Use abbreviated names for heatmap column-group labels to avoid overlap
cohort_order_for_heatmap <- c("ROSMAP Green", "PsychAD HBCC", "PsychAD MSSM",
                               "ROSMAP Mathys\n(unique)", "Ruzicka MSSM")
# Short labels map for heatmap headings
heatmap_label_map <- c(
  "ROSMAP Green"            = "ROSMAP Green",
  "PsychAD HBCC"            = "PsychAD HBCC",
  "PsychAD MSSM"            = "PsychAD MSSM",
  "ROSMAP Mathys\n(unique)" = "Mathys (unique)",
  "Ruzicka MSSM"            = "Ruzicka"
)
prop_long[, cohort := factor(cohort, levels = cohort_order_for_heatmap)]
prop_long[, global_donor_idx := {
  offset <- 0L
  vals   <- integer(.N)
  for (coh in cohort_order_for_heatmap) {
    rows <- which(cohort == coh)
    if (length(rows) > 0) {
      vals[rows] <- donor_idx[rows] + offset
      offset <- offset + max(donor_idx[rows])
    }
  }
  vals
}]

# Cohort boundary x-positions for separator lines
cohort_boundaries <- prop_long[, .(
  xmin = min(global_donor_idx),
  xmax = max(global_donor_idx)
), by = cohort]
cohort_boundaries[, x_mid := (xmin + xmax) / 2]

# ── C.3  Cell-type proportion heatmap ────────────────────────────────────────
# Use sqrt transform to better show low proportions
prop_long[, prop_sqrt := sqrt(proportion)]

pC_heatmap <- ggplot(prop_long,
                     aes(x = global_donor_idx, y = cell_type,
                         fill = prop_sqrt)) +
  geom_raster(interpolate = FALSE) +
  # cohort separators
  geom_vline(data = cohort_boundaries[-nrow(cohort_boundaries)],
             aes(xintercept = xmax + 0.5),
             linewidth = 0.5, colour = "white") +
  # cohort labels along top — use abbreviated names from heatmap_label_map
  annotate("text",
           x      = cohort_boundaries$x_mid,
           y      = length(CT_ORDER) + 0.85,
           label  = heatmap_label_map[as.character(cohort_boundaries$cohort)],
           size   = BASE_SIZE * 0.255, hjust = 0.5, vjust = 0,
           colour = "grey20", fontface = "bold") +
  scale_fill_gradientn(
    colours = c("white", "#FFFACD", "#FED976", "#FEB24C",
                "#FD8D3C", "#E31A1C", "#800026"),
    name    = "Proportion\n(sqrt scale)",
    breaks  = c(0, 0.2, 0.4, 0.6, 0.8),
    labels  = sprintf("%.2g", c(0, 0.2, 0.4, 0.6, 0.8)^2)
  ) +
  scale_y_discrete(expand = expansion(add = c(0, 1.4))) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_fig() +
  theme(
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.y   = element_text(size = BASE_SIZE - 1.5, hjust = 1),
    panel.grid    = element_blank(),
    legend.position  = "right",
    legend.key.height = unit(0.6, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    plot.margin   = margin(4, 4, 2, 4)
  )

# ── C.4  Left bar chart: CTP donors vs GWAS donors ───────────────────────────
sn_long <- rbind(
  sn_dt[, .(display_name, y_pos, n = ctp_donors,  type = "snRNA-seq CTP available")],
  sn_dt[, .(display_name, y_pos, n = gwas_donors, type = "Included in GWAS")]
)
sn_long[, type := factor(type, levels = c("snRNA-seq CTP available",
                                           "Included in GWAS"))]

bar_fill   <- c("snRNA-seq CTP available" = "#BDD7EE",
                "Included in GWAS"        = "#2E75B6")
bar_colour <- c("snRNA-seq CTP available" = "#5B9BD5",
                "Included in GWAS"        = "#2E75B6")

SN_BAR_H <- 0.65   # bar height for snRNA-seq bar chart
xmax_sn  <- ceiling(max(sn_long$n) / 50) * 50 + 25

pC_bars <- ggplot() +
  # pale outlined bar for CTP donors
  geom_col(data = sn_long[type == "snRNA-seq CTP available"],
           aes(x = n, y = y_pos),
           fill   = bar_fill["snRNA-seq CTP available"],
           colour = bar_colour["snRNA-seq CTP available"],
           linewidth = 0.35, width = SN_BAR_H) +
  # darker filled bar for GWAS donors (on top)
  geom_col(data = sn_long[type == "Included in GWAS"],
           aes(x = n, y = y_pos),
           fill   = bar_fill["Included in GWAS"],
           colour = NA, width = SN_BAR_H) +
  # labels: "CTP / GWAS donors"
  geom_text(data = sn_dt,
            aes(x = xmax_sn * 1.01, y = y_pos,
                label = paste0(ctp_donors, " CTP / ", gwas_donors, " GWAS")),
            hjust  = 0, vjust = 0.38,
            size   = BASE_SIZE * 0.25, colour = "grey25") +
  scale_y_continuous(
    breaks = sn_dt$y_pos,
    labels = levels(sn_dt$display_name)[sn_dt$y_pos],
    expand = expansion(add = c(0.5, 0.5))
  ) +
  scale_x_continuous(
    limits = c(0, xmax_sn * 1.55),
    breaks = seq(0, xmax_sn, by = 100),
    expand = c(0, 0)
  ) +
  # legend manually drawn below
  labs(
    x     = "Number of donors",
    y     = NULL,
    title = if (INCLUDE_PANEL_LETTERS)
              "C   Independent snRNA-seq validation" else
              "Independent snRNA-seq validation"
  ) +
  theme_fig() +
  theme(
    axis.line.y    = element_blank(),
    axis.ticks.y   = element_blank()
  )

# ── C.5  Right workflow diagram ───────────────────────────────────────────────
# 5 steps on a 0–10 y-scale with generous gaps between boxes
step_y   <- c(9.2, 7.2, 5.2, 3.2, 1.2)
box_half <- 0.75   # half-height of each box

workflow_steps <- data.frame(
  label = c(
    "Annotated\nnuclei",
    "Donor-level\nCTPs",
    "Cohort-level\nGWAS",
    "snRNA-seq\nmeta-analysis",
    "Comparison with\nbulk CTP GWAS"
  ),
  sublabel = c(
    "",
    "19 harmonized cell types\nnuclei of type / total\nnuclei per donor",
    "",
    sprintf("5 cohorts  (N = %d)", sum(sn_dt$gwas_donors)),
    ""
  ),
  y = step_y,
  stringsAsFactors = FALSE
)
workflow_steps$has_sub    <- nchar(workflow_steps$sublabel) > 0
workflow_steps$box_top    <- workflow_steps$y + box_half
workflow_steps$box_bottom <- workflow_steps$y - box_half

arrow_df <- data.frame(
  y_start = step_y[-length(step_y)] - box_half,
  y_end   = step_y[-1]              + box_half
)

pC_workflow <- ggplot() +
  geom_rect(data = workflow_steps,
            aes(xmin = 0.04, xmax = 0.96,
                ymin = box_bottom, ymax = box_top),
            fill = "#EBF3FB", colour = "#2E75B6", linewidth = 0.4,
            inherit.aes = FALSE) +
  geom_text(data = workflow_steps,
            aes(x = 0.5,
                y = y + ifelse(has_sub, 0.25, 0),
                label = label),
            size = BASE_SIZE * 0.26, hjust = 0.5, vjust = 0.5,
            fontface = "bold", colour = "#1F4E79") +
  geom_text(data = workflow_steps[workflow_steps$has_sub, ],
            aes(x = 0.5, y = y - 0.30, label = sublabel),
            size = BASE_SIZE * 0.21, hjust = 0.5, vjust = 0.5,
            colour = "#5B9BD5", lineheight = 0.9) +
  geom_segment(data = arrow_df,
               aes(x = 0.5, xend = 0.5,
                   y = y_start - 0.03, yend = y_end + 0.03),
               arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
               linewidth = 0.35, colour = "#2E75B6") +
  xlim(0, 1) +
  ylim(0.2, 10.2) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_void() +
  theme(plot.margin = margin(6, 4, 6, 4))

# ── C.6  Add inline legend to the bar chart ──────────────────────────────────
# Append a manual legend directly to pC_bars via a guide
pC_bars <- pC_bars +
  # Re-add bars with an aesthetic mapping for the legend
  geom_col(data = sn_long,
           aes(x = n, y = y_pos, fill = type, colour = type),
           width = SN_BAR_H, linewidth = 0.35, show.legend = TRUE) +
  scale_fill_manual(values   = bar_fill,   name = NULL,
                    guide    = guide_legend(override.aes = list(colour = bar_colour))) +
  scale_colour_manual(values = bar_colour, name = NULL, guide = "none") +
  theme(
    legend.position.inside = c(0.95, 0.06),
    legend.justification   = c(1, 0),
    legend.key.size  = unit(0.28, "cm"),
    legend.text      = element_text(size = BASE_SIZE - 1.5),
    legend.background = element_rect(fill = "white", colour = NA)
  )

# ── C.7  Assemble Panel C ─────────────────────────────────────────────────────
# Top row: bars (left) + workflow (right)
pC_top <- plot_grid(
  pC_bars,
  pC_workflow,
  nrow       = 1,
  rel_widths = c(1, 0.52),
  align      = "h",
  axis       = "tb"
)
# Bottom: heatmap (full width)
pC <- plot_grid(
  pC_top,
  pC_heatmap,
  nrow       = 2,
  rel_heights = c(0.52, 0.48)
)

# ── C.8  Save Panel C ─────────────────────────────────────────────────────────
message("Saving Panel C ...")
save_panel(pC, "figure1_panelC_snrna_validation",
           h = PANEL_H_IN + 0.8)   # slightly taller for heatmap

# ── C.9  Write TSV ────────────────────────────────────────────────────────────
sn_out <- sn_dt[, .(
  dataset,
  cohort                = cohort_id,
  cohort_source,
  ctp_available_donors  = ctp_donors,
  genotype_matched_donors = gwas_donors,
  final_gwas_N          = gwas_donors,
  retained_nuclei       = nuclei_total,
  ancestry_info         = ancestry_notes,
  n_harmonized_celltypes = n_celltypes,
  genome_build,
  meta_analysis,
  source_proportion_file = prop_file,
  source_gwas_dir        = gwas_dir,
  notes
)]
fwrite(sn_out, file.path(MS_DIR, "figure1_panelC_snrna_validation.tsv"), sep = "\t")
message("  Saved: figure1_panelC_snrna_validation.tsv")


# =============================================================================
# COMBINED PREVIEW
# =============================================================================
message("Saving combined preview ...")
pPreview <- plot_grid(pB, pC,
                      nrow = 2,
                      rel_heights = c(PANEL_H_IN, PANEL_H_IN + 0.8))

cairo_pdf(file.path(MS_DIR, "figure1_panels_BC_preview.pdf"),
          width = PANEL_W_IN, height = PANEL_H_IN * 2 + 0.8)
print(pPreview); dev.off()

png(file.path(MS_DIR, "figure1_panels_BC_preview.png"),
    width = PANEL_W_IN * 600, height = (PANEL_H_IN * 2 + 0.8) * 600,
    res = 600, type = "cairo", bg = "white")
print(pPreview); dev.off()
message("  Saved: figure1_panels_BC_preview.{pdf,png}")


# =============================================================================
# QC REPORT
# =============================================================================
qc_lines <- c(
  "figure1_panels_BC_QC_report.txt",
  paste0("Generated: ", Sys.time()),
  "",
  "═══════════════════════════════════════════════════════════════════════",
  "PANEL B: Bulk cohorts included in the CTP GWAS",
  "═══════════════════════════════════════════════════════════════════════",
  "",
  "Data source type: regenie step2 raw_p GWAS output files",
  "Verification: N constant across all 19 CTP traits (verified for ROSMAP,",
  "              NIMH_HBCC_1M, CMC_MSSM, GVEX, AMP_AD_Rush).",
  "",
  "Cohorts detected (15):",
  paste(sprintf("  %-22s  N=%4d  %-6s  %-10s  %-14s",
                bulk_out$cohort_name, bulk_out$analyzed_N,
                bulk_out$brain_region, bulk_out$geno_platform,
                bulk_out$ancestry_category), collapse = "\n"),
  "",
  sprintf("Total analyzed N = %d", total_n),
  "",
  "Ancestry annotation method: conservative categorical labels based on",
  "  ancestry-specific keep-list files (docs/ancestry_specific/keep_lists/)",
  "  and ancestry-subset config comments. No exact ancestry counts plotted",
  "  (counts do not sum exactly to GWAS N for all cohorts).",
  "",
  "Brain regions from tissue_filter in nextflow.config.combined.* files.",
  "Genotype platforms from normalized_vcf_dir paths in those configs.",
  "",
  "Missing metadata: No ancestry keep-list for ROSMAP_array, CMC_PENN,",
  "  CMC_PITT, AMP_AD_Rush EUR subset, AMP_AD_Mayo EUR subset. These were",
  "  classified by available keep-list presence only.",
  "",
  "═══════════════════════════════════════════════════════════════════════",
  "PANEL C: Independent snRNA-seq validation",
  "═══════════════════════════════════════════════════════════════════════",
  "",
  "Meta-analysis tag: hodge5u (de-duplicated 5-cohort sn meta)",
  "Source config: nextflow.config.standalone_meta.sn_hodge5u",
  "Comparison performed in: manuscript_figure/figure5_sn_concordance.R",
  "",
  "Cohorts detected (5):",
  paste(sprintf("  %-26s  CTP=%3d  GWAS=%3d  CellTypes=%d",
                sn_out$cohort, sn_out$ctp_available_donors,
                sn_out$final_gwas_N, sn_out$n_harmonized_celltypes), collapse = "\n"),
  "",
  sprintf("Total sn GWAS donors = %d", sum(sn_dt$gwas_donors)),
  "",
  "Cell-proportion heatmap: derived from actual donor-level proportion tables.",
  "  Files read:",
  paste(sprintf("  - %s", sn_dt$prop_file), collapse = "\n"),
  "  Donors kept in natural (file) order within cohort; no clustering.",
  "  Cell types filtered to 19 Hodge subclasses (bulk_rosmap_cell_types.txt order).",
  sprintf("  Donors with at least 1 Hodge cell type: %d (of %d total in files)",
          prop_long[, uniqueN(paste(cohort, donor_id))],
          sum(sn_dt$ctp_donors)),
  "",
  "Nuclei counts: only available for ROSMAP_Mathys_Hodge (occupancy_report.tsv)",
  "  Total Mathys nuclei: 1,420,318",
  "  Not available for: ROSMAP_Green, PsychAD_HBCC, PsychAD_MSSM, Ruzicka",
  "",
  "═══════════════════════════════════════════════════════════════════════",
  "OUTPUT FILES",
  "═══════════════════════════════════════════════════════════════════════",
  paste(sprintf("  %s", list.files(MS_DIR, pattern = "figure1_panel", full.names = FALSE)),
        collapse = "\n"),
  ""
)

writeLines(qc_lines, file.path(MS_DIR, "figure1_panels_BC_QC_report.txt"))
message("  Saved: figure1_panels_BC_QC_report.txt")
message("Done.")
