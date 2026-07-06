#!/usr/bin/env Rscript
# Regional plot for a single locus — bulk panel on top, sn panel below.
# If the SNP is not in sn, shows bulk-only with a note.
# Works for any cell type and marker regardless of concordance status.
#
# Usage:
#   Rscript plot_regional_single_locus.R \
#     --cell_type Oligodendrocyte --marker chr7:73043561:G:A
#   Rscript plot_regional_single_locus.R \
#     --cell_type Oligodendrocyte --marker chr7:12263378:G:T

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(optparse)
})

option_list <- list(
  make_option("--cell_type", type = "character", default = "Oligodendrocyte"),
  make_option("--marker",    type = "character", default = "chr7:73043561:G:A"),
  make_option("--sn_cell_type", type = "character", default = "",
    help = "sn cell type name (auto-inferred if empty)"),
  make_option("--bulk_dir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_13cohorts_design_matrix"),
  make_option("--sn_dir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_sn_rosmap_msbb_hbcc_alias15"),
  make_option("--outdir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity_design_matrix/top_hits/regional_plots"),
  make_option("--window_kb", type = "integer", default = 500)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

HALF <- opt$window_kb * 1000L

# ---- cell type crosswalk ----------------------------------------------------
crosswalk <- c(
  Astrocyte = "ASTRO", Endothelial = "ENDO", IT = "L23IT",
  "L4.IT" = "L23IT", "L5.6.IT.Car3" = "L23IT",
  "L5.6.NP" = "L56NP", "L5.ET" = "L5ET", L6b = "L6B", "L6.CT" = "L6CT",
  LAMP5 = "LAMP5LHX6", Microglia = "PVM", Oligodendrocyte = "OLIGO",
  OPC = "OPC", PAX6 = "L23IT", Pericyte = "ENDO",
  PVALB = "PVALB", SST = "SST", VIP = "VIP", VLMC = "VLMC"
)
sn_ct <- if (opt$sn_cell_type != "") opt$sn_cell_type else crosswalk[opt$cell_type]
if (is.na(sn_ct)) sn_ct <- NULL

# ---- parse marker -----------------------------------------------------------
parts   <- strsplit(opt$marker, ":", fixed = TRUE)[[1]]
chrom   <- parts[1]
pos     <- as.integer(parts[2])
start   <- pmax(0L, pos - HALF)
end     <- pos + HALF

# ---- helpers ----------------------------------------------------------------
find_tbl <- function(dir, ct) {
  pat <- paste0("^", gsub("([.])", "\\\\.", ct), "_meta_analysis_.*\\.tbl$")
  f   <- list.files(dir, pattern = pat, full.names = TRUE)
  f   <- f[!grepl("1\\.tbl$", basename(f))]
  if (length(f) == 0) return(NA_character_)
  f[1]
}

read_region <- function(path, chrom, start, end) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  dt <- fread(path, sep = "\t", data.table = TRUE, verbose = FALSE)
  setnames(dt, tolower(gsub("[^A-Za-z0-9]", "_", names(dt))))
  if ("markername" %in% names(dt) && !("marker" %in% names(dt)))
    setnames(dt, "markername", "marker")
  if ("p_value" %in% names(dt) && !("pvalue" %in% names(dt)))
    setnames(dt, "p_value", "pvalue")
  needed <- c("marker", "allele1", "allele2", "effect", "stderr", "pvalue")
  if (!all(needed %in% names(dt))) return(NULL)
  dt <- dt[, needed, with = FALSE]
  dt <- dt[nchar(allele1) == 1 & nchar(allele2) == 1]
  parts2 <- tstrsplit(dt$marker, ":", fixed = TRUE)
  dt[, chrom_col := parts2[[1]]]
  dt[, pos2      := as.integer(parts2[[2]])]
  dt[, effect    := as.numeric(effect)]
  dt[, pvalue    := as.numeric(pvalue)]
  dt[chrom_col == chrom & pos2 >= start & pos2 <= end]
}

lz_colours <- colorRampPalette(
  c("#357EBD", "#46B8DA", "#5CB85C", "#F0AD4E", "#D9534F"))(100)

make_panel <- function(dt, lead_pos, half_win, gwas_label, lead_p_label) {
  dt <- copy(dt)
  dt[, proxy       := 1 - pmin(abs(pos2 - lead_pos) / half_win, 1)]
  dt[, neg_log10_p := -log10(pmax(pvalue, 1e-300))]
  dt[, is_lead     := (pos2 == lead_pos)]

  ggplot(dt[is_lead == FALSE],
         aes(x = pos2 / 1e6, y = neg_log10_p, colour = proxy)) +
    geom_point(size = 1.5, alpha = 0.75) +
    geom_point(data = dt[is_lead == TRUE],
               aes(x = pos2 / 1e6, y = neg_log10_p),
               colour = "#CC0000", fill = "#CC0000",
               size = 4, shape = 23) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
               colour = "grey25", linewidth = 0.45) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dotted",
               colour = "grey50", linewidth = 0.4) +
    scale_colour_gradientn(colours = lz_colours, limits = c(0, 1),
                           name = "Distance\nproxy",
                           guide = guide_colourbar(barheight = 3)) +
    scale_x_continuous(labels = function(x) sprintf("%.3f", x),
                       expand = expansion(mult = 0.03)) +
    labs(title = paste0(gwas_label, "  (lead p = ", lead_p_label, ")"),
         x = paste0(chrom, " position (Mb)"),
         y = expression(-log[10](p))) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank(),
          legend.position = "right")
}

# ---- read regions -----------------------------------------------------------
bulk_path <- find_tbl(opt$bulk_dir, opt$cell_type)
bulk_reg  <- read_region(bulk_path, chrom, start, end)

sn_reg <- NULL
if (!is.null(sn_ct) && !is.na(sn_ct)) {
  sn_path <- find_tbl(opt$sn_dir, sn_ct)
  sn_reg  <- read_region(sn_path, chrom, start, end)
}

if (is.null(bulk_reg) || nrow(bulk_reg) == 0)
  stop("No bulk data found in region ", chrom, ":", start, "-", end)

# get lead SNP p-value from bulk
lead_row  <- bulk_reg[pos2 == pos]
bulk_lead_p <- if (nrow(lead_row) > 0) signif(lead_row$pvalue[1], 3) else "N/A"

p_bulk <- make_panel(bulk_reg, pos, HALF,
                     paste0("Bulk GWAS — ", opt$cell_type), bulk_lead_p)

# ---- combine panels ---------------------------------------------------------
if (!is.null(sn_reg) && nrow(sn_reg) > 0) {
  lead_sn   <- sn_reg[pos2 == pos]
  sn_lead_p <- if (nrow(lead_sn) > 0) signif(lead_sn$pvalue[1], 3) else "not found"
  p_sn <- make_panel(sn_reg, pos, HALF,
                     paste0("snRNA-seq GWAS — ", sn_ct), sn_lead_p)
  combined <- plot_grid(p_bulk, p_sn, ncol = 1, align = "v",
                        axis = "lr", rel_heights = c(1, 1))
  n_panels <- 2
} else {
  message("Lead SNP not found in sn GWAS — showing bulk only")
  combined  <- p_bulk
  n_panels  <- 1
  sn_lead_p <- "not found in sn"
}

title_grob <- ggdraw() +
  draw_label(
    paste0("Regional association plot  |  ", chrom, ":", format(pos, big.mark=","),
           "  (", opt$cell_type, ")\n",
           "Bulk p = ", bulk_lead_p,
           "   sn p = ", if (exists("sn_lead_p")) sn_lead_p else "N/A",
           "   Window: ±", opt$window_kb, " kb",
           "\nRed diamond = lead SNP  |  Colour = distance to lead (proxy for LD)"),
    fontface = "bold", size = 9, x = 0.02, hjust = 0, vjust = 1
  )

final <- plot_grid(title_grob, combined, ncol = 1,
                   rel_heights = c(0.12 * n_panels, 1))

safe_marker <- gsub(":", "_", opt$marker)
out_png <- file.path(opt$outdir,
  paste0(opt$cell_type, "_", sn_ct, "_", safe_marker, ".png"))
ggsave(out_png, final,
       width = 9, height = if (n_panels == 2) 8 else 5, dpi = 150)
message("Saved: ", out_png)
