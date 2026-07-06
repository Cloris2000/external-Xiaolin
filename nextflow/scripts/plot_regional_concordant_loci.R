#!/usr/bin/env Rscript
# Regional (LocusZoom-style) plots for every concordant independent locus.
#
# For each row in independent_loci_concordant.tsv, extracts the ±500 kb
# region from both the bulk and sn .tbl files and draws a two-panel
# regional association plot (bulk on top, sn on bottom).
# SNPs are coloured by their r² proxy to the lead SNP using physical
# distance as a proxy (no LD panel needed):
#   red   = lead SNP itself
#   gradient blue→orange = distance from lead (closer = warmer)
#
# Outputs: one PNG per locus saved to <outdir>/regional_plots/
#
# Usage:
#   Rscript plot_regional_concordant_loci.R
#   Rscript plot_regional_concordant_loci.R \
#     --manifest /path/to/independent_loci_concordant.tsv \
#     --bulk_dir  /path/to/bulk_meta \
#     --sn_dir    /path/to/sn_meta   \
#     --outdir    /path/to/top_hits

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# ---- options ----------------------------------------------------------------
option_list <- list(
  make_option("--manifest", type = "character",
    default = paste0("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/",
                     "results/sn_bulk_meta_similarity_design_matrix/top_hits/",
                     "independent_loci_concordant.tsv")),
  make_option("--bulk_dir", type = "character",
    default = paste0("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/",
                     "results/meta_analysis_13cohorts_design_matrix")),
  make_option("--sn_dir", type = "character",
    default = paste0("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/",
                     "results/meta_analysis_sn_rosmap_msbb_hbcc_alias15")),
  make_option("--outdir", type = "character",
    default = paste0("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/",
                     "results/sn_bulk_meta_similarity_design_matrix/top_hits")),
  make_option("--window_kb", type = "integer", default = 500,
    help = "Half-window around lead SNP in kb [default %default]"),
  make_option("--max_loci", type = "integer", default = 200,
    help = "Safety cap: process at most this many loci [default %default]")
)
opt  <- parse_args(OptionParser(option_list = option_list))
HALF <- opt$window_kb * 1000L
out_regional <- file.path(opt$outdir, "regional_plots")
dir.create(out_regional, recursive = TRUE, showWarnings = FALSE)

# ---- load manifest ----------------------------------------------------------
if (!file.exists(opt$manifest)) {
  stop("Manifest not found: ", opt$manifest,
       "\nRun bulk_top_hits_sn_direction.R first.")
}
manifest <- fread(opt$manifest, sep = "\t")
message("Manifest loaded: ", nrow(manifest), " concordant loci")
if (nrow(manifest) > opt$max_loci) {
  message("  Capping to top ", opt$max_loci, " loci (by bulk_p)")
  manifest <- manifest[order(bulk_p)][seq_len(opt$max_loci)]
}

# ---- helpers ----------------------------------------------------------------
find_tbl <- function(dir, cell_type) {
  pat <- paste0("^", gsub("([.])", "\\\\.", cell_type),
                "_meta_analysis_.*\\.tbl$")
  f <- list.files(dir, pattern = pat, full.names = TRUE)
  f <- f[!grepl("1\\.tbl$", basename(f))]
  if (length(f) == 0) return(NA_character_)
  f[1]
}

# Read one chromosome region from a tbl file using fread's select + filter
read_region <- function(path, chrom, start, end) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  dt <- fread(path, sep = "\t", data.table = TRUE, verbose = FALSE)
  setnames(dt, tolower(gsub("[^A-Za-z0-9]", "_", names(dt))))
  col_map <- list(marker = "markername", allele1 = "allele1",
                  allele2 = "allele2", effect = "effect",
                  stderr = "stderr", pvalue = "p_value")
  for (nm in names(col_map)) {
    old <- col_map[[nm]]
    if (old %in% names(dt) && !(nm %in% names(dt))) setnames(dt, old, nm)
  }
  needed <- c("marker", "allele1", "allele2", "effect", "stderr", "pvalue")
  if (!all(needed %in% names(dt))) return(NULL)
  dt <- dt[, needed, with = FALSE]
  dt <- dt[nchar(allele1) == 1 & nchar(allele2) == 1]
  parts  <- tstrsplit(dt$marker, ":", fixed = TRUE)
  dt[, chrom_col := parts[[1]]]
  dt[, pos       := as.integer(parts[[2]])]
  dt[, effect    := as.numeric(effect)]
  dt[, stderr    := as.numeric(stderr)]
  dt[, pvalue    := as.numeric(pvalue)]
  dt[chrom_col == chrom & pos >= start & pos <= end]
}

# Colour SNPs by distance to lead (used as proxy for LD in absence of panel)
# Returns a numeric 0–1 score: 1 = lead SNP, 0 = far edge of window
distance_proxy <- function(pos_vec, lead_pos, half_win) {
  d <- abs(pos_vec - lead_pos)
  1 - pmin(d / half_win, 1)
}

# Colour palette mimicking LocusZoom LD colouring (blue → cyan → green →
# yellow → orange → red), but driven by distance instead of r²
lz_colours <- colorRampPalette(
  c("#357EBD", "#46B8DA", "#5CB85C", "#F0AD4E", "#D9534F")
)(100)

# ---- main loop --------------------------------------------------------------
n_ok <- 0L

# Cache open tbl file handles to avoid re-reading the same file repeatedly
# (each cell type may appear multiple times across loci)
tbl_cache <- list()

get_region_cached <- function(dir, cell_type, chrom, start, end, cache_env) {
  key <- paste0(dir, "|", cell_type)
  if (!exists(key, envir = cache_env)) {
    path <- find_tbl(dir, cell_type)
    if (is.na(path)) {
      assign(key, NULL, envir = cache_env)
    } else {
      message("    caching full tbl for ", cell_type, " ...")
      dt_full <- fread(path, sep = "\t", data.table = TRUE, verbose = FALSE)
      setnames(dt_full, tolower(gsub("[^A-Za-z0-9]", "_", names(dt_full))))
      col_map <- list(marker = "markername", allele1 = "allele1",
                      allele2 = "allele2", effect = "effect",
                      stderr = "stderr", pvalue = "p_value")
      for (nm in names(col_map)) {
        old <- col_map[[nm]]
        if (old %in% names(dt_full) && !(nm %in% names(dt_full)))
          setnames(dt_full, old, nm)
      }
      needed <- c("marker", "allele1", "allele2", "effect", "stderr", "pvalue")
      if (!all(needed %in% names(dt_full))) {
        assign(key, NULL, envir = cache_env)
      } else {
        dt_full <- dt_full[, needed, with = FALSE]
        dt_full <- dt_full[nchar(allele1) == 1 & nchar(allele2) == 1]
        parts   <- tstrsplit(dt_full$marker, ":", fixed = TRUE)
        dt_full[, chrom_col := parts[[1]]]
        dt_full[, pos       := as.integer(parts[[2]])]
        dt_full[, effect    := as.numeric(effect)]
        dt_full[, stderr    := as.numeric(stderr)]
        dt_full[, pvalue    := as.numeric(pvalue)]
        setkeyv(dt_full, c("chrom_col", "pos"))
        assign(key, dt_full, envir = cache_env)
      }
    }
  }
  full <- get(key, envir = cache_env)
  if (is.null(full)) return(NULL)
  full[chrom_col == chrom & pos >= start & pos <= end]
}

cache_env <- new.env(parent = emptyenv())

for (k in seq_len(nrow(manifest))) {
  row <- manifest[k]

  bulk_ct <- row$bulk_cell_type
  sn_ct   <- row$sn_cell_type
  chrom   <- row$chrom
  lead_p  <- row$pos
  start   <- pmax(0L, lead_p - HALF)
  end     <- lead_p + HALF
  label   <- row$locus_label

  safe_label <- gsub("[^A-Za-z0-9_]", "_", label)
  out_png <- file.path(out_regional,
    sprintf("%s_%s_%s.png", bulk_ct, sn_ct, safe_label))

  if (file.exists(out_png)) {
    message("  [skip] ", basename(out_png), " already exists")
    n_ok <- n_ok + 1L
    next
  }

  message("\n[", k, "/", nrow(manifest), "] ",
          bulk_ct, " / ", sn_ct, "  — ", chrom, ":", start, "-", end)

  bulk_reg <- get_region_cached(opt$bulk_dir, bulk_ct, chrom, start, end, cache_env)
  sn_reg   <- get_region_cached(opt$sn_dir,   sn_ct,   chrom, start, end, cache_env)

  if (is.null(bulk_reg) || nrow(bulk_reg) == 0) {
    message("  skipping – no bulk data in region"); next
  }

  # compute distance-proxy colour score
  bulk_reg[, proxy := distance_proxy(pos, lead_p, HALF)]
  bulk_reg[, neg_log10_p := -log10(pmax(pvalue, 1e-300))]
  bulk_reg[, is_lead := (pos == lead_p)]

  make_panel <- function(dt, gwas_label, lead_pos, chrom_str, lead_bulk_p, lead_sn_p) {
    dt <- copy(dt)
    dt[, proxy        := distance_proxy(pos, lead_pos, HALF)]
    dt[, neg_log10_p  := -log10(pmax(pvalue, 1e-300))]
    dt[, is_lead      := (pos == lead_pos)]
    dt[, colour_score := proxy]

    # non-lead points coloured by distance
    p <- ggplot(dt[is_lead == FALSE],
                aes(x = pos / 1e6, y = neg_log10_p, colour = colour_score)) +
      geom_point(size = 1.4, alpha = 0.75) +
      # lead SNP on top in red
      geom_point(data = dt[is_lead == TRUE],
                 aes(x = pos / 1e6, y = neg_log10_p),
                 colour = "#CC0000", size = 3.5, shape = 23,
                 fill = "#CC0000") +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
                 colour = "grey30", linewidth = 0.4) +
      geom_hline(yintercept = -log10(1e-5), linetype = "dotted",
                 colour = "grey50", linewidth = 0.4) +
      scale_colour_gradientn(
        colours = lz_colours,
        limits  = c(0, 1),
        name    = "Distance\nproxy\n(1=lead)",
        guide   = guide_colourbar(barheight = 4)
      ) +
      scale_x_continuous(
        labels = function(x) sprintf("%.2f", x),
        expand = expansion(mult = 0.03)
      ) +
      labs(
        title = gwas_label,
        x     = paste0(chrom_str, " position (Mb)"),
        y     = expression(-log[10](p))
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title       = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        legend.position  = "right"
      )
    p
  }

  p_bulk <- make_panel(bulk_reg,
    gwas_label  = paste0("Bulk GWAS — ", bulk_ct,
                         "  (lead p = ", signif(row$bulk_p, 3), ")"),
    lead_pos    = lead_p,
    chrom_str   = chrom,
    lead_bulk_p = row$bulk_p,
    lead_sn_p   = row$sn_p)

  if (!is.null(sn_reg) && nrow(sn_reg) > 0) {
    p_sn <- make_panel(sn_reg,
      gwas_label  = paste0("snRNA-seq GWAS — ", sn_ct,
                           "  (lead p = ", signif(row$sn_p, 3), ")"),
      lead_pos    = lead_p,
      chrom_str   = chrom,
      lead_bulk_p = row$bulk_p,
      lead_sn_p   = row$sn_p)
    # stack bulk (top) and sn (bottom)
    combined <- cowplot::plot_grid(p_bulk, p_sn,
                                   ncol  = 1,
                                   align = "v",
                                   axis  = "lr",
                                   rel_heights = c(1, 1))
    title_grob <- cowplot::ggdraw() +
      cowplot::draw_label(
        paste0("Regional plot — ", label,
               "  |  bulk: ", bulk_ct, "  /  sn: ", sn_ct,
               "\nbulk effect = ", round(row$bulk_effect, 3),
               " (SE ", round(row$bulk_stderr, 3), ")",
               "   sn effect = ", round(row$sn_effect, 3),
               " (SE ", round(row$sn_stderr, 3), ")",
               "   direction: concordant"),
        fontface = "bold", size = 9, x = 0.02, hjust = 0)
    final_plot <- cowplot::plot_grid(title_grob, combined,
                                     ncol = 1,
                                     rel_heights = c(0.08, 1))
  } else {
    message("  no sn data in region — plotting bulk only")
    final_plot <- p_bulk +
      labs(title = paste0("Bulk GWAS — ", bulk_ct,
                          "  (lead p = ", signif(row$bulk_p, 3), ")",
                          "\n[sn data not available in region]"))
  }

  tryCatch(
    ggsave(out_png, final_plot, width = 9, height = if (!is.null(sn_reg) && nrow(sn_reg) > 0) 8 else 5,
           dpi = 150),
    error = function(e) message("  ERROR saving plot: ", e$message)
  )
  message("  saved: ", basename(out_png))
  n_ok <- n_ok + 1L
}

message("\n=== Done ===  ", n_ok, " regional plots saved to: ", out_regional)
