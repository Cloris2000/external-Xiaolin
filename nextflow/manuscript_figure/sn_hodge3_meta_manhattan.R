#!/usr/bin/env Rscript
# =============================================================================
# Standard per-cell-type Manhattan plots for sn Hodge 3-cohort meta-analysis
# using topr::manhattanExtra() with nearest-gene annotation.
#
# Inputs:  results/meta_analysis_sn_hodge3/*_meta_analysis_*.tbl  (final only)
# Outputs:
#   manuscript_figure/sn_hodge3_meta_manhattan/{CellType}_manhattan.png  (x19)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(topr)
})

META_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/meta_analysis_sn_hodge3"
OUT_DIR  <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"
PLOT_DIR <- file.path(OUT_DIR, "sn_hodge3_meta_manhattan")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

PVAL_GWS  <- 5e-8
PVAL_SUGG <- 1e-5
# Label peaks at suggestive threshold (few GWS hits in this meta)
ANNOTATE  <- 1e-5

tbl_files <- list.files(META_DIR, pattern = "_meta_analysis_.*\\.tbl$", full.names = TRUE)
tbl_files <- tbl_files[!grepl("1\\.tbl$|\\.tbl\\.info$", tbl_files)]
if (!length(tbl_files)) stop("No METAL .tbl files in ", META_DIR)

# Stable cell-type order (Hodge subclasses)
ct_order <- c(
  "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.ET", "L5.6.IT.Car3",
  "L5.6.NP", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC",
  "Oligodendrocyte", "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC"
)
names(tbl_files) <- sub("_meta_analysis_.*$", "", basename(tbl_files))
tbl_files <- tbl_files[intersect(ct_order, names(tbl_files))]
cat("Plotting", length(tbl_files), "cell types with topr::manhattanExtra()\n")

read_meta_for_topr <- function(f) {
  dt <- fread(
    f,
    select = c("MarkerName", "P-value"),
    showProgress = FALSE
  )
  setnames(dt, c("snp", "P"))
  parts <- tstrsplit(dt$snp, ":", fixed = TRUE)
  dt[, CHROM := as.integer(sub("^[Cc][Hh][Rr]", "", parts[[1L]]))]
  dt[, POS := as.integer(parts[[2L]])]
  dt <- dt[CHROM %in% 1:22 & !is.na(POS) & !is.na(P) & P > 0 & P <= 1]
  dt <- dt[!grepl("*", snp, fixed = TRUE)]
  # topr prefers plain data.frame with CHROM/POS/P
  as.data.frame(dt[, .(CHROM, POS, P)])
}

make_one_manhattan <- function(df, cell_type) {
  # ymax with a little headroom for labels
  max_nlp <- max(-log10(df$P), na.rm = TRUE)
  ymax <- max(10, ceiling(max_nlp) + 1.5)

  p <- manhattanExtra(
    df,
    genome_wide_thresh = PVAL_GWS,
    suggestive_thresh  = PVAL_SUGG,
    flank_size         = 1e6,
    region_size        = 1e6,
    annotate           = ANNOTATE,
    sign_thresh        = c(PVAL_SUGG, PVAL_GWS),
    ymax               = ymax,
    show_legend        = FALSE
  )
  # topr returns ggplot; add cell-type title
  p + labs(title = cell_type) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 7)
    )
}

for (i in seq_along(tbl_files)) {
  ct <- names(tbl_files)[i]
  cat(sprintf("[%2d/%2d] %s\n", i, length(tbl_files), ct))
  df <- read_meta_for_topr(tbl_files[[i]])
  cat("  variants:", nrow(df),
      "  minP=", format(min(df$P), scientific = TRUE), "\n")
  p <- tryCatch(
    make_one_manhattan(df, ct),
    error = function(e) {
      cat("  manhattanExtra failed:", conditionMessage(e),
          " — falling back to topr::manhattan()\n")
      topr::manhattan(
        df, annotate = ANNOTATE, region_size = 1e6,
        sign_thresh = c(PVAL_SUGG, PVAL_GWS),
        show_legend = FALSE
      ) + labs(title = ct) +
        theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
    }
  )

  out_png <- file.path(PLOT_DIR, paste0(ct, "_manhattan.png"))
  ggsave(out_png, p, width = 10, height = 3.2, dpi = 200, bg = "white")
  cat("  wrote", out_png, "\n")
  rm(df, p); gc(verbose = FALSE)
}

cat("Done. 19 PNGs in:", PLOT_DIR, "\n")
