#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_meta_results.R
#
# Three visualizations for the 12-cohort multi-ancestry GWAS meta-analysis:
#
#   1. Manhattan + QQ plots   — per cell type, with λ annotation
#   2. I² heterogeneity histogram — per cell type + 19-panel overview
#   3. Cross-cell-type heatmap   — -log10(P) at p < 1e-5 loci
#
# Usage (single cell type):
#   Rscript plot_meta_results.R \
#     --results_dir  /path/to/meta_analysis_all_cohorts \
#     --output_dir   /path/to/plots \
#     --cell_type    Astrocyte          # optional; omit to process all
#     --p_thresh     1e-5               # locus threshold for heatmap
#     --gw_thresh    5e-8               # genome-wide significance line
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(pheatmap)
  library(qqman)
  library(cowplot)
})

# ---- CLI args ---------------------------------------------------------------
option_list <- list(
  make_option("--results_dir", type = "character", default = NULL),
  make_option("--output_dir",  type = "character", default = NULL),
  make_option("--cell_type",   type = "character", default = NULL,
              help = "Single cell type to process (default: all)"),
  make_option("--p_thresh",    type = "double",    default = 1e-5,
              help = "P-value threshold for heatmap loci [default: 1e-5]"),
  make_option("--gw_thresh",   type = "double",    default = 5e-8,
              help = "Genome-wide significance threshold [default: 5e-8]"),
  make_option("--min_r2",      type = "double",    default = 1e2,
              help = "Minimum r^2 for LD clumping (unused placeholder)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$results_dir)) stop("--results_dir is required")
if (is.null(opt$output_dir))  stop("--output_dir is required")

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$output_dir, "manhattan"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$output_dir, "het"),       showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$output_dir, "heatmap"),   showWarnings = FALSE, recursive = TRUE)

GW   <- opt$gw_thresh
SUGG <- opt$p_thresh

# ---- Discover .tbl files ----------------------------------------------------
tbl_files <- list.files(opt$results_dir, pattern = "\\.tbl$", full.names = TRUE)
# Exclude stale *1.tbl files (previous partial runs)
tbl_files <- tbl_files[!grepl("1\\.tbl$", tbl_files)]

if (length(tbl_files) == 0) stop("No .tbl files found in: ", opt$results_dir)

# Extract cell type from filename (everything before _meta_analysis_)
get_cell_type <- function(path) {
  bn <- basename(path)
  sub("_meta_analysis_.*\\.tbl$", "", bn)
}

cell_types_all <- sapply(tbl_files, get_cell_type)
names(tbl_files) <- cell_types_all

if (!is.null(opt$cell_type)) {
  if (!(opt$cell_type %in% cell_types_all))
    stop("Cell type '", opt$cell_type, "' not found. Available: ",
         paste(sort(cell_types_all), collapse = ", "))
  tbl_files  <- tbl_files[opt$cell_type]
  cell_types <- opt$cell_type
} else {
  cell_types <- sort(cell_types_all)
}

cat(sprintf("Processing %d cell type(s): %s\n", length(cell_types),
            paste(cell_types, collapse = ", ")))

# ---- Helper: read & parse one .tbl ------------------------------------------
read_tbl <- function(path, cell_type) {
  cat(sprintf("  Reading %s ...\n", basename(path)))
  dt <- fread(path, header = TRUE, showProgress = FALSE) %>% as.data.frame()
  dt %>%
    tidyr::separate(MarkerName, into = c("CHROM", "POS", "REF", "ALT"),
                    sep = ":", remove = FALSE, extra = "drop") %>%
    rename(ID = MarkerName, P = `P-value`) %>%
    mutate(
      cell_type = cell_type,
      CHROM = suppressWarnings(as.integer(gsub("^chr", "", CHROM))),
      POS   = as.numeric(POS),
      P     = as.numeric(P),
      HetISq   = if ("HetISq"   %in% names(.)) as.numeric(HetISq)   else NA_real_,
      HetPVal  = if ("HetPVal"  %in% names(.)) as.numeric(HetPVal)  else NA_real_
    ) %>%
    filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1)
}

# ============================================================================
# 1. MANHATTAN + QQ  (per cell type)
# ============================================================================
make_manhattan_qq <- function(df, cell_type, out_dir) {

  # -- genomic inflation
  chi2   <- qchisq(1 - df$P, 1)
  lambda <- round(median(chi2, na.rm = TRUE) / qchisq(0.5, 1), 4)

  qqman_df <- df %>%
    filter(!is.na(CHROM)) %>%
    rename(CHR = CHROM, BP = POS, SNP = ID)

  # Manhattan
  man_file <- file.path(out_dir, "manhattan",
                        paste0(cell_type, "_manhattan.png"))
  png(man_file, width = 5400, height = 1800, res = 300)
  qqman::manhattan(
    qqman_df,
    suggestiveline = -log10(SUGG),
    genomewideline = -log10(GW),
    col            = c("#2171B5", "#6BAED6"),
    cex            = 0.45,
    cex.axis       = 0.8,
    main           = paste("Meta-Analysis Manhattan:", cell_type),
    sub            = sprintf("λ = %.4f  |  GW sig (p<%.0e): %d  |  Suggestive (p<%.0e): %d",
                             lambda,
                             GW, sum(df$P < GW,   na.rm = TRUE),
                             SUGG, sum(df$P < SUGG, na.rm = TRUE))
  )
  dev.off()

  # QQ
  qq_file <- file.path(out_dir, "manhattan",
                       paste0(cell_type, "_qq.png"))
  png(qq_file, width = 1600, height = 1600, res = 300)
  qqman::qq(
    df$P,
    main = paste("QQ Plot:", cell_type),
    sub  = sprintf("λ = %.4f", lambda),
    cex  = 0.5,
    col  = "#2171B5"
  )
  dev.off()

  list(lambda = lambda, n_gw = sum(df$P < GW, na.rm = TRUE),
       n_sugg = sum(df$P < SUGG, na.rm = TRUE))
}

# ============================================================================
# 2. I² HETEROGENEITY HISTOGRAM  (per cell type + combined overview)
# ============================================================================
make_het_histogram <- function(df, cell_type, out_dir) {

  if (all(is.na(df$HetISq))) {
    cat(sprintf("  [%s] No HetISq data — skipping heterogeneity plot\n", cell_type))
    return(NULL)
  }

  het_df <- df %>% filter(!is.na(HetISq))
  pct50 <- round(100 * mean(het_df$HetISq > 50, na.rm = TRUE), 1)
  pct75 <- round(100 * mean(het_df$HetISq > 75, na.rm = TRUE), 1)
  med   <- round(median(het_df$HetISq, na.rm = TRUE), 1)

  p <- ggplot(het_df, aes(x = HetISq)) +
    geom_histogram(binwidth = 2, fill = "#4292C6", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 50, color = "#E6550D", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 75, color = "#A63603", linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = 52, y = Inf, label = "I²=50%", color = "#E6550D",
             hjust = 0, vjust = 1.5, size = 3) +
    annotate("text", x = 77, y = Inf, label = "I²=75%", color = "#A63603",
             hjust = 0, vjust = 1.5, size = 3) +
    labs(
      title    = paste("Heterogeneity (I²) distribution:", cell_type),
      subtitle = sprintf("Median I²=%.1f%%  |  >50%%: %.1f%%  |  >75%%: %.1f%%  |  n=%s",
                         med, pct50, pct75, format(nrow(het_df), big.mark = ",")),
      x = "I² (%)",
      y = "Number of variants"
    ) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))

  out_file <- file.path(out_dir, "het",
                        paste0(cell_type, "_het_histogram.png"))
  ggsave(out_file, p, width = 7, height = 4, dpi = 300)

  list(cell_type = cell_type, median_ISq = med, pct_gt50 = pct50, pct_gt75 = pct75,
       n_variants = nrow(het_df))
}

# ============================================================================
# 3. CROSS-CELL-TYPE HEATMAP
# ============================================================================
# Collect top loci from all cell types, then plot -log10(P) matrix

build_heatmap_matrix <- function(all_top, cell_types_used) {

  # Cluster nearby variants into loci (1 Mb window)
  all_top <- all_top %>%
    arrange(CHROM, POS) %>%
    mutate(locus = NA_character_)

  loci   <- list()
  prev_c <- -1L; prev_p <- -1L; locus_id <- 0L
  for (i in seq_len(nrow(all_top))) {
    c <- all_top$CHROM[i]; p <- all_top$POS[i]
    if (c != prev_c || (p - prev_p) > 1e6) {
      locus_id <- locus_id + 1L
      all_top$locus[i] <- sprintf("chr%d:%d", c, round(p / 1e6))
    } else {
      all_top$locus[i] <- all_top$locus[i - 1]
    }
    prev_c <- c; prev_p <- p
  }

  # For each (locus, cell_type), take minimum P
  mat_long <- all_top %>%
    group_by(locus, cell_type) %>%
    summarise(min_P = min(P, na.rm = TRUE), .groups = "drop")

  # Pivot to wide matrix
  mat <- mat_long %>%
    pivot_wider(names_from = cell_type, values_from = min_P, values_fill = 1) %>%
    column_to_rownames("locus")

  # Keep columns in consistent order
  mat <- mat[, intersect(cell_types_used, colnames(mat)), drop = FALSE]

  # Convert to -log10(P)
  log_mat <- -log10(as.matrix(mat))
  log_mat[is.infinite(log_mat)] <- 0
  log_mat
}

make_cross_celltype_heatmap <- function(all_data, out_dir) {

  cat(sprintf("  Building cross-cell-type heatmap (p < %.0e)...\n", SUGG))

  all_top <- all_data %>%
    filter(P < SUGG, !is.na(CHROM))

  if (nrow(all_top) == 0) {
    cat("  No variants pass threshold — skipping heatmap\n")
    return(invisible(NULL))
  }

  cat(sprintf("  %d variant×cell-type rows pass threshold\n", nrow(all_top)))

  cell_types_used <- sort(unique(all_top$cell_type))
  log_mat <- build_heatmap_matrix(all_top, cell_types_used)

  n_loci <- nrow(log_mat)
  cat(sprintf("  %d loci in heatmap\n", n_loci))

  # Determine height: scale with number of loci, floor at 5 inches
  h <- max(5, min(0.25 * n_loci + 3, 40))

  breaks <- seq(0, max(log_mat, na.rm = TRUE) + 0.1, length.out = 101)
  colors <- colorRampPalette(c("white", "#C6DBEF", "#2171B5", "#084594", "#08306B"))(100)

  # Annotation: number of cell types with p < gw_thresh per locus
  row_ann <- data.frame(
    GW_sig_cell_types = rowSums(log_mat >= -log10(GW)),
    row.names = rownames(log_mat)
  )
  ann_colors <- list(
    GW_sig_cell_types = colorRampPalette(c("white", "#FD8D3C", "#D94701"))(
      max(row_ann$GW_sig_cell_types) + 1)
  )

  out_file <- file.path(out_dir, "heatmap", "cross_celltype_heatmap.png")
  png(out_file, width = max(2400, 200 * length(cell_types_used)),
      height = max(1800, 40 * n_loci + 400), res = 200)
  pheatmap(
    log_mat,
    color            = colors,
    breaks           = breaks,
    cluster_rows     = n_loci > 2,
    cluster_cols     = TRUE,
    annotation_row   = row_ann,
    annotation_colors = ann_colors,
    fontsize_row     = max(4, min(9, 200 / n_loci)),
    fontsize_col     = 9,
    border_color     = NA,
    main             = sprintf(
      "Cross-cell-type meta-analysis signal  (−log₁₀P, threshold p<%.0e)", SUGG),
    legend_breaks    = c(0, 3, 5, 7.3, max(log_mat, na.rm = TRUE)),
    legend_labels    = c("0", "p<1e-3", "p<1e-5", "p<5e-8",
                         sprintf("%.1f", max(log_mat, na.rm = TRUE)))
  )
  dev.off()
  cat(sprintf("  Saved heatmap: %s\n", out_file))

  # Also save the underlying table
  tbl_out <- as.data.frame(log_mat) %>%
    tibble::rownames_to_column("locus") %>%
    arrange(desc(rowMeans(across(where(is.numeric)))))
  write.table(tbl_out,
              file = file.path(out_dir, "heatmap", "cross_celltype_matrix.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ============================================================================
# MAIN — iterate over cell types
# ============================================================================
summary_rows    <- list()
het_stats_list  <- list()
all_data_list   <- list()

for (ct in cell_types) {
  cat(sprintf("\n=== %s ===\n", ct))
  path <- tbl_files[ct]
  df   <- tryCatch(read_tbl(path, ct), error = function(e) {
    cat(sprintf("  ERROR reading %s: %s\n", ct, e$message)); NULL
  })
  if (is.null(df)) next

  # -- Manhattan + QQ
  stats <- make_manhattan_qq(df, ct, opt$output_dir)
  summary_rows[[ct]] <- c(cell_type = ct,
                          n_variants  = nrow(df),
                          lambda      = stats$lambda,
                          n_gw        = stats$n_gw,
                          n_sugg      = stats$n_sugg)

  # -- Heterogeneity histogram
  het_info <- make_het_histogram(df, ct, opt$output_dir)
  if (!is.null(het_info)) het_stats_list[[ct]] <- het_info

  # -- Accumulate top hits for heatmap
  all_data_list[[ct]] <- df %>%
    filter(P < SUGG) %>%
    dplyr::select(cell_type, CHROM, POS, P)

  # Explicit cleanup for very large tables
  rm(df)
  invisible(gc(verbose = FALSE))
}

# ---- Combined I² overview (lightweight summary panel) -------------------
# Only run cross-cell-type summaries in all-cell-type mode.
if (is.null(opt$cell_type) && length(het_stats_list) > 1) {
  cat("\nGenerating combined heterogeneity overview panel...\n")
  het_summary <- bind_rows(lapply(het_stats_list, as.data.frame))
  het_summary$cell_type <- factor(het_summary$cell_type, levels = het_summary$cell_type)

  p1 <- ggplot(het_summary, aes(x = reorder(cell_type, median_ISq), y = median_ISq)) +
    geom_col(fill = "#4292C6") +
    coord_flip() +
    labs(title = "Median I² by Cell Type", x = "Cell Type", y = "Median I² (%)") +
    theme_bw(base_size = 10)

  p2 <- ggplot(het_summary, aes(x = reorder(cell_type, pct_gt50), y = pct_gt50)) +
    geom_col(fill = "#E6550D") +
    coord_flip() +
    labs(title = "Percent Variants with I² > 50%", x = "Cell Type", y = "Percent (%)") +
    theme_bw(base_size = 10)

  p3 <- ggplot(het_summary, aes(x = reorder(cell_type, pct_gt75), y = pct_gt75)) +
    geom_col(fill = "#A63603") +
    coord_flip() +
    labs(title = "Percent Variants with I² > 75%", x = "Cell Type", y = "Percent (%)") +
    theme_bw(base_size = 10)

  combined_het <- (p1 | p2 | p3) +
    plot_annotation(
      title = "I² Heterogeneity Summary Across Cell Types",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  tryCatch({
    ggsave(
      file.path(opt$output_dir, "het", "all_celltypes_het_overview.png"),
      combined_het,
      width = 18, height = 8, dpi = 220, limitsize = FALSE
    )
  }, error = function(e) {
    cat("  WARNING: combined het overview ggsave failed (patchwork/ggplot2 guide mismatch):", e$message, "\n")
    cat("  Saving panels individually instead...\n")
    ggsave(file.path(opt$output_dir, "het", "all_celltypes_median_I2.png"),    p1, width=7, height=6, dpi=220)
    ggsave(file.path(opt$output_dir, "het", "all_celltypes_pct_gt50.png"),     p2, width=7, height=6, dpi=220)
    ggsave(file.path(opt$output_dir, "het", "all_celltypes_pct_gt75.png"),     p3, width=7, height=6, dpi=220)
    cat("  Saved 3 individual panels.\n")
  })

  write.table(
    het_summary,
    file = file.path(opt$output_dir, "het", "all_celltypes_het_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("  Saved combined heterogeneity overview + summary table.\n")
}

# ---- Cross-cell-type heatmap -------------------------------------------
# Only run cross-cell-type summaries in all-cell-type mode.
if (is.null(opt$cell_type) && length(all_data_list) > 0) {
  all_data <- bind_rows(all_data_list)
  make_cross_celltype_heatmap(all_data, opt$output_dir)
}

# ---- Summary table ---------------------------------------------------------
if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, lapply(summary_rows, function(x) as.data.frame(t(x)))) %>%
    mutate(across(c(n_variants, n_gw, n_sugg), as.integer),
           lambda = as.numeric(lambda))
  out_sum <- file.path(opt$output_dir, "meta_analysis_summary.tsv")
  write.table(summary_df, out_sum, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n=== Summary across all cell types ===\n")
  print(summary_df)
  cat(sprintf("\nSummary table saved to: %s\n", out_sum))
}

cat("\nAll plots complete.\n")
cat(sprintf("Output directory: %s\n", opt$output_dir))
