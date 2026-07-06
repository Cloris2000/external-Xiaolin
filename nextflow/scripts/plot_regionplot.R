#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_regionplot.R
#
# Draw a locus-zoom region plot for a given cell type centered on a lead SNP,
# coloring points by direction of effect (positive β vs negative β).
#
# Uses topr::regionplot() with two datasets (positive / negative effect),
# so the color legend shows effect direction rather than cell type.
#
# Usage:
#   Rscript plot_regionplot.R \
#     --tbl_file    /path/to/CellType_meta_analysis_...tbl \
#     --cell_type   Oligodendrocyte \
#     --lead_snp    chr7:12252119:G:GT \
#     --region_size 1000000 \
#     --output_dir  /path/to/plots/regionplot \
#     --gw_thresh   5e-8
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(topr)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option("--tbl_file",    type = "character", default = NULL),
  make_option("--cell_type",   type = "character", default = NULL),
  make_option("--lead_snp",    type = "character", default = NULL,
              help = "Lead SNP MarkerName e.g. chr7:12252119:G:GT"),
  make_option("--region_size", type = "double",    default = 1e6,
              help = "Half-window in bp around lead SNP [default 1Mb total = ±500kb]"),
  make_option("--output_dir",  type = "character", default = NULL),
  make_option("--gw_thresh",   type = "double",    default = 5e-8)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$tbl_file))   stop("--tbl_file is required")
if (is.null(opt$cell_type))  stop("--cell_type is required")
if (is.null(opt$lead_snp))   stop("--lead_snp is required")
if (is.null(opt$output_dir)) stop("--output_dir is required")

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

CT       <- opt$cell_type
LEAD     <- opt$lead_snp
GW       <- opt$gw_thresh
HALF_WIN <- opt$region_size / 2

# ---- Parse lead SNP coordinates --------------------------------------------
lead_parts <- strsplit(LEAD, ":")[[1]]
lead_chr   <- lead_parts[1]                          # e.g. "chr7"
lead_pos   <- as.numeric(lead_parts[2])              # e.g. 12252119
chr_num    <- as.integer(gsub("^chr", "", lead_chr)) # e.g. 7

region_start <- lead_pos - HALF_WIN
region_end   <- lead_pos + HALF_WIN
cat(sprintf("[%s] Region: %s:%d-%d\n", CT, lead_chr, region_start, region_end))

# ---- Fast region extraction via grep ----------------------------------------
# Build grep pattern matching chrN:POS in the window
cat(sprintf("[%s] Extracting region from %s ...\n", CT, basename(opt$tbl_file)))

# Read header
header <- fread(opt$tbl_file, nrows = 0)
col_names <- names(header)

# Use fread with grep filter for speed (reads only matching lines)
# Pattern: chrN: followed by digits in range — use R range filter after reading chr
# Fastest: read full chrom grep then filter by position in R
tmp_cmd <- sprintf("grep '^%s:' '%s'", lead_chr, opt$tbl_file)
region_raw <- fread(cmd = tmp_cmd, header = FALSE, col.names = col_names,
                    showProgress = FALSE)

cat(sprintf("[%s] Chr %s variants loaded: %d\n", CT, lead_chr, nrow(region_raw)))

# ---- Parse and filter to window --------------------------------------------
df <- region_raw %>%
  tidyr::separate(MarkerName, into = c("CHR_", "POS_", "REF_", "ALT_"),
                  sep = ":", remove = FALSE, extra = "drop") %>%
  rename(ID = MarkerName, P = `P-value`, Effect = Effect) %>%
  mutate(
    CHROM = as.integer(gsub("^chr", "", CHR_)),
    POS   = as.numeric(POS_),
    P     = as.numeric(P),
    Effect = as.numeric(Effect)
  ) %>%
  filter(!is.na(POS), POS >= region_start, POS <= region_end,
         !is.na(P), P > 0, P <= 1) %>%
  dplyr::select(ID, CHROM, POS, P, Effect)

cat(sprintf("[%s] Variants in region: %d\n", CT, nrow(df)))

if (nrow(df) == 0) stop("No variants in region after filtering.")

# ---- Genomic inflation in this region (info only) --------------------------
chi2   <- qchisq(1 - df$P, 1)
lambda <- round(median(chi2, na.rm = TRUE) / qchisq(0.5, 1), 4)
n_gw   <- sum(df$P < GW, na.rm = TRUE)
cat(sprintf("[%s] λ=%.4f  GW hits in region=%d\n", CT, lambda, n_gw))

# ---- Split by effect direction ---------------------------------------------
# topr::regionplot requires the lead variant to be in the FIRST dataset.
# Determine which group contains the lead SNP and put it first.
df_pos <- df %>% filter(Effect >= 0) %>% dplyr::select(CHROM, POS, P, ID, Effect)
df_neg <- df %>% filter(Effect <  0) %>% dplyr::select(CHROM, POS, P, ID, Effect)
cat(sprintf("[%s] Positive β: %d variants | Negative β: %d variants\n",
            CT, nrow(df_pos), nrow(df_neg)))

lead_in_pos <- LEAD %in% df_pos$ID
if (lead_in_pos) {
  datasets      <- list(df_pos, df_neg)
  colors        <- c("#D73027", "#4575B4")  # red first = positive β
  leg_labels    <- c("Positive effect (beta > 0)", "Negative effect (beta < 0)")
} else {
  datasets      <- list(df_neg, df_pos)
  colors        <- c("#4575B4", "#D73027")  # blue first = negative β
  leg_labels    <- c("Negative effect (beta < 0)", "Positive effect (beta > 0)")
}
cat(sprintf("[%s] Lead SNP in %s group → placed first in dataset list\n",
            CT, if (lead_in_pos) "positive" else "negative"))

# ---- Build region plot with topr -------------------------------------------
cat(sprintf("[%s] Generating region plot...\n", CT))

out_file <- file.path(opt$output_dir,
                      sprintf("%s_%s_regionplot.png", CT, gsub(":", "_", LEAD)))

plots <- regionplot(
  datasets,
  variant           = LEAD,
  region_size       = opt$region_size,
  color             = colors,
  legend_labels     = leg_labels,
  legend_name       = "Effect direction",
  sign_thresh       = GW,
  sign_thresh_color = "black",
  ntop              = 5,
  title             = sprintf("%s - locus %s:%s-%sMb  (lambda=%.4f, %d GW-sig variants)",
                              CT, lead_chr,
                              format(round(region_start/1e6, 2), nsmall=2),
                              format(round(region_end/1e6,   2), nsmall=2),
                              lambda, n_gw),
  build             = 38,
  show_overview     = TRUE,
  extract_plots     = TRUE   # return list(main_plot, overview_plot, gene_plot)
)

# Assemble with cowplot and save via ggsave (ggplot objects, not grid)
combined <- cowplot::plot_grid(
  plots$main_plot,
  plots$overview_plot,
  plots$gene_plot,
  ncol    = 1,
  rel_heights = c(7, 1.25, 2)
)

ggsave(out_file, combined, width = 14, height = 9, dpi = 300, limitsize = FALSE)
cat(sprintf("[%s] Saved: %s\n", CT, out_file))
