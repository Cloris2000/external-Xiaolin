#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_meta_annotated.R
#
# Generate annotated Manhattan plots (gene labels at significant loci) and
# QQ plots for one cell type from the 15-cohort meta-analysis.
# Uses topr::manhattanExtra() for automatic nearest-gene annotation.
#
# Usage (called once per cell type, e.g. via SLURM array):
#   Rscript plot_meta_annotated.R \
#     --results_dir /path/to/meta_analysis_15cohorts \
#     --output_dir  /path/to/plots/annotated \
#     --cell_type   Astrocyte \
#     --gw_thresh   5e-8 \
#     --sugg_thresh 1e-5
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(topr)
  library(ggplot2)
  library(ggplot2)
})

option_list <- list(
  make_option("--results_dir", type = "character", default = NULL),
  make_option("--output_dir",  type = "character", default = NULL),
  make_option("--cell_type",   type = "character", default = NULL),
  make_option("--gw_thresh",   type = "double",    default = 5e-8),
  make_option("--sugg_thresh", type = "double",    default = 1e-5)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$results_dir)) stop("--results_dir is required")
if (is.null(opt$output_dir))  stop("--output_dir is required")
if (is.null(opt$cell_type))   stop("--cell_type is required")

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

GW   <- opt$gw_thresh
SUGG <- opt$sugg_thresh
CT   <- opt$cell_type

# ---- Find the .tbl file for this cell type ----------------------------------
tbl_files <- list.files(opt$results_dir, pattern = "\\.tbl$", full.names = TRUE)
# Exclude the *1.tbl intermediate files (keep symlinks named without trailing 1)
tbl_files  <- tbl_files[!grepl("1\\.tbl$", tbl_files)]
names(tbl_files) <- sub("_meta_analysis_.*\\.tbl$", "", basename(tbl_files))

if (!(CT %in% names(tbl_files)))
  stop(sprintf("Cell type '%s' not found. Available: %s",
               CT, paste(sort(names(tbl_files)), collapse = ", ")))

path <- tbl_files[CT]
cat(sprintf("[%s] Reading %s ...\n", CT, basename(path)))

# ---- Read & parse -----------------------------------------------------------
dt <- fread(path, header = TRUE, showProgress = FALSE) %>% as.data.frame()

df <- dt %>%
  tidyr::separate(MarkerName, into = c("CHROM", "POS", "REF", "ALT"),
                  sep = ":", remove = FALSE, extra = "drop") %>%
  rename(ID = MarkerName, P = `P-value`) %>%
  mutate(
    CHROM = suppressWarnings(as.integer(gsub("^chr", "", CHROM))),
    POS   = as.numeric(POS),
    P     = as.numeric(P)
  ) %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1)

cat(sprintf("[%s] %s variants loaded\n", CT, format(nrow(df), big.mark = ",")))

# topr expects columns named CHROM, POS, P (already correct)

# ---- Genomic inflation ------------------------------------------------------
chi2   <- qchisq(1 - df$P, 1)
lambda <- round(median(chi2, na.rm = TRUE) / qchisq(0.5, 1), 4)
n_gw   <- sum(df$P < GW,   na.rm = TRUE)
n_sugg <- sum(df$P < SUGG, na.rm = TRUE)
cat(sprintf("[%s] lambda=%.4f  GW hits=%d  Suggestive hits=%d\n",
            CT, lambda, n_gw, n_sugg))

# ---- manhattanExtra plot ----------------------------------------------------
man_file <- file.path(opt$output_dir, paste0(CT, "_manhattan_annotated.png"))
cat(sprintf("[%s] Saving annotated Manhattan -> %s\n", CT, basename(man_file)))

p_man <- manhattanExtra(
  df,
  genome_wide_thresh  = GW,
  suggestive_thresh   = SUGG,
  title               = sprintf("Meta-Analysis Manhattan: %s  (λ=%.4f | GW=%d | Sugg=%d)",
                                CT, lambda, n_gw, n_sugg),
  flank_size          = 1e6,
  region_size         = 1e6
)

ggsave(man_file, p_man, width = 18, height = 6, dpi = 300, limitsize = FALSE)
cat(sprintf("[%s] Manhattan saved.\n", CT))

# ---- QQ plot ----------------------------------------------------------------
qq_file <- file.path(opt$output_dir, paste0(CT, "_qq.png"))
cat(sprintf("[%s] Saving QQ -> %s\n", CT, basename(qq_file)))

p_qq <- qqtopr(
  df,
  title = sprintf("QQ Plot: %s  (λ=%.4f)", CT, lambda)
)

ggsave(qq_file, p_qq, width = 6, height = 6, dpi = 300)
cat(sprintf("[%s] QQ saved.\n", CT))

cat(sprintf("[%s] Done.\n", CT))
