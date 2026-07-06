#!/usr/bin/env Rscript
# Merge incremental coloc results into existing per-cell-type tables.
# Removes rows for replaced diseases and appends new results.
#
# Usage (single cell type):
#   Rscript merge_coloc_results.R \
#     --existing_dir results/coloc/coloc_results_full \
#     --new_file     results/coloc/coloc_results_updated/VIP_coloc_results.tsv \
#     --replace_diseases BD_Mullins2021,MDD_Wray2018
#
# Usage (all files in new_dir):
#   Rscript merge_coloc_results.R --existing_dir ... --new_dir ...

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--existing_dir", type = "character",
              help = "Directory with current coloc result TSVs"),
  make_option("--new_dir", type = "character", default = NULL,
              help = "Directory with new partial coloc result TSVs"),
  make_option("--new_file", type = "character", default = NULL,
              help = "Single new coloc result TSV to merge (preferred for array jobs)"),
  make_option("--replace_diseases", type = "character", default = NULL,
              help = "Comma-separated disease labels to remove from existing results"),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory [default: existing_dir]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$existing_dir))
if (is.null(opt$output_dir)) opt$output_dir <- opt$existing_dir

replace_diseases <- character(0)
if (!is.null(opt$replace_diseases) && nzchar(opt$replace_diseases)) {
  replace_diseases <- trimws(strsplit(opt$replace_diseases, ",", fixed = TRUE)[[1]])
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(opt$new_file) && nzchar(opt$new_file)) {
  new_files <- opt$new_file
} else if (!is.null(opt$new_dir) && nzchar(opt$new_dir)) {
  new_files <- list.files(opt$new_dir, pattern = "_coloc_results\\.tsv$", full.names = TRUE)
} else {
  stop("Provide --new_file or --new_dir")
}
if (length(new_files) == 0) stop("No new coloc result files found")

for (new_file in new_files) {
  cell_type <- sub("_coloc_results\\.tsv$", "", basename(new_file))
  existing_file <- file.path(opt$existing_dir, basename(new_file))
  out_file <- file.path(opt$output_dir, basename(new_file))

  new_df <- fread(new_file, sep = "\t", data.table = FALSE)
  if (file.exists(existing_file)) {
    old_df <- fread(existing_file, sep = "\t", data.table = FALSE)
    if (length(replace_diseases) > 0) {
      old_df <- old_df %>% filter(!disease %in% replace_diseases)
    }
    merged <- bind_rows(old_df, new_df) %>% arrange(desc(PP.H4))
  } else {
    merged <- new_df %>% arrange(desc(PP.H4))
  }

  fwrite(merged, out_file, sep = "\t", quote = FALSE)
  cat(sprintf("[%s] merged %d rows -> %s\n", cell_type, nrow(merged), out_file))
}

cat("Merge complete.\n")
