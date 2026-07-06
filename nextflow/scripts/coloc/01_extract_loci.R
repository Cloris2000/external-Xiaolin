#!/usr/bin/env Rscript
# Extract independent suggestive loci (p < 1e-5) from 15-cohort meta-analysis
# for each cell type, using position-based greedy clumping (500 kb windows).
#
# Outputs (per cell type):
#   results/coloc/loci/{CellType}_loci.tsv              -- lead SNP table
#   results/coloc/loci/{CellType}_locus_data.tsv.gz     -- SNPs for coloc input
#       (suggestive-only by default; all SNPs in window with --full_window_snps)
#
# Usage:
#   Rscript 01_extract_loci.R --cell_type VIP \
#       --tbl_file results/meta_analysis_15cohorts/VIP_meta_analysis_...tbl \
#       --output_dir results/coloc/loci \
#       --p_threshold 1e-5 \
#       --window_kb 500

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--cell_type",    type = "character", help = "Cell type name"),
  make_option("--tbl_file",     type = "character", help = "METAL .tbl file"),
  make_option("--output_dir",   type = "character", default = "results/coloc/loci"),
  make_option("--p_threshold",  type = "double",    default = 1e-5,
              help = "Suggestive p-value threshold [default: 1e-5]"),
  make_option("--window_kb",    type = "integer",   default = 500,
              help = "Half-window size in kb around lead SNP [default: 500]"),
  make_option("--full_window_snps", type = "logical", default = FALSE,
              help = "Include all SNPs in each locus window (recommended for coloc) [default: FALSE]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$cell_type), !is.null(opt$tbl_file))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

window_bp <- opt$window_kb * 1000L

standardize_tbl <- function(dt) {
  dt <- as.data.table(dt)
  setnames(dt, c("MarkerName", "Effect", "StdErr", "P-value"),
           c("snp", "beta", "se", "p"), skip_absent = TRUE)
  dt[, p := as.numeric(p)]
  dt <- dt[!is.na(p) & p > 0 & p <= 1]
  dt[, c("chr", "pos", "ref", "alt") := tstrsplit(snp, ":", fixed = TRUE)]
  dt[, chr := sub("^chr", "", chr, ignore.case = TRUE)]
  dt[, pos := as.integer(pos)]
  dt[!is.na(pos)]
}

read_suggestive_tbl <- function(tbl_file, p_threshold) {
  cmd <- sprintf(
    "awk -F'\\t' 'NR==1 || ($10 < %g && $10 > 0)' %s",
    p_threshold, shQuote(tbl_file)
  )
  cat(sprintf("  Reading suggestive SNPs via awk (p < %g)...\n", p_threshold))
  standardize_tbl(fread(cmd = cmd, sep = "\t", showProgress = FALSE))
}

read_full_tbl <- function(tbl_file) {
  cat(sprintf("  Reading full summary statistics...\n"))
  standardize_tbl(fread(tbl_file, sep = "\t", showProgress = FALSE))
}

read_chr_window_tbl <- function(tbl_file, chr, win_start, win_end) {
  cmd <- sprintf(
    paste0(
      "awk -F'\\t' 'NR==1 {print; next} ",
      "{ split($1,a,\":\"); c=a[1]; gsub(/^chr/,\"\",c); pos=a[2]+0; ",
      "if(c==\"%s\" && pos>=%d && pos<=%d) print }' %s"
    ),
    chr, win_start, win_end, shQuote(tbl_file)
  )
  standardize_tbl(fread(cmd = cmd, sep = "\t", showProgress = FALSE))
}

build_full_window_locus_data <- function(tbl_file, leads) {
  leads_dt <- as.data.table(leads)[, .(locus_id, chr, window_start, window_end)]
  chunks <- lapply(split(leads_dt, leads_dt$chr), function(leads_chr) {
    chr_val <- leads_chr$chr[1]
    win_start <- min(leads_chr$window_start)
    win_end <- max(leads_chr$window_end)
    cat(sprintf("  chr%s: %d loci, positions %s-%s\n",
                chr_val, nrow(leads_chr),
                format(win_start, big.mark = ","),
                format(win_end, big.mark = ",")))
    raw_ov <- read_chr_window_tbl(tbl_file, chr_val, win_start, win_end)
    if (nrow(raw_ov) == 0) return(NULL)
    raw_ov[, `:=`(start = pos, end = pos)]
    leads_ov <- leads_chr[, .(chr, start = window_start, end = window_end, locus_id)]
    setkey(leads_ov, chr, start, end)
    out <- foverlaps(raw_ov, leads_ov, type = "within", mult = "all")
    out[, c("start", "end") := NULL]
    out
  })
  bind_rows(chunks)
}

cat(sprintf("[%s] Reading %s\n", opt$cell_type, opt$tbl_file))

if (isTRUE(opt$full_window_snps)) {
  sugg <- read_suggestive_tbl(opt$tbl_file, opt$p_threshold)
} else {
  sugg <- read_full_tbl(opt$tbl_file)
  sugg <- sugg[p < opt$p_threshold]
}

cat(sprintf("[%s] Suggestive (p < %g): %d\n", opt$cell_type, opt$p_threshold, nrow(sugg)))

if (nrow(sugg) == 0) {
  cat(sprintf("[%s] No suggestive variants. Writing empty outputs.\n", opt$cell_type))
  leads_out <- file.path(opt$output_dir, paste0(opt$cell_type, "_loci.tsv"))
  write.table(data.frame(), leads_out, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(status = 0)
}

sugg_sorted <- as.data.frame(sugg)[order(sugg$p), ]
sugg_sorted$locus_id <- NA_character_
locus_counter <- 0

for (i in seq_len(nrow(sugg_sorted))) {
  if (!is.na(sugg_sorted$locus_id[i])) next
  locus_counter <- locus_counter + 1
  lead_chr <- sugg_sorted$chr[i]
  lead_pos <- sugg_sorted$pos[i]
  locus_id <- sprintf("%s_chr%s_%d", opt$cell_type, lead_chr, lead_pos)

  same_chr <- sugg_sorted$chr == lead_chr
  in_window <- abs(sugg_sorted$pos - lead_pos) <= window_bp
  assign_mask <- same_chr & in_window & is.na(sugg_sorted$locus_id)
  sugg_sorted$locus_id[assign_mask] <- locus_id
}

cat(sprintf("[%s] Independent loci: %d\n", opt$cell_type, locus_counter))

leads <- sugg_sorted %>%
  as.data.frame() %>%
  group_by(locus_id) %>%
  slice_min(p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    cell_type    = opt$cell_type,
    window_start = pmax(1L, pos - window_bp),
    window_end   = pos + window_bp,
    n_sugg_snps  = as.integer(table(sugg_sorted$locus_id)[locus_id])
  ) %>%
  select(cell_type, locus_id, chr, lead_pos = pos, lead_snp = snp,
         lead_p = p, window_start, window_end, n_sugg_snps) %>%
  arrange(chr, lead_pos)

leads_out <- file.path(opt$output_dir, paste0(opt$cell_type, "_loci.tsv"))
fwrite(leads, leads_out, sep = "\t", quote = FALSE)
cat(sprintf("[%s] Lead loci written to %s\n", opt$cell_type, leads_out))

locus_data_out <- file.path(opt$output_dir, paste0(opt$cell_type, "_locus_data.tsv.gz"))
if (isTRUE(opt$full_window_snps)) {
  cat(sprintf("[%s] Building full-window locus data (all SNPs per locus)...\n", opt$cell_type))
  locus_data <- build_full_window_locus_data(opt$tbl_file, leads)
  cat(sprintf("[%s] Full-window SNPs: %d (across %d loci)\n",
              opt$cell_type, nrow(locus_data), nrow(leads)))
} else {
  locus_data <- sugg_sorted
}
fwrite(locus_data, locus_data_out, sep = "\t", quote = FALSE)
cat(sprintf("[%s] Locus data written to %s\n", opt$cell_type, locus_data_out))

cat(sprintf("[%s] Done. %d loci across chromosomes: %s\n",
            opt$cell_type, locus_counter,
            paste(sort(unique(leads$chr)), collapse = ", ")))
