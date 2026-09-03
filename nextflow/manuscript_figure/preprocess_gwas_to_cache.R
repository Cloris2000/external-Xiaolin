# =============================================================================
# preprocess_gwas_to_cache.R
#
# One-time conversion of METAL .tbl files to a fast binary cache.
# Run this ONCE; subsequent runs of figure3_gwas_discovery.R will be
# ~10–40× faster.
#
# Priority order for cache format:
#   1. Parquet (arrow package) — ~50 MB/file, reads in <1 s
#   2. Compressed TSV (data.table) — ~200 MB/file, reads in ~10 s
#
# Output directory: CACHE_DIR (one file per cell type)
# =============================================================================

# =============================================================================
# SECTION 1 — Paths
# =============================================================================

GWAS_META_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/meta_analysis_15cohorts"
CACHE_DIR     <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure/gwas_cache"

# Columns to keep in cache (chr/pos parsed; only what figure3 needs)
KEEP_COLS_BASE  <- c("chr", "pos", "beta", "se", "p", "eaf")  # non-focus CTs
KEEP_COLS_FOCUS <- c("chr", "pos", "beta", "se", "p", "eaf",
                     "ea", "nea")                               # VIP, L5.6.IT.Car3

FOCUS_CTS <- c("VIP", "L5.6.IT.Car3", "Microglia")

# Parallel workers for conversion (one worker per file)
N_CORES <- min(6L, parallel::detectCores(logical = FALSE))

# =============================================================================
# SECTION 2 — Packages
# =============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

# Try to install / load arrow
HAS_ARROW <- requireNamespace("arrow", quietly = TRUE)
if (!HAS_ARROW) {
  cat("arrow not found — attempting install (this takes 1–3 minutes)...\n")
  tryCatch({
    install.packages("arrow", repos = "https://cloud.r-project.org",
                     quiet = TRUE)
    HAS_ARROW <- requireNamespace("arrow", quietly = TRUE)
  }, error = function(e) message("arrow install failed: ", conditionMessage(e)))
}

if (HAS_ARROW) {
  suppressPackageStartupMessages(library(arrow))
  CACHE_FORMAT <- "parquet"
  CACHE_EXT    <- ".parquet"
  cat("Cache format: Parquet (arrow", as.character(packageVersion("arrow")), ")\n")
} else {
  CACHE_FORMAT <- "tsv_gz"
  CACHE_EXT    <- ".tsv.gz"
  cat("Cache format: compressed TSV (data.table fallback)\n")
}

# =============================================================================
# SECTION 3 — Setup
# =============================================================================

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

per_files <- list.files(GWAS_META_DIR, pattern = "\\.tbl$", full.names = TRUE)
per_files <- per_files[!grepl("1\\.tbl$", per_files)]
per_files <- per_files[!grepl("\\.info$", per_files)]

cat("Found", length(per_files), "METAL files to process.\n")
cat("Output directory:", CACHE_DIR, "\n")
cat("Using", N_CORES, "parallel workers.\n\n")

# =============================================================================
# SECTION 4 — Worker function
# =============================================================================

convert_one_file <- function(f, cache_dir, cache_ext, cache_format,
                              focus_cts, keep_base, keep_focus) {
  ct_raw <- sub("_meta_.*", "", basename(f))

  # Determine output path
  out_path <- file.path(cache_dir, paste0(ct_raw, cache_ext))
  if (file.exists(out_path)) {
    return(list(ct = ct_raw, status = "skipped (already exists)",
                path = out_path, elapsed = 0))
  }

  t0 <- proc.time()

  # Read all needed columns
  is_focus <- any(sapply(focus_cts, function(fc)
    grepl(gsub("\\.", "\\\\.", fc), ct_raw, ignore.case = TRUE)))

  sel_cols <- c("MarkerName", "Effect", "StdErr", "P-value", "Freq1")
  if (is_focus) sel_cols <- c(sel_cols, "Allele1", "Allele2")

  dt <- tryCatch(
    data.table::fread(f, sep = "\t", showProgress = FALSE, select = sel_cols),
    error = function(e) {
      return(list(ct = ct_raw, status = paste("ERROR:", conditionMessage(e)),
                  path = NA, elapsed = NA))
    }
  )
  if (is.list(dt) && !is.data.table(dt)) return(dt)  # propagate error

  # Rename
  setnames(dt,
           old = c("MarkerName", "Effect", "StdErr", "P-value", "Freq1"),
           new = c("snp",        "beta",   "se",     "p",       "eaf"),
           skip_absent = TRUE)
  if ("Allele1" %in% names(dt)) setnames(dt, "Allele1", "ea",  skip_absent = TRUE)
  if ("Allele2" %in% names(dt)) setnames(dt, "Allele2", "nea", skip_absent = TRUE)

  # Parse MarkerName → chr (integer), pos (integer)
  parts    <- data.table::tstrsplit(dt$snp, ":", fixed = TRUE, keep = 1:2)
  dt[, chr := as.integer(sub("^[Cc][Hh][Rr]", "", parts[[1L]]))]
  dt[, pos := as.integer(parts[[2L]])]

  # Filter
  dt <- dt[chr %in% 1:22 & !is.na(pos) & !is.na(p) & !is.na(beta) & !is.na(se)
           & p > 0 & p <= 1 & !grepl(":\\*:", snp, fixed = FALSE)]

  # Drop snp column (already split; saves ~40% file size)
  dt[, snp := NULL]

  # Keep only needed columns
  keep <- if (is_focus) keep_focus else keep_base
  keep <- intersect(keep, names(dt))
  dt   <- dt[, ..keep]

  # Write cache
  if (cache_format == "parquet") {
    arrow::write_parquet(dt, out_path, compression = "snappy")
  } else {
    data.table::fwrite(dt, out_path, sep = "\t", compress = "gzip")
  }

  elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
  list(ct = ct_raw, status = "converted", path = out_path, elapsed = elapsed)
}

# =============================================================================
# SECTION 5 — Run conversions in parallel
# =============================================================================

cat("Starting conversion...\n")
t_start <- proc.time()

results <- parallel::mclapply(
  per_files,
  convert_one_file,
  cache_dir    = CACHE_DIR,
  cache_ext    = CACHE_EXT,
  cache_format = CACHE_FORMAT,
  focus_cts    = FOCUS_CTS,
  keep_base    = KEEP_COLS_BASE,
  keep_focus   = KEEP_COLS_FOCUS,
  mc.cores     = N_CORES,
  mc.preschedule = FALSE
)

total_elapsed <- round((proc.time() - t_start)[["elapsed"]], 0)

# =============================================================================
# SECTION 6 — Summary
# =============================================================================

cat(sprintf("\nConversion complete in %d s (%.1f min).\n\n",
            total_elapsed, total_elapsed / 60))

for (r in results) {
  if (is.list(r) && "ct" %in% names(r)) {
    cat(sprintf("  %-20s  %s  (%.1f s)\n",
                r$ct,
                r$status,
                if (is.numeric(r$elapsed)) r$elapsed else 0))
  }
}

cat("\nCache files written to:", CACHE_DIR, "\n")
cache_files <- list.files(CACHE_DIR, full.names = TRUE)
total_size  <- sum(file.info(cache_files)$size, na.rm = TRUE)
cat(sprintf("Total cache size: %.1f MB across %d files\n",
            total_size / 1e6, length(cache_files)))

cat("\n")
cat("Next step: set GWAS_CACHE_DIR in figure3_gwas_discovery.R:\n")
cat("  GWAS_CACHE_DIR <- \"", CACHE_DIR, "\"\n", sep = "")
cat("The figure script will automatically use the cache on the next run.\n")
