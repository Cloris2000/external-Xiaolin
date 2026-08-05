# =============================================================================
# Figure 3. Genome-wide association analysis of bulk-derived cell-type
#           proportion phenotypes
#
# Main figure panels:
#   A — Genome-wide signal burden across cell types   (stacked bar)
#   B — VIP Miami plot          (signed -log10 p: up = positive beta)
#   C — L5.6.IT.Car3 Miami plot
#   D — Microglia Miami plot
#
# Layout:  A / B / C / D  (each full width)
#
# Supplementary figure (GWAS diagnostics, saved separately):
#   A — Genomic inflation (lambda GC) across all cell types
#   B–D — QQ plots for the three focus cell types
#
# Miami plots replace Manhattan plots so the direction of effect is visible:
# a SNP raising the cell-type proportion (beta > 0) is plotted above the axis,
# one lowering it (beta < 0) below.  Diagnostics moved to the supplement so the
# main panels can use the full page width with legible text.
#
# Input files (define paths in SECTION 1):
#   GWAS_META_FILE        — combined meta-analysis (all cell types, optional)
#   GWAS_META_DIR         — fallback directory of per-cell-type files
#   INDEPENDENT_LOCI_FILE — pre-clumped independent loci (optional)
#   CELLTYPE_GROUP_FILE   — cell-type → broad class mapping (optional)
#   GWAS_QC_SUMMARY_FILE  — per-cell-type QC summary (optional)
#
# Output:
#   figure3_gwas_discovery.png / .pdf
#   figure3_panelA_signal_burden.tsv
#   figure3_panelB_lambda_gc.tsv
#   figure3_gws_loci_deduplicated.tsv
# =============================================================================

# =============================================================================
# SECTION 1 — File paths and parameters  (edit here)
# =============================================================================

OUT_DIR  <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure"

# --- Primary combined GWAS file ---
# Set to NULL (or a non-existent path) to fall back to per-cell-type directory.
# If you later create a combined file, set this path.
GWAS_META_FILE <- file.path(OUT_DIR, "gwas_meta_all_celltypes.tsv")

# --- Per-cell-type METAL output directory ---
# Files are named  {CellType}_meta_analysis_..._array.tbl  (METAL format).
# The script reads all *.tbl files that do NOT end in "1.tbl".
GWAS_META_DIR <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts"

# --- Pre-processed cache directory (OPTIONAL — strongly recommended) ---
# Run preprocess_gwas_to_cache.R once to create Parquet or compressed-TSV
# files with chr/pos already parsed.  Set this path to enable the fast
# read path; leave NULL to fall back to reading the raw METAL .tbl files.
# Expected files: {CellType}.parquet  OR  {CellType}.tsv.gz
GWAS_CACHE_DIR <- file.path(
  "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manuscript_figure",
  "gwas_cache"
)

# --- Optional: pre-clumped independent loci ---
# Expected columns: cell_type, chr, lead_pos, lead_snp, lead_p, lead_beta, locus_id
# Optional columns: nearest_gene, locus_label, start, end
INDEPENDENT_LOCI_FILE <- file.path(OUT_DIR, "independent_loci.tsv")

# --- Optional: cell-type grouping ---
# Expected columns: cell_type, broad_class
CELLTYPE_GROUP_FILE <- file.path(OUT_DIR, "celltype_grouping.tsv")

# --- Optional: per-cell-type QC summary ---
# Expected columns: cell_type, lambda_gc, n_snps, n_genomewide, n_suggestive, min_p
GWAS_QC_SUMMARY_FILE <- file.path(OUT_DIR, "gwas_qc_summary.tsv")

# --- Representative / focus cell types ---
# Selected from coloc PP.H4 >= 0.5 hits (MDD2025 / BD bip2024 / SCZ):
#   VIP         — MDD × TMEM106B chr7 (PP.H4 ≈ 1.0); strong VIP-specific GWAS peak
#   L5.6.IT.Car3 — BD/SCZ × CACNA1C chr12 (PP.H4 0.88/0.63); cell-type-specific, not shared
#   Microglia   — MDD × TMEM106B + BD × PRKN chr6 (PP.H4 1.0/0.62)
# (VLMC/LAMP5/etc. share the same TMEM106B signal — redundant for Manhattan panels.)
representative_cell_types <- c("VIP", "L5.6.IT.Car3", "Microglia")

# --- Canonical biological cell-type order ---
CELLTYPE_ORDER_CANONICAL <- c(
  "Oligodendrocyte", "OPC", "Astrocyte", "Microglia",
  "Endothelial", "Pericyte", "VLMC",
  "IT", "L4 IT", "L5 6 NP", "L5 ET",
  "L6 CT", "L5 6 IT Car3", "L6B",
  "Pvalb", "SST", "VIP", "Lamp5", "Pax6"
)

# --- GWAS parameters ---
PVAL_GWS       <- 5e-8    # genome-wide significance threshold
PVAL_SUGG      <- 1e-5    # suggestive threshold
CLUMP_WINDOW   <- 500e3   # ±500 kb greedy clumping window (bp)
NEG_LOG10_CAP  <- 20      # cap for plotting only
BIN_SIZE_MB    <- 2       # genomic bin size for Panel B heatmap (Mb)
BIN_SIZE_BP    <- BIN_SIZE_MB * 1e6

# --- Manhattan labelling ---
# Inf = label all independent loci (after ±500 kb clumping) per focus cell type
MANHATTAN_N_LABELS      <- Inf   # genome-wide significant (bold italic)
MANHATTAN_N_LABELS_SUGG <- Inf   # suggestive-only (plain italic)

# --- Parallel workers ---
# Each worker reads one cell-type file independently.
# Set to 1 to disable parallelism (e.g. for debugging).
# 6 workers on a 12-core machine gives each worker 2 cores of effective
# memory bandwidth, reducing contention vs. 8 workers.
N_CORES <- min(6L, parallel::detectCores(logical = FALSE))

# --- Figure dimensions ---
# Narrower than the previous 18x24 layout: dropping the side-by-side QQ plots
# frees the full width for each Miami panel, so the same text renders larger
# relative to the figure.
FIG_WIDTH  <- 14
FIG_HEIGHT <- 21   # room for Panel A tag/legend rows + Miami legend
FIG_DPI    <- 300

# --- Supplementary (diagnostics) figure dimensions ---
SUPP_WIDTH  <- 14
SUPP_HEIGHT <- 11

# =============================================================================
# SECTION 2 — Package loading
# =============================================================================

required_pkgs <- c("data.table", "tidyverse", "ggplot2", "patchwork",
                   "cowplot", "ggrepel", "scales", "RColorBrewer", "parallel", "topr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  Installing: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

HAS_JANITOR <- requireNamespace("janitor", quietly = TRUE)
if (HAS_JANITOR) {
  suppressPackageStartupMessages(library(janitor))
} else {
  message("Note: 'janitor' not found; using base-R fallback for name cleaning.")
}

HAS_VIRIDIS <- requireNamespace("viridis", quietly = TRUE)
if (HAS_VIRIDIS) suppressPackageStartupMessages(library(viridis))

cat("Packages loaded.\n")

# A Miami panel holds ~15M points. Kept as vectors, the PDF never finishes
# writing, so point layers are rasterised while axes and text stay vector.
HAS_RASTR <- requireNamespace("ggrastr", quietly = TRUE)
if (!HAS_RASTR)
  message("  ggrastr unavailable: PDF output will be skipped (PNG still written).")

raster_points <- function(layer, dpi = 300) {
  if (HAS_RASTR) ggrastr::rasterise(layer, dpi = dpi) else layer
}

# =============================================================================
# SECTION 3 — Helper functions
# =============================================================================

clean_names_fn <- function(df) {
  if (HAS_JANITOR) return(janitor::clean_names(df))
  nms <- sub("^_+|_+$", "",
             tolower(gsub("[^a-zA-Z0-9]+", "_",
                          gsub("([a-z])([A-Z])", "\\1_\\2", names(df)))))
  names(df) <- nms; df
}

make_clean_name <- function(x) {
  sub("^_+|_+$", "",
      tolower(gsub("[^a-zA-Z0-9]+", "_",
                   gsub("([a-z])([A-Z])", "\\1_\\2", x))))
}

# Normalise a cell-type string to a lowercase underscore key for matching
normalise_ct <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]+", "_", trimws(x)))
}

# Build a lookup: normalised canonical name → display label
build_ct_lookup <- function(canonical) {
  setNames(canonical, normalise_ct(canonical))
}

# Map raw cell_type strings onto canonical display labels
map_to_canonical <- function(raw, lookup) {
  k <- normalise_ct(raw)
  m <- lookup[k]
  ifelse(is.na(m), raw, m)   # fall back to raw if no match
}

# Placeholder panel for missing data
make_placeholder <- function(msg, title = "") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 3.5, hjust = 0.5, vjust = 0.5,
             color = "grey50", fontface = "italic") +
    labs(title = title) +
    theme_void() +
    theme(
      panel.border    = element_rect(color = "grey75", fill = NA, linewidth = 0.5),
      plot.title      = element_text(size = 10, face = "bold"),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

# Greedy ±window clumping within a single cell-type × chromosome group
# Returns row indices (1-based within the dt) that survive clumping.
greedy_clump <- function(dt, p_col = "p", pos_col = "pos", window = CLUMP_WINDOW) {
  # dt must already be sorted by p ascending
  keep <- logical(nrow(dt))
  used <- rep(FALSE, nrow(dt))
  pos  <- dt[[pos_col]]
  for (i in seq_len(nrow(dt))) {
    if (!used[i]) {
      keep[i] <- TRUE
      # mark all within window as used
      diffs <- abs(pos - pos[i])
      used[diffs <= window] <- TRUE
    }
  }
  which(keep)
}

# Derive independent loci by greedy clumping from full GWAS data
derive_loci <- function(gwas_dt, p_thresh = PVAL_GWS, window = CLUMP_WINDOW) {
  sig <- gwas_dt[p < p_thresh]
  if (nrow(sig) == 0L) return(data.table())
  setorder(sig, p)
  # Keep EVERY independent lead on a chromosome (not just the strongest).
  # Using snp[1]/pos[1] here previously collapsed multi-peak chromosomes
  # to a single locus and left many suggestive Miami peaks unlabeled.
  result <- sig[, {
    idx <- greedy_clump(.SD, window = window)
    .SD[idx, .(lead_snp = snp, lead_pos = pos,
               lead_p = p, lead_beta = beta)]
  }, by = .(cell_type, chr)]
  result[, locus_id := paste0(cell_type, "_", chr, "_", lead_pos)]
  result
}

# Lambda GC from p-values
compute_lambda <- function(p) {
  p  <- p[!is.na(p) & p > 0 & p < 1]
  if (length(p) < 100L) return(NA_real_)
  chisq <- qchisq(p, df = 1, lower.tail = FALSE)
  round(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE), 4)
}

# Chromosome cumulative offset for Manhattan layout
make_chr_offsets <- function(gwas_dt) {
  chr_max <- gwas_dt[, .(chr_len = max(pos, na.rm = TRUE)), by = chr]
  setorder(chr_max, chr)
  chr_max[, offset := cumsum(as.numeric(shift(chr_len, fill = 0)))]
  chr_max
}

# Common ggplot2 theme
theme_fig3 <- function(base_size = 9) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      margin = margin(b = 4)),
      axis.text        = element_text(size = base_size - 1),
      axis.title       = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1.5),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.key.size  = unit(0.35, "cm"),
      strip.text       = element_text(size = base_size - 1),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# Parse a METAL .tbl file into standardised columns.
# METAL columns: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, P-value
# MarkerName format: chrN:POS:REF:ALT  (e.g. chr11:12541586:A:G)
parse_metal_tbl <- function(path, cell_type_label) {
  dt <- data.table::fread(path, sep = "\t", data.table = TRUE,
                           showProgress = FALSE)
  # Rename METAL columns → standard names
  setnames(dt,
           old = c("MarkerName", "Allele1", "Allele2", "Freq1",
                   "Effect", "StdErr", "P-value"),
           new = c("snp",        "ea",      "nea",     "eaf",
                   "beta",       "se",      "p"),
           skip_absent = TRUE)

  # Parse MarkerName: "chrN:POS:REF:ALT"
  parts <- data.table::tstrsplit(dt$snp, ":", fixed = TRUE)
  dt[, chr := sub("^[Cc][Hh][Rr]", "", parts[[1L]])]
  dt[, pos := as.integer(parts[[2L]])]

  # Keep SNPs only (skip indels / multiallelic with "*" in alleles)
  dt <- dt[!grepl("\\*", snp, fixed = FALSE)]

  # Attach cell type
  dt[, cell_type := cell_type_label]

  # Return only the columns needed downstream
  keep <- c("cell_type", "chr", "pos", "snp", "beta", "se", "p",
            intersect(c("ea", "nea", "eaf"), names(dt)))
  dt[, ..keep]
}

cat("Helper functions defined.\n")

# =============================================================================
# SECTION 4 — Streaming read of GWAS summary statistics
#
# Each cell-type file (~14.8M rows) is processed in isolation to avoid
# loading all ~280M rows into memory simultaneously.
#
# Per file, we accumulate:
#   gwas_focus   — full SNP data for focus cell types only (Panels C & D)
#   gwas_sig     — variants with p < P_KEEP across all cell types (Panels A, E)
#   binned_list  — max -log10(p) per 2 Mb bin × cell type (Panel B)
#   chr_len_list — max position per chromosome (for layout coordinates)
#   lambda_list  — lambda GC per cell type (for QQ inset)
# =============================================================================

P_KEEP <- 1e-4   # retain variants with p < this threshold for Panel E lookups

# RDS checkpoint: save/load assembled data to skip 45-min loading on re-runs.
# Delete this file to force a fresh load from raw or cache files.
DATA_CHECKPOINT <- file.path(OUT_DIR, "figure3_data_checkpoint.rds")

ct_lookup <- build_ct_lookup(CELLTYPE_ORDER_CANONICAL)

if (file.exists(DATA_CHECKPOINT)) {
  cat("Loading data from checkpoint:", DATA_CHECKPOINT, "\n")
  chk       <- readRDS(DATA_CHECKPOINT)
  binned    <- chk$binned
  chr_stats <- chk$chr_stats
  gwas_sig  <- chk$gwas_sig
  gwas_focus <- chk$gwas_focus
  qc_summary <- chk$qc_summary
  ct_order_use <- chk$ct_order_use
  cat("Checkpoint loaded. Skipping file I/O.\n")
} else {

# ---- 4a. Optional: broad class grouping ----
ct_group <- NULL
if (file.exists(CELLTYPE_GROUP_FILE)) {
  ct_group <- data.table::fread(CELLTYPE_GROUP_FILE)
  ct_group <- clean_names_fn(ct_group)
  ct_group[, cell_type := map_to_canonical(cell_type, ct_lookup)]
  cat("Cell-type grouping loaded.\n")
}

# ---- 4b. Optional: QC summary ----
qc_summary <- NULL
if (file.exists(GWAS_QC_SUMMARY_FILE)) {
  qc_summary <- data.table::fread(GWAS_QC_SUMMARY_FILE)
  qc_summary <- clean_names_fn(qc_summary)
  qc_summary[, cell_type := map_to_canonical(cell_type, ct_lookup)]
  cat("QC summary loaded.\n")
}

# ---- 4c. Locate per-cell-type files ----
# Priority: combined file > pre-processed cache > raw METAL .tbl files

USE_CACHE  <- FALSE
CACHE_FORMAT_DETECTED <- NULL

if (file.exists(GWAS_META_FILE)) {
  per_files <- GWAS_META_FILE
  cat("Using combined GWAS file:", GWAS_META_FILE, "\n")

} else if (!is.null(GWAS_CACHE_DIR) && dir.exists(GWAS_CACHE_DIR)) {
  # Check for Parquet or compressed-TSV cache files
  pq_files  <- list.files(GWAS_CACHE_DIR, pattern = "\\.parquet$", full.names = TRUE)
  tsv_files <- list.files(GWAS_CACHE_DIR, pattern = "\\.tsv\\.gz$", full.names = TRUE)

  if (length(pq_files) >= 5L && requireNamespace("arrow", quietly = TRUE)) {
    per_files             <- pq_files
    USE_CACHE             <- TRUE
    CACHE_FORMAT_DETECTED <- "parquet"
    suppressPackageStartupMessages(library(arrow))
    cat("Using Parquet cache (", length(per_files), "files) — fast path enabled.\n")
  } else if (length(tsv_files) >= 5L) {
    per_files             <- tsv_files
    USE_CACHE             <- TRUE
    CACHE_FORMAT_DETECTED <- "tsv_gz"
    cat("Using compressed-TSV cache (", length(per_files), "files) — fast path enabled.\n")
  } else {
    cat("Cache directory exists but no cache files found; falling back to raw .tbl files.\n")
  }
}

if (!USE_CACHE && !file.exists(GWAS_META_FILE)) {
  cat("Reading raw METAL .tbl files from:\n ", GWAS_META_DIR, "\n")
  if (!dir.exists(GWAS_META_DIR))
    stop("GWAS_META_DIR not found: ", GWAS_META_DIR)
  per_files <- list.files(GWAS_META_DIR, pattern = "\\.tbl$", full.names = TRUE)
  per_files <- per_files[!grepl("1\\.tbl$", per_files)]
  per_files <- per_files[!grepl("\\.info$",  per_files)]
  if (length(per_files) == 0L)
    stop("No .tbl files in GWAS_META_DIR: ", GWAS_META_DIR)
  cat("  Found", length(per_files), "cell-type files.\n")
}

# ---- 4d. Helper: standardise a single loaded data.table ----
standardise_gwas_dt <- function(dt, ct_label, ct_lookup_map) {
  if (!"cell_type" %in% names(dt)) dt[, cell_type := ct_label]

  # Rename METAL columns if present
  setnames(dt,
           old = c("MarkerName", "Allele1", "Allele2", "Freq1",
                   "Effect",     "StdErr",  "P-value"),
           new = c("snp",        "ea",      "nea",     "eaf",
                   "beta",       "se",      "p"),
           skip_absent = TRUE)

  # Parse MarkerName to chr / pos if needed
  if ("snp" %in% names(dt) && (!"chr" %in% names(dt) || !"pos" %in% names(dt))) {
    parts <- data.table::tstrsplit(dt$snp, ":", fixed = TRUE)
    dt[, chr := sub("^[Cc][Hh][Rr]", "", parts[[1L]])]
    dt[, pos := as.integer(parts[[2L]])]
  }

  # Strip "chr" prefix; keep autosomes 1–22
  dt[, chr := as.integer(sub("^[Cc][Hh][Rr]", "", as.character(chr)))]

  # Coerce numerics
  for (col in intersect(c("pos", "beta", "se", "p"), names(dt)))
    dt[, (col) := as.numeric(get(col))]

  # Filter
  dt <- dt[chr %in% 1:22 & !is.na(pos) & !is.na(p) & !is.na(beta) & !is.na(se)]
  dt <- dt[p > 0 & p <= 1]

  # Remove indels / multiallelic markers with "*"
  if ("snp" %in% names(dt))
    dt <- dt[!grepl("*", snp, fixed = TRUE)]

  # Map to canonical cell-type label
  dt[, cell_type := map_to_canonical(cell_type, ct_lookup_map)]

  dt
}

# ---- 4e. Parallel processing — one worker per cell-type file ----
focus_norm_keys <- normalise_ct(representative_cell_types)

# Worker function: reads one file, returns compressed summaries.
# Performance design:
#   1. Read only the 4 essential columns (MarkerName, Effect, StdErr, P-value)
#      for non-focus cell types — reduces fread work by ~70%.
#      For focus cell types (VIP, L5.6.IT.Car3) read all columns.
#   2. Map cell_type label ONCE per file (not 14.4 M times).
#   3. Compute lambda from a 100k-row sample, not all rows.
process_one_file <- function(f,
                              gwas_meta_file   = GWAS_META_FILE,
                              ct_lookup_map    = ct_lookup,
                              focus_keys       = focus_norm_keys,
                              bin_size         = BIN_SIZE_BP,
                              neg_log10_cap    = NEG_LOG10_CAP,
                              p_keep           = P_KEEP,
                              compute_lam      = is.null(qc_summary),
                              use_cache        = USE_CACHE,
                              cache_fmt        = CACHE_FORMAT_DETECTED) {
  ct_raw       <- if (identical(f, gwas_meta_file)) "combined"
                  else sub("_meta_.*|\\.tbl$|\\.parquet$|\\.tsv\\.gz$", "", basename(f))
  ct_canonical <- map_to_canonical(ct_raw, ct_lookup_map)
  is_focus     <- normalise_ct(ct_canonical) %in% focus_keys

  # ---- Read file (cache or raw METAL) ----
  dt <- tryCatch({
    if (use_cache && cache_fmt == "parquet") {
      # Parquet: chr/pos already parsed, no string work needed
      as.data.table(arrow::read_parquet(f))
    } else if (use_cache && cache_fmt == "tsv_gz") {
      # Pre-split TSV: numeric columns, no MarkerName parsing
      data.table::fread(f, showProgress = FALSE)
    } else {
      # Raw METAL .tbl: read selected columns only
      core_cols  <- c("MarkerName", "Effect", "StdErr", "P-value")
      extra_cols <- c("Allele1", "Allele2", "Freq1")
      sel_cols   <- if (is_focus) c(core_cols, extra_cols) else core_cols
      data.table::fread(f, sep = "\t", showProgress = FALSE, select = sel_cols)
    }
  }, error = function(e) {
    message("WARNING: failed to read ", basename(f), ": ", conditionMessage(e))
    NULL
  })
  if (is.null(dt)) return(NULL)

  # ---- Standardise columns (raw METAL only) ----
  if (!use_cache) {
    setnames(dt,
             old = c("MarkerName","Effect","StdErr","P-value"),
             new = c("snp",       "beta",  "se",    "p"),
             skip_absent = TRUE)
    if ("Allele1" %in% names(dt)) setnames(dt, "Allele1", "ea",  skip_absent = TRUE)
    if ("Allele2" %in% names(dt)) setnames(dt, "Allele2", "nea", skip_absent = TRUE)
    if ("Freq1"   %in% names(dt)) setnames(dt, "Freq1",   "eaf", skip_absent = TRUE)

    # Parse MarkerName: "chrN:POS:REF:ALT"
    parts <- data.table::tstrsplit(dt$snp, ":", fixed = TRUE, keep = 1:2)
    dt[, chr := as.integer(sub("^[Cc][Hh][Rr]", "", parts[[1L]]))]
    dt[, pos := as.integer(parts[[2L]])]

    if (!is.numeric(dt$beta)) dt[, beta := as.numeric(beta)]
    if (!is.numeric(dt$se))   dt[, se   := as.numeric(se)]
    if (!is.numeric(dt$p))    dt[, p    := as.numeric(p)]

    # Filter
    dt <- dt[chr %in% 1:22 & !is.na(pos) & !is.na(p) & !is.na(beta) & !is.na(se)
             & p > 0 & p <= 1]
    dt <- dt[!grepl(":\\*:", snp, fixed = FALSE)]
  }

  n_snps <- nrow(dt)

  # ---- Attach canonical cell-type label (single scalar, fast) ----
  dt[, cell_type := ct_canonical]

  # ---- Lambda GC: sample 100k rows to avoid qchisq on 14.4M values ----
  lam_dt <- if (compute_lam) {
    p_samp <- if (n_snps > 100000L) sample(dt$p, 100000L) else dt$p
    lam    <- compute_lambda(p_samp)
    data.table(cell_type = ct_canonical, lambda_gc = lam)
  } else NULL

  # ---- Chromosome lengths ----
  chr_len_dt <- dt[, .(chr_max = max(pos, na.rm = TRUE)), by = chr]

  # ---- Binned heatmap: max -log10(p) per 2 Mb bin ----
  dt[, bin := floor(pos / bin_size)]
  bnd <- dt[, .(max_nlp = pmin(neg_log10_cap, max(-log10(p), na.rm = TRUE))),
            by = .(chr, bin)]
  bnd[, cell_type := ct_canonical]

  # ---- Significant subset (p < P_KEEP) for loci clumping and Panel E ----
  sig_cols <- intersect(c("cell_type","chr","pos","snp","beta","se","p"), names(dt))
  sig      <- dt[p < p_keep, ..sig_cols]

  # ---- Full data for focus cell types (Panels C & D Manhattan) ----
  focus_full <- if (is_focus) {
    keep <- intersect(c("cell_type","chr","pos","snp","beta","se","p","ea","nea","eaf"),
                      names(dt))
    dt[, ..keep]
  } else NULL

  list(
    ct        = ct_canonical,
    n_snps    = n_snps,
    binned    = bnd,
    sig       = sig,
    focus     = focus_full,
    chr_len   = chr_len_dt,
    lambda    = lam_dt
  )
}

cat("Launching", N_CORES, "parallel workers for", length(per_files), "files...\n")
t0 <- proc.time()

results_par <- parallel::mclapply(
  per_files,
  process_one_file,
  use_cache      = USE_CACHE,
  cache_fmt      = CACHE_FORMAT_DETECTED,
  mc.cores       = N_CORES,
  mc.preschedule = FALSE
)

elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("Parallel read complete in %.0f s (%.1f min).\n",
            elapsed, elapsed / 60))

# Drop any failed workers
failed <- sapply(results_par, is.null)
if (any(failed))
  message("WARNING: ", sum(failed), " file(s) failed to load.")
results_par <- results_par[!failed]

# Unpack results into named lists
binned_list   <- lapply(results_par, `[[`, "binned")
gwas_sig_list <- lapply(results_par, `[[`, "sig")
chr_len_list  <- lapply(results_par, `[[`, "chr_len")
lambda_list   <- Filter(Negate(is.null), lapply(results_par, `[[`, "lambda"))
gwas_focus_list <- Filter(Negate(is.null), lapply(results_par, `[[`, "focus"))

for (r in results_par)
  cat("  ", r$ct, ":", r$n_snps, "SNPs\n")

# ---- 4f. Assemble accumulated data ----

# Binned heatmap table
binned <- data.table::rbindlist(binned_list, use.names = TRUE)
binned[, bin_mid := (bin + 0.5) * BIN_SIZE_BP]

# Chromosome reference table (union of all cell types, take max per chr)
chr_all <- data.table::rbindlist(chr_len_list, use.names = TRUE)
chr_stats <- chr_all[, .(chr_max = max(chr_max)), by = chr]
setorder(chr_stats, chr)
chr_stats <- chr_stats[chr %in% 1:22]
chr_stats[, offset  := cumsum(as.numeric(shift(chr_max, fill = 0)))]
chr_stats[, chr_mid := offset + as.numeric(chr_max) / 2]

# Significant variants
gwas_sig <- data.table::rbindlist(gwas_sig_list, fill = TRUE, use.names = TRUE)
gwas_sig[, neg_log10_p := -log10(p)]

# Full focus data
gwas_focus <- if (length(gwas_focus_list) > 0L) {
  data.table::rbindlist(gwas_focus_list, fill = TRUE, use.names = TRUE)
} else {
  data.table()
}

if (nrow(gwas_focus) > 0L) gwas_focus[, neg_log10_p := -log10(p)]

# Lambda GC table (if not loaded from file)
if (is.null(qc_summary) && length(lambda_list) > 0L) {
  qc_summary <- data.table::rbindlist(lambda_list, use.names = TRUE)
  cat("Lambda GC computed for", nrow(qc_summary), "cell types.\n")
}

# ---- 4g. Cell-type display order ----
ct_present   <- unique(c(binned$cell_type, gwas_sig$cell_type))
ct_order_use <- CELLTYPE_ORDER_CANONICAL[CELLTYPE_ORDER_CANONICAL %in% ct_present]
ct_extra     <- setdiff(ct_present, ct_order_use)
ct_order_use <- c(ct_order_use, sort(ct_extra))

binned[,   cell_type := factor(cell_type, levels = ct_order_use)]
gwas_sig[, cell_type := factor(cell_type, levels = ct_order_use)]
if (nrow(gwas_focus) > 0L)
  gwas_focus[, cell_type := factor(cell_type, levels = ct_order_use)]

cat("Cell types in data:", paste(ct_order_use, collapse = ", "), "\n")
cat("Significant variants (p <", P_KEEP, "):", nrow(gwas_sig), "\n")
cat("Focus cell-type rows:", nrow(gwas_focus), "\n")

  # Save checkpoint for fast re-runs
  saveRDS(list(binned      = binned,
               chr_stats   = chr_stats,
               gwas_sig    = gwas_sig,
               gwas_focus  = gwas_focus,
               qc_summary  = qc_summary,
               ct_order_use = ct_order_use),
          DATA_CHECKPOINT)
  cat("Data checkpoint saved:", DATA_CHECKPOINT, "\n")

} # end if/else checkpoint block

# =============================================================================
# SECTION 5 — Read or derive independent loci
# =============================================================================

loci <- NULL

if (file.exists(INDEPENDENT_LOCI_FILE)) {
  cat("Reading independent loci file:", INDEPENDENT_LOCI_FILE, "\n")
  loci <- data.table::fread(INDEPENDENT_LOCI_FILE)
  loci <- clean_names_fn(loci)
  loci[, chr := as.integer(sub("^[Cc][Hh][Rr]", "", as.character(chr)))]
  loci <- loci[!is.na(chr) & chr %in% 1:22]
  loci[, cell_type := map_to_canonical(cell_type, ct_lookup)]
  for (col in intersect(c("lead_pos", "lead_p", "lead_beta"), names(loci)))
    loci[, (col) := as.numeric(get(col))]
  cat("  Independent loci loaded:", nrow(loci), "\n")
} else {
  # Derive by greedy ±500 kb clumping from gwas_sig (p < 5e-8)
  cat("Independent loci file not found; deriving by greedy ±",
      CLUMP_WINDOW / 1e3, "kb clumping (p < 5e-8) from significant subset.\n")
  loci <- derive_loci(gwas_sig, p_thresh = PVAL_GWS)
  if (nrow(loci) > 0L) {
    cat("  Derived GWS loci:", nrow(loci), "\n")
  } else {
    cat("  No genome-wide significant variants found at p < 5e-8.\n")
  }
}

# Derive suggestive loci (p < 1e-5) for Panel A stacked bar
loci_sugg <- derive_loci(gwas_sig, p_thresh = PVAL_SUGG)

# Apply cell-type factor and add nearest_gene placeholder where absent
for (lt in list(loci, loci_sugg)) {
  if (!is.null(lt) && nrow(lt) > 0L) {
    lt[, cell_type := factor(map_to_canonical(as.character(cell_type), ct_lookup),
                              levels = ct_order_use)]
    if (!"nearest_gene" %in% names(lt))
      lt[, nearest_gene := paste0(chr, ":", lead_pos)]
  }
}
if (!is.null(loci)     && nrow(loci)     > 0L) {
  loci[, cell_type := factor(map_to_canonical(as.character(cell_type), ct_lookup),
                               levels = ct_order_use)]
  if (!"nearest_gene" %in% names(loci))
    loci[, nearest_gene := paste0(chr, ":", lead_pos)]
}
if (nrow(loci_sugg) > 0L) {
  loci_sugg[, cell_type := factor(map_to_canonical(as.character(cell_type), ct_lookup),
                                    levels = ct_order_use)]
  if (!"nearest_gene" %in% names(loci_sugg))
    loci_sugg[, nearest_gene := paste0(chr, ":", lead_pos)]
}

cat("GWS loci:", if (!is.null(loci) && nrow(loci) > 0L) nrow(loci) else 0,
    "| Suggestive loci:", nrow(loci_sugg), "\n")

# =============================================================================
# SECTION 6 — Focus cell-type subsets for Manhattan / QQ panels
# =============================================================================

# Helper: flexible match of a focus label against gwas$cell_type
match_focus <- function(label, ct_vec) {
  k <- normalise_ct(label)
  m <- ct_vec[normalise_ct(as.character(ct_vec)) == k]
  if (length(m) == 0L) NULL else m[1]
}

focus_matched <- lapply(representative_cell_types, match_focus,
                        ct_vec = ct_order_use)
names(focus_matched) <- representative_cell_types

cat("Focus cell types matched:\n")
for (nm in names(focus_matched))
  cat("  ", nm, "->", if (is.null(focus_matched[[nm]])) "NOT FOUND" else focus_matched[[nm]], "\n")

# ---- 6b. Load full GWAS for focus cell types missing from checkpoint ----
load_one_focus_file <- function(ct_canonical, ct_lookup_map) {
  pattern <- paste0("^", gsub("\\.", "\\\\.", ct_canonical))
  f <- list.files(GWAS_META_DIR, pattern = "\\.tbl$", full.names = TRUE)
  f <- f[!grepl("1\\.tbl$", f) & !grepl("\\.info$", f)]
  f <- f[grepl(pattern, basename(f), ignore.case = TRUE)]
  if (length(f) == 0L) {
    alt <- gsub(" ", ".", ct_canonical)
    pattern2 <- paste0("^", gsub("\\.", "\\\\.", alt))
    f <- list.files(GWAS_META_DIR, pattern = "\\.tbl$", full.names = TRUE)
    f <- f[!grepl("1\\.tbl$", f) & !grepl("\\.info$", f)]
    f <- f[grepl(pattern2, basename(f), ignore.case = TRUE)]
  }
  if (length(f) == 0L) return(NULL)
  f <- f[1]

  cat("  On-demand load for", ct_canonical, ":", basename(f), "\n")
  sel_cols <- c("MarkerName", "Effect", "StdErr", "P-value", "Allele1", "Allele2", "Freq1")
  dt <- data.table::fread(f, sep = "\t", showProgress = FALSE, select = sel_cols)
  setnames(dt,
           old = c("MarkerName", "Effect", "StdErr", "P-value", "Freq1",
                   "Allele1", "Allele2"),
           new = c("snp", "beta", "se", "p", "eaf", "ea", "nea"),
           skip_absent = TRUE)
  parts <- data.table::tstrsplit(dt$snp, ":", fixed = TRUE, keep = 1:2)
  dt[, chr := as.integer(sub("^[Cc][Hh][Rr]", "", parts[[1L]]))]
  dt[, pos := as.integer(parts[[2L]])]
  dt[, beta := as.numeric(beta)]
  dt[, se   := as.numeric(se)]
  dt[, p    := as.numeric(p)]
  dt <- dt[chr %in% 1:22 & !is.na(pos) & !is.na(p) & !is.na(beta) & !is.na(se)
           & p > 0 & p <= 1 & !grepl(":\\*:", snp, fixed = FALSE)]
  dt[, cell_type := ct_canonical]
  keep <- intersect(c("cell_type", "chr", "pos", "snp", "beta", "se", "p",
                      "ea", "nea", "eaf"), names(dt))
  dt[, ..keep]
}

focus_present <- if (nrow(gwas_focus) > 0L) unique(as.character(gwas_focus$cell_type)) else character()
for (fl in representative_cell_types) {
  ct <- focus_matched[[fl]]
  if (is.null(ct)) next
  if (!(ct %in% focus_present)) {
    loaded <- load_one_focus_file(ct, ct_lookup)
    if (!is.null(loaded) && nrow(loaded) > 0L) {
      loaded[, neg_log10_p := -log10(p)]
      loaded[, cell_type := factor(ct, levels = ct_order_use)]
      gwas_focus <- rbind(gwas_focus, loaded, fill = TRUE)
      focus_present <- c(focus_present, ct)
      cat("  Added", nrow(loaded), "SNPs for focus cell type:", ct, "\n")
    }
  }
}

# Refresh checkpoint if new focus cell types were loaded on demand
if (file.exists(DATA_CHECKPOINT) && exists("focus_present")) {
  chk <- readRDS(DATA_CHECKPOINT)
  old_n <- nrow(chk$gwas_focus)
  if (nrow(gwas_focus) > old_n) {
    chk$gwas_focus <- gwas_focus
    saveRDS(chk, DATA_CHECKPOINT)
    cat("Checkpoint updated with additional focus cell-type data.\n")
  }
}

get_focus_gwas <- function(label) {
  ct <- focus_matched[[label]]
  if (is.null(ct) || nrow(gwas_focus) == 0L) return(NULL)
  gwas_focus[as.character(cell_type) == ct]
}

# =============================================================================
# SECTION 7 — Colours and style constants
# =============================================================================

# Alternating chromosome colours (Manhattan plots)
CHR_COLS <- c("#3A6EA5", "#88A8C3")  # dark-blue / steel-blue

# Accent colour (supplementary lambda panel / legacy label helpers)
FOCUS_COL <- "#D55E00"   # vermilion (colour-blind safe)

# Shared significance colours — Panel A bars AND Miami plot points use the
# same palette so genome-wide / suggestive tiers are visually consistent.
COL_GWS  <- "#D55E00"    # genome-wide (p < 5e-8)
COL_SUGG <- "#3A6EA5"    # suggestive  (p < 1e-5)
COL_NS   <- "grey65"     # not significant (Miami legend only)

# =============================================================================
# SECTION 8 — Panel A: Genome-wide signal burden across cell types
# =============================================================================
cat("\n--- Building Panel A ---\n")

build_panelA_data <- function(loci_gws, loci_sg, ct_ord) {
  make_summary <- function(dt, label) {
    if (is.null(dt) || nrow(dt) == 0L)
      return(data.table(cell_type = factor(ct_ord, levels = ct_ord), n = 0L, tier = label))
    cnt <- dt[, .N, by = cell_type]
    full <- data.table(cell_type = factor(ct_ord, levels = ct_ord))
    cnt  <- cnt[full, on = "cell_type"]
    cnt[is.na(N), N := 0L]
    cnt[, tier := label]
    setnames(cnt, "N", "n")
    cnt
  }
  gws_dt  <- make_summary(loci_gws, "Genome-wide (p < 5e-8)")
  sugg_dt <- make_summary(loci_sg,  "Suggestive (p < 1e-5)")
  # Suggestive count = loci at p<1e-5 minus those at p<5e-8 (incremental)
  gws_map  <- setNames(gws_dt$n,  as.character(gws_dt$cell_type))
  sugg_dt[, n := pmax(0L, n - gws_map[as.character(cell_type)])]
  rbind(gws_dt, sugg_dt)
}

panelA_dat <- build_panelA_data(
  if (!is.null(loci) && nrow(loci) > 0L) loci else NULL,
  if (nrow(loci_sugg) > 0L) loci_sugg else NULL,
  ct_order_use
)

gws_totals <- panelA_dat[tier == "Genome-wide (p < 5e-8)",
                          setNames(n, as.character(cell_type))]
panelA_order <- ct_order_use[order(-gws_totals[ct_order_use])]
panelA_dat[, cell_type := factor(as.character(cell_type), levels = rev(panelA_order))]

pA <- ggplot(panelA_dat,
             aes(x = n, y = cell_type,
                 fill = factor(tier, levels = c("Suggestive (p < 1e-5)",
                                                "Genome-wide (p < 5e-8)")))) +
  geom_col(width = 0.65, position = position_stack()) +
  scale_fill_manual(
    values = c("Genome-wide (p < 5e-8)" = COL_GWS,
               "Suggestive (p < 1e-5)"  = COL_SUGG),
    guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = scales::pretty_breaks(n = 5)) +
  labs(title = NULL,
       x = "Number of independent loci", y = NULL) +
  theme_fig3(base_size = 13) +
  theme(
    plot.margin  = margin(8, 12, 4, 4),
    axis.text.y  = element_text(size = 12, colour = "black", face = "plain",
                               margin = margin(r = 6)),
    legend.position = "none"
  )

# Annotate total counts at bar ends (before wrapping into cowplot)
totals_by_ct <- panelA_dat[, .(total = sum(n)), by = cell_type]
pA_bars <- pA +
  geom_text(data = totals_by_ct[total > 0],
            aes(x = total, y = cell_type, label = total),
            inherit.aes = FALSE,
            hjust = -0.2, size = 4.2, color = "grey30")

# Panel A tag + legend drawn outside the bar plot so they cannot collide
# with y-axis labels or each other (ggplot legend spacing is unreliable).
panelA_legend <- ggplot(
  data.frame(
    key = factor(c("sugg", "gws"), levels = c("sugg", "gws")),
    lab = c("Suggestive (p < 1e-5)", "Genome-wide (p < 5e-8)"),
    x = c(1.0, 5.0),
    y = 1
  ),
  aes(x = x, y = y)
) +
  geom_point(aes(colour = key), shape = 15, size = 6.0) +
  geom_text(aes(label = lab), hjust = 0, nudge_x = 0.22,
            size = 4.0, colour = "grey20") +
  scale_colour_manual(values = c(sugg = COL_SUGG, gws = COL_GWS), guide = "none") +
  coord_cartesian(xlim = c(0.7, 8.0), ylim = c(0.6, 1.4), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(2, 8, 4, 8))

# Put "A" on its own row above the bars — never overlaps y-axis labels
panelA_tag <- cowplot::ggdraw() +
  cowplot::draw_label("A", fontface = "bold", size = 16,
                      x = 0.02, y = 0.5, hjust = 0, vjust = 0.5)

pA <- cowplot::plot_grid(
  panelA_tag,
  pA_bars,
  panelA_legend,
  ncol = 1,
  # Tag row must be tall enough that the bold "A" does not overflow
  # downward onto "Oligodendrocyte"
  rel_heights = c(0.14, 1, 0.16)
)

data.table::fwrite(panelA_dat, file.path(OUT_DIR, "figure3_panelA_signal_burden.tsv"),
                   sep = "\t")
cat("Panel A data saved.\n")

# =============================================================================
# SECTION 9 — Panel B: Genomic inflation (lambda GC) across cell types
# =============================================================================
cat("--- Building Panel B ---\n")

if (is.null(qc_summary) || nrow(qc_summary) == 0L) {
  qc_summary <- data.table(cell_type = ct_order_use, lambda_gc = NA_real_)
}

lambda_dat <- qc_summary[, .(cell_type, lambda_gc)]
full_ct    <- data.table(cell_type = ct_order_use)
lambda_dat <- merge(full_ct, lambda_dat, by = "cell_type", all.x = TRUE)
lambda_dat[, cell_type := factor(as.character(cell_type), levels = ct_order_use)]

focus_norm2 <- normalise_ct(representative_cell_types)
lambda_dat[, is_focus := normalise_ct(as.character(cell_type)) %in% focus_norm2]

lam_rng <- range(lambda_dat$lambda_gc, na.rm = TRUE)
lam_pad <- max(0.004, diff(lam_rng) * 0.15)
# Include 1.0 on the axis (reference line); all lambdas are < 1 here.
xlim_b <- c(lam_rng[1] - lam_pad, max(lam_rng[2] + lam_pad, 1.005))

# Build colour/face vectors aligned to lambda-sorted order (bottom→top on axis)
lambda_sorted_cts <- as.character(lambda_dat[order(lambda_gc), cell_type])
pB_ylab_col  <- ifelse(normalise_ct(lambda_sorted_cts) %in% focus_norm2, FOCUS_COL, "black")
pB_ylab_face <- ifelse(normalise_ct(lambda_sorted_cts) %in% focus_norm2, "bold",    "plain")

pB <- ggplot(lambda_dat,
             aes(x = lambda_gc,
                 y = reorder(cell_type, lambda_gc),
                 colour = is_focus)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point(size = 3.5) +
  geom_segment(aes(x = 1, xend = lambda_gc,
                   y = cell_type, yend = cell_type),
               colour = "grey75", linewidth = 0.6) +
  scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = FOCUS_COL),
                      guide = "none") +
  scale_x_continuous(limits = xlim_b, expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A. Genomic inflation across cell types",
       subtitle = expression(lambda[GC] ~ "(median" ~ chi^2 / "expected)"),
       x = expression(lambda[GC]), y = NULL) +
  theme_fig3(base_size = 13) +
  theme(axis.text.y = element_text(
    size   = 12,
    colour = pB_ylab_col,
    face   = pB_ylab_face
  ))

data.table::fwrite(lambda_dat, file.path(OUT_DIR, "figure3_panelB_lambda_gc.tsv"),
                   sep = "\t")
cat("Panel B data saved.\n")

# =============================================================================
# SECTION 10 — Manhattan + QQ helpers (topr, side-by-side layout)
# =============================================================================

prepare_topr_df <- function(dt_ct) {
  data.frame(CHROM = dt_ct$chr, POS = dt_ct$pos, P = dt_ct$p)
}

# Remove topr auto gene labels (GeomTextRepel layers)
strip_topr_labels <- function(p) {
  p$layers <- Filter(function(l) class(l$geom)[1] != "GeomTextRepel", p$layers)
  p
}

# Label top independent loci only (from clumped loci table or top SNPs)
add_top_locus_labels <- function(p, dt_ct, loci_dt, cell_label,
                                 n_label = MANHATTAN_N_LABELS) {
  if (n_label <= 0L) return(p)

  pts <- data.table::rbindlist(lapply(p$layers, function(l) {
    if (class(l$geom)[1] == "GeomPoint") as.data.table(l$data) else NULL
  }), fill = TRUE)
  if (nrow(pts) == 0L) return(p)

  label_src <- NULL
  if (!is.null(loci_dt) && nrow(loci_dt) > 0L) {
    label_src <- loci_dt[as.character(cell_type) == cell_label]
    if (nrow(label_src) > 0L) {
      setorder(label_src, lead_p)
      label_src <- head(label_src, n_label)
      label_src[, lbl := paste0(chr, ":", lead_pos)]
      label_src[, match_pos := lead_pos]
    }
  }
  if (is.null(label_src) || nrow(label_src) == 0L) {
    label_src <- dt_ct[order(p)][1:min(n_label, .N)]
    label_src[, lbl := paste0(chr, ":", pos)]
    label_src[, match_pos := pos]
  }

  label_pts <- label_src[, {
    sub <- pts[as.integer(CHROM) == chr]
    if (nrow(sub) == 0L) NULL else {
      idx <- which.min(abs(sub$POS - match_pos))
      data.table(x = sub$POS_orig[idx], y = sub$log10p[idx], lbl = lbl[1])
    }
  }, by = seq_len(nrow(label_src))]

  if (nrow(label_pts) == 0L) return(p)

  p + ggrepel::geom_label_repel(
    data        = label_pts,
    aes(x = x, y = y, label = lbl),
    inherit.aes = FALSE,
    size        = 2.6,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size  = 0.3,
    segment.colour = "grey40",
    fill          = alpha("white", 0.85),
    label.size    = 0.2,
    colour        = FOCUS_COL,
    max.overlaps  = 20
  )
}

make_qq_plot <- function(dt_ct, qc_dt = NULL, cell_label = NULL,
                         panel_letter = NULL) {
  if (is.null(dt_ct) || nrow(dt_ct) == 0L) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "No QQ data", size = 4))
  }

  lambda_val <- NA_real_
  if (!is.null(qc_dt) && !is.null(cell_label)) {
    r <- qc_dt[as.character(cell_type) == cell_label]
    if (nrow(r) > 0L && "lambda_gc" %in% names(r))
      lambda_val <- r$lambda_gc[1]
  }
  if (is.na(lambda_val))
    lambda_val <- compute_lambda(dt_ct$p)

  n_snps   <- nrow(dt_ct)
  observed <- sort(-log10(dt_ct$p), decreasing = TRUE)
  expected <- -log10((seq_len(n_snps)) / (n_snps + 1))
  idx      <- unique(round(seq(1, n_snps, length.out = min(5000L, n_snps))))
  qq_dt    <- data.table(exp = expected[idx], obs = observed[idx])

  ttl <- if (is.null(panel_letter)) paste0("QQ plot: ", cell_label)
         else paste0(panel_letter, ". ", cell_label)

  ggplot(qq_dt, aes(x = exp, y = obs)) +
    geom_point(size = 0.8, colour = "grey40", alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0, colour = "firebrick", linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("lambda[GC]==%.3f", lambda_val),
             parse = TRUE, hjust = -0.08, vjust = 1.4, size = 4.2) +
    labs(title = ttl,
         x = expression("Expected " * -log[10](italic(p))),
         y = expression("Observed " * -log[10](italic(p)))) +
    theme_fig3(base_size = 12) +
    theme(plot.title  = element_text(size = 13, face = "bold"),
          plot.margin = margin(6, 8, 6, 8))
}

# Nearest-gene symbols for a set of lead SNPs. Returns a character vector the
# same length as chr/pos, falling back to "chr:pos" if topr lookup fails.
nearest_gene_labels <- function(chr, pos, p) {
  fallback <- paste0(chr, ":", pos)
  out <- tryCatch({
    ann <- topr::annotate_with_nearest_gene(
      data.frame(CHROM = as.character(chr), POS = as.integer(pos), P = as.numeric(p))
    )
    g <- if ("Gene_Symbol" %in% names(ann)) as.character(ann$Gene_Symbol) else NULL
    if (is.null(g) || length(g) != length(fallback)) fallback else g
  }, error = function(e) {
    cat("    (nearest-gene lookup failed:", conditionMessage(e), "- using chr:pos)\n")
    fallback
  })
  out[is.na(out) | out == ""] <- fallback[is.na(out) | out == ""]
  out
}

# Deterministic non-overlapping label placement in data coordinates.
# Stronger signals are placed first; later labels shift in y (then x) until
# their text boxes no longer collide. Fast at draw time (no ggrepel).
# Boxes are deliberately conservative (a bit larger than the glyphs) so
# dense same-chromosome clusters do not stack on top of each other.
layout_gene_labels <- function(lab, x_range, y_ceiling,
                               char_h = 0.72, y_pad = 0.18) {
  if (is.null(lab) || nrow(lab) == 0L) return(lab)
  lab <- data.table::copy(as.data.table(lab))
  # GWS first, then stronger (smaller) p
  if (!"tier" %in% names(lab)) lab[, tier := "sugg"]
  lab[, tier_ord := ifelse(tier == "gws", 0L, 1L)]
  setorder(lab, tier_ord, p)

  x_span <- max(diff(range(x_range, na.rm = TRUE)), 1)
  # ~1.6% of genome width per character (half-width of text box)
  char_w <- x_span * 0.016
  x_pad  <- x_span * 0.004
  lab[, half_w := pmax(nchar(as.character(gene)), 2L) * char_w * 0.5]
  lab[, half_h := char_h]
  lab[, side := ifelse(signed_y >= 0, 1L, -1L)]
  lab[, x_lab := as.numeric(x)]
  lab[, y_lab := signed_y + side * 1.6]

  placed_x <- numeric(0)
  placed_y <- numeric(0)
  placed_w <- numeric(0)
  placed_h <- numeric(0)

  overlaps <- function(xl, yl, wi, hi) {
    if (!length(placed_x)) return(FALSE)
    any(abs(xl - placed_x) < (wi + placed_w + x_pad) &
        abs(yl - placed_y) < (hi + placed_h + y_pad))
  }

  # Fine vertical ladder + wider horizontal fan for crowded peaks
  dy_steps <- seq(0, y_ceiling, by = char_h * 2 + y_pad)
  dx_steps <- c(0, 1, -1, 2, -2, 3.5, -3.5, 5, -5, 7, -7) * x_span * 0.012

  for (i in seq_len(nrow(lab))) {
    xi <- lab$x[i]
    yi0 <- lab$signed_y[i] + lab$side[i] * 1.6
    wi <- lab$half_w[i]
    hi <- lab$half_h[i]
    side <- lab$side[i]
    found <- FALSE
    for (dx in dx_steps) {
      for (dy in dy_steps) {
        xl <- xi + dx
        yl <- yi0 + side * dy
        # stay on the same side of the Miami midline
        if (side > 0 && yl < 0.5) next
        if (side < 0 && yl > -0.5) next
        if (abs(yl) > y_ceiling) next
        if (overlaps(xl, yl, wi, hi)) next
        lab$x_lab[i] <- xl
        lab$y_lab[i] <- yl
        placed_x <- c(placed_x, xl)
        placed_y <- c(placed_y, yl)
        placed_w <- c(placed_w, wi)
        placed_h <- c(placed_h, hi)
        found <- TRUE
        break
      }
      if (found) break
    }
    if (!found) {
      # last resort: walk outward until free, then clamp
      for (k in seq_len(40L)) {
        xl <- xi + ((k %% 2L) * 2L - 1L) * (k %/% 2L) * x_span * 0.01
        yl <- yi0 + side * (char_h * 2.1 * k)
        yl <- max(min(yl, y_ceiling), -y_ceiling)
        if (!overlaps(xl, yl, wi, hi) &&
            ((side > 0 && yl >= 0.5) || (side < 0 && yl <= -0.5))) {
          lab$x_lab[i] <- xl
          lab$y_lab[i] <- yl
          placed_x <- c(placed_x, xl)
          placed_y <- c(placed_y, yl)
          placed_w <- c(placed_w, wi)
          placed_h <- c(placed_h, hi)
          found <- TRUE
          break
        }
      }
      if (!found) {
        yl <- max(min(yi0 + side * y_ceiling * 0.85, y_ceiling), -y_ceiling)
        lab$x_lab[i] <- xi
        lab$y_lab[i] <- yl
        placed_x <- c(placed_x, xi)
        placed_y <- c(placed_y, yl)
        placed_w <- c(placed_w, wi)
        placed_h <- c(placed_h, hi)
      }
    }
  }
  # Center glyphs on the laid-out point so boxes match rendered text
  lab[, vjust := 0.5]
  lab[, hjust := 0.5]
  lab
}

# Miami plot: Manhattan with the sign of beta applied to -log10(p), so SNPs
# that raise the cell-type proportion sit above the axis and those that lower
# it sit below.  All SNPs are plotted (no thinning, no y-axis truncation).
make_miami_plot <- function(dt_ct, cell_label, loci_dt = NULL,
                            loci_sugg_dt = NULL,
                            n_label = MANHATTAN_N_LABELS,
                            n_label_sugg = MANHATTAN_N_LABELS_SUGG,
                            show_legend = FALSE,
                            panel_tag = NULL) {
  if (is.null(dt_ct) || nrow(dt_ct) == 0L) {
    return(make_placeholder(paste("No GWAS data for", cell_label)))
  }

  dt <- data.table::copy(as.data.table(dt_ct))
  dt[, chr := suppressWarnings(as.integer(chr))]
  dt <- dt[chr %in% 1:22 & !is.na(pos) & !is.na(p) & p > 0 & !is.na(beta)]
  if (nrow(dt) == 0L) {
    return(make_placeholder(paste("No GWAS data for", cell_label)))
  }

  # ---- Orient every effect to the same allele ----------------------------
  # METAL reports Effect relative to Allele1, but Allele1 is not a consistent
  # allele across SNPs (it equals ALT only ~54% of the time), so the raw sign
  # of beta is arbitrary and the up/down split would be meaningless: at the
  # chr7 TMEM106B locus in VIP the raw betas split 17+/19-, whereas after
  # re-orientation they are 35+/1-. Re-orient to the ALT allele parsed from the
  # MarkerName (chr:pos:REF:ALT).
  if (all(c("snp", "ea") %in% names(dt))) {
    alt_allele <- toupper(data.table::tstrsplit(dt$snp, ":", fixed = TRUE,
                                                keep = 4L)[[1L]])
    flip <- !is.na(alt_allele) & toupper(dt$ea) != alt_allele
    dt[flip, beta := -beta]
    cat("    oriented to ALT allele: flipped", sum(flip), "of", nrow(dt),
        "effects\n")
    rm(alt_allele, flip)
  } else {
    warning("No allele columns for ", cell_label,
            "; Miami plot direction would be arbitrary.", call. = FALSE)
  }

  # Signed significance. Cap magnitude for plotting only.
  dt[, signed_y := -log10(p) * sign(beta)]
  dt[signed_y >  NEG_LOG10_CAP, signed_y :=  NEG_LOG10_CAP]
  dt[signed_y < -NEG_LOG10_CAP, signed_y := -NEG_LOG10_CAP]

  # Cumulative genomic coordinate
  off <- make_chr_offsets(dt)
  dt  <- merge(dt, off[, .(chr, offset)], by = "chr")
  dt[, x := as.numeric(pos) + offset]

  ticks <- dt[, .(center = (min(x) + max(x)) / 2), by = chr][order(chr)]

  # Colour key: alternating grey for null SNPs, blue suggestive, orange GWS
  dt[, pt_col := ifelse(
    p < PVAL_GWS,  "gws",
    ifelse(p < PVAL_SUGG, "sugg",
           ifelse(chr %% 2L == 1L, "ns_odd", "ns_even")))]
  # Draw null points first so significant SNPs sit on top
  dt[, draw_order := match(pt_col, c("ns_odd", "ns_even", "sugg", "gws"))]
  setorder(dt, draw_order)

  # Thin only non-significant SNPs for draw/save speed (keep every 40th).
  # All suggestive + GWS points are retained — scientific peaks unchanged.
  n_full <- nrow(dt)
  keep_ns <- dt$pt_col %in% c("sugg", "gws")
  ns_idx  <- which(!keep_ns)
  if (length(ns_idx) > 0L) {
    keep_ns[ns_idx[seq(1L, length(ns_idx), by = 40L)]] <- TRUE
  }
  dt <- dt[keep_ns]
  cat("  ", cell_label, "Miami SNPs:", nrow(dt),
      "(from", n_full, "; null thinned 1/40 for rendering)\n")

  thr_gws  <- -log10(PVAL_GWS)
  thr_sugg <- -log10(PVAL_SUGG)

  # Gene labels: GWS = bold italic; suggestive-only = plain italic
  pick_leads <- function(loci_src, p_lo, p_hi, n_keep) {
    if (!is.finite(n_keep)) n_keep <- .Machine$integer.max
    if (n_keep <= 0L) return(NULL)
    lead <- NULL
    if (!is.null(loci_src) && nrow(loci_src) > 0L) {
      lead <- as.data.table(loci_src)[
        as.character(cell_type) == cell_label &
          lead_p < p_hi & lead_p >= p_lo
      ]
      if (nrow(lead) > 0L) {
        setorder(lead, lead_p)
        if (n_keep < nrow(lead)) lead <- head(lead, n_keep)
        lead <- lead[, .(chr = as.integer(chr), pos = as.integer(lead_pos),
                         p = as.numeric(lead_p))]
      }
    }
    if (is.null(lead) || nrow(lead) == 0L) {
      cand <- dt[p < p_hi & p >= p_lo]
      if (nrow(cand) > 0L) {
        setorder(cand, p)
        lead <- cand[, .SD[greedy_clump(.SD, window = CLUMP_WINDOW)], by = chr
                     ][order(p)]
        if (n_keep < nrow(lead)) lead <- lead[seq_len(n_keep)]
        lead <- lead[, .(chr, pos, p)]
      }
    }
    if (is.null(lead) || nrow(lead) == 0L) return(NULL)
    lead[, gene := nearest_gene_labels(chr, pos, p)]
    lead <- merge(lead, off[, .(chr, offset)], by = "chr")
    lead[, x := as.numeric(pos) + offset]
    lead <- merge(lead, dt[, .(chr, pos, signed_y)], by = c("chr", "pos"),
                  all.x = TRUE)
    n_drop <- sum(is.na(lead$signed_y))
    if (n_drop > 0L)
      cat("    warning:", n_drop, "lead(s) not found in Miami SNPs; dropped\n")
    lead <- lead[!is.na(signed_y)]
    if (nrow(lead) == 0L) NULL else lead
  }

  lead_gws  <- pick_leads(loci_dt,      p_lo = 0,        p_hi = PVAL_GWS,  n_keep = n_label)
  lead_sugg <- pick_leads(loci_sugg_dt, p_lo = PVAL_GWS, p_hi = PVAL_SUGG, n_keep = n_label_sugg)

  if (!is.null(lead_gws))
    cat("    labelled GWS genes:", paste(lead_gws$gene, collapse = ", "), "\n")
  if (!is.null(lead_sugg))
    cat("    labelled suggestive genes:", paste(lead_sugg$gene, collapse = ", "), "\n")

  # Layout all labels together so GWS and suggestive never collide.
  # Collapse identical nearest-gene × chromosome × effect-side to one label
  # (annotate ×N when several independent leads share a gene) — this is what
  # made panels unreadable (e.g. CNBD1/IZUMO3/BLM repeated 4–5×).
  lab_all <- rbindlist(list(
    if (!is.null(lead_gws))  cbind(lead_gws,  tier = "gws")  else NULL,
    if (!is.null(lead_sugg)) cbind(lead_sugg, tier = "sugg") else NULL
  ), fill = TRUE)
  y_data <- max(thr_gws * 1.15, max(abs(dt$signed_y), na.rm = TRUE) * 1.12)
  y_ceiling <- max(y_data * 1.8, y_data + 5)
  x_range <- range(dt$x, na.rm = TRUE)
  if (nrow(lab_all) > 0L) {
    lab_all[, side := ifelse(signed_y >= 0, 1L, -1L)]
    lab_all[, tier_ord := ifelse(tier == "gws", 0L, 1L)]
    lab_all[, gene_key := gene]
    setorder(lab_all, tier_ord, p)
    n_before <- nrow(lab_all)
    lab_all <- lab_all[, {
      n <- .N
      g <- gene_key[1L]
      .(pos = pos[1L], x = x[1L], signed_y = signed_y[1L],
        p = p[1L], tier = tier[1L],
        gene = if (n > 1L) paste0(g, "\u00d7", n) else g)
    }, by = .(chr, gene_key, side)]
    cat("    label collapse:", n_before, "leads ->", nrow(lab_all),
        "unique gene\u00d7chr\u00d7side labels\n")
    lab_all <- layout_gene_labels(lab_all, x_range = x_range, y_ceiling = y_ceiling)
    y_max <- max(y_data, max(abs(lab_all$y_lab), na.rm = TRUE) + 0.9)
  } else {
    y_max <- y_data
  }

  p <- ggplot(dt, aes(x = x, y = signed_y, colour = pt_col)) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
    geom_hline(yintercept = c(thr_gws, -thr_gws),
               linetype = "dashed", colour = "grey35", linewidth = 0.4) +
    geom_hline(yintercept = c(thr_sugg, -thr_sugg),
               linetype = "dotted", colour = "grey55", linewidth = 0.35) +
    # Raster layer: hide legend (ggrastr keeps legend dots tiny).
    raster_points(geom_point(size = 0.55, alpha = 0.85, show.legend = FALSE)) +
    scale_colour_manual(
      values = c(ns_odd = "grey82", ns_even = COL_NS,
                 sugg = COL_SUGG, gws = COL_GWS),
      breaks = c("gws", "sugg", "ns_even"),
      guide = "none"
    ) +
    scale_x_continuous(breaks = ticks$center, labels = ticks$chr,
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(
      limits = c(-y_max, y_max),
      labels = function(v) format(abs(v), trim = TRUE),
      breaks = scales::pretty_breaks(n = 7),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(title = NULL, tag = panel_tag,
         x = "Chromosome",
         y = expression(-log[10](italic(p)))) +
    theme_fig3(base_size = 13) +
    theme(axis.text.x  = element_text(size = 10),
          axis.text.y  = element_text(size = 11),
          plot.margin  = margin(12, 8, 6, 12),
          plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
          plot.tag.position = "topleft",
          # Legend drawn separately below (ggrastr + tiny points break in-plot legend)
          legend.position = "none")

  p <- p +
    annotate("text", x = -Inf, y = Inf, label = "higher proportion",
             hjust = -0.04, vjust = 1.4, size = 3.6, colour = "grey35") +
    annotate("text", x = -Inf, y = -Inf, label = "lower proportion",
             hjust = -0.04, vjust = -0.9, size = 3.6, colour = "grey35") +
    annotate("text", x = Inf, y = Inf, label = cell_label,
             hjust = 1.05, vjust = 1.4, size = 4.5,
             fontface = "bold", colour = "grey20")

  # ggrepel only sees the ~30–50 label rows (not 15M GWAS points), so it is
  # safe/fast. Pre-layout provides the initial nudge; repulsion finishes the job.
  add_repel <- function(plot, lab, face, size, colour, seg_col) {
    if (is.null(lab) || nrow(lab) == 0L) return(plot)
    plot + ggrepel::geom_text_repel(
      data = lab,
      aes(x = x, y = signed_y, label = gene),
      inherit.aes = FALSE,
      size = size,
      fontface = face,
      colour = colour,
      nudge_x = lab$x_lab - lab$x,
      nudge_y = lab$y_lab - lab$signed_y,
      box.padding = 0.45,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.size = 0.25,
      segment.colour = seg_col,
      force = 3.5,
      force_pull = 0.4,
      max.overlaps = Inf,
      max.time = 12,
      max.iter = 40000,
      seed = 1,
      xlim = c(-Inf, Inf),
      ylim = c(-y_max * 0.98, y_max * 0.98)
    )
  }
  if (nrow(lab_all) > 0L) {
    p <- add_repel(p, lab_all[tier == "sugg"],
                   face = "italic", size = 3.1,
                   colour = "grey30", seg_col = "grey60")
    p <- add_repel(p, lab_all[tier == "gws"],
                   face = "bold.italic", size = 3.8,
                   colour = "grey10", seg_col = "grey40")
  }

  p
}

# =============================================================================
# SECTION 11 — Panels C–E: Focus cell-type Manhattan + QQ rows
# =============================================================================

focus_label <- function(key) {
  if (!is.null(focus_matched[[key]])) focus_matched[[key]] else key
}

cat("--- Building Panel B (VIP Miami) ---\n")
miami_vip <- make_miami_plot(
  dt_ct        = get_focus_gwas("VIP"),
  cell_label   = focus_label("VIP"),
  loci_dt      = loci,
  loci_sugg_dt = loci_sugg,
  panel_tag    = "B"
)

cat("--- Building Panel C (L5.6.IT.Car3 Miami) ---\n")
miami_car3 <- make_miami_plot(
  dt_ct        = get_focus_gwas("L5.6.IT.Car3"),
  cell_label   = focus_label("L5.6.IT.Car3"),
  loci_dt      = loci,
  loci_sugg_dt = loci_sugg,
  panel_tag    = "C"
)

cat("--- Building Panel D (Microglia Miami) ---\n")
miami_micro <- make_miami_plot(
  dt_ct        = get_focus_gwas("Microglia"),
  cell_label   = focus_label("Microglia"),
  loci_dt      = loci,
  loci_sugg_dt = loci_sugg,
  panel_tag    = "D"
)

# Shared Miami colour legend as a real plot (cowplot::get_legend returns
# zeroGrob under ggplot2 >= 3.5). Large dots + explicit text, no guide box.
miami_legend_df <- data.frame(
  key = factor(c("gws", "sugg", "ns"), levels = c("gws", "sugg", "ns")),
  lab = c("Genome-wide (p < 5e-8)", "Suggestive (p < 1e-5)", "Not significant"),
  x   = c(1.0, 3.6, 6.0),
  y   = 1
)
miami_legend <- ggplot(miami_legend_df, aes(x = x, y = y)) +
  geom_point(aes(colour = key), size = 8.5) +
  geom_text(aes(label = lab), hjust = 0, nudge_x = 0.30,
            size = 4.4, colour = "grey20") +
  scale_colour_manual(
    values = c(gws = COL_GWS, sugg = COL_SUGG, ns = COL_NS),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0.7, 8.2), ylim = c(0.5, 1.5), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(4, 10, 8, 10))

# ---- Supplementary QQ plots (same three focus cell types) ----
cat("--- Building supplementary QQ plots ---\n")
qq_vip <- make_qq_plot(get_focus_gwas("VIP"), qc_summary,
                       focus_label("VIP"), panel_letter = "B")
qq_car3 <- make_qq_plot(get_focus_gwas("L5.6.IT.Car3"), qc_summary,
                        focus_label("L5.6.IT.Car3"), panel_letter = "C")
qq_micro <- make_qq_plot(get_focus_gwas("Microglia"), qc_summary,
                         focus_label("Microglia"), panel_letter = "D")

# =============================================================================
# SECTION 12 — Export deduplicated GWS loci (supplementary table, not in figure)
# =============================================================================
cat("--- Exporting deduplicated GWS loci table ---\n")

if (!is.null(loci) && nrow(loci) > 0L) {
  loci_export <- copy(loci)
  setorder(loci_export, chr, lead_pos)
  loci_export[, region_id := {
    rid <- integer(.N)
    cur <- 0L
    last_chr <- NA_integer_
    last_pos <- -Inf
    for (i in seq_len(.N)) {
      if (is.na(loci_export$chr[i]) || loci_export$chr[i] != last_chr ||
          loci_export$lead_pos[i] - last_pos > CLUMP_WINDOW) {
        cur <- cur + 1L
        last_chr <- loci_export$chr[i]
      }
      rid[i] <- cur
      last_pos <- loci_export$lead_pos[i]
    }
    rid
  }]
  loci_dedup <- loci_export[, .(
    lead_pos  = lead_pos[which.min(lead_p)],
    lead_p    = min(lead_p),
    lead_beta = lead_beta[which.min(lead_p)],
    lead_snp  = lead_snp[which.min(lead_p)],
    n_cell_types = .N,
    cell_types   = paste(sort(unique(as.character(cell_type))), collapse = "; ")
  ), by = .(chr, region_id)]
  setorder(loci_dedup, chr, lead_pos)
  loci_dedup[, locus_label := paste0("chr", chr, ":", lead_pos)]
  data.table::fwrite(loci_dedup,
                     file.path(OUT_DIR, "figure3_gws_loci_deduplicated.tsv"),
                     sep = "\t")
  cat("  Deduplicated GWS regions:", nrow(loci_dedup),
      "(from", nrow(loci), "cell-type-specific loci)\n")
}

# =============================================================================
# SECTION 13 — Assemble figure
# =============================================================================
cat("\n--- Assembling figure ---\n")

cat("  Panel grobs built OK.\n")

# Main figure: Panel A already includes its tag+legend; Miami tags via plot.tag;
# shared Miami colour legend is a separate large-dot plot under panel D.
miami_stack <- cowplot::plot_grid(
  miami_vip, miami_car3, miami_micro,
  ncol = 1,
  rel_heights = c(1.20, 1.20, 1.20)
)
fig3 <- cowplot::plot_grid(
  pA, miami_stack, miami_legend,
  ncol        = 1,
  rel_heights = c(1.55, 3.35, 0.42)
)

# Supplementary diagnostics figure: lambda GC across all cell types, then the
# three focus-cell-type QQ plots side by side.
supp_qq_row <- cowplot::plot_grid(
  qq_vip, qq_car3, qq_micro,
  nrow = 1, align = "h", axis = "tb"
)

fig3_supp <- cowplot::plot_grid(
  ggplotGrob(pB + theme(plot.margin = margin(6, 8, 6, 8))),
  supp_qq_row,
  ncol        = 1,
  rel_heights = c(1.5, 1.0)
)

# =============================================================================
# SECTION 14 — Save outputs
# =============================================================================
cat("--- Saving figure ---\n")

out_png <- file.path(OUT_DIR, "figure3_gwas_discovery.png")
out_svg <- file.path(OUT_DIR, "figure3_gwas_discovery.svg")

cowplot::save_plot(out_png, fig3,
                  base_width  = FIG_WIDTH,
                  base_height = FIG_HEIGHT,
                  dpi         = FIG_DPI,
                  bg          = "white")
cat("  Wrote PNG:", out_png, "\n")

# SVG optional — full Miami panels are slow to vectorise; set WRITE_SVG=1 to enable.
WRITE_SVG <- identical(Sys.getenv("WRITE_SVG", unset = "0"), "1")
HAS_SVG <- WRITE_SVG && requireNamespace("svglite", quietly = TRUE)
if (HAS_SVG) {
  cowplot::save_plot(out_svg, fig3,
                    base_width  = FIG_WIDTH,
                    base_height = FIG_HEIGHT,
                    bg          = "white")
  cat("  Wrote SVG:", out_svg, "\n")
} else {
  cat("  Skipping SVG (set WRITE_SVG=1 to enable).\n")
}

# --- Supplementary diagnostics figure (PNG only) ---
supp_png <- file.path(OUT_DIR, "figure3_supp_gwas_diagnostics.png")

cowplot::save_plot(supp_png, fig3_supp,
                  base_width  = SUPP_WIDTH,
                  base_height = SUPP_HEIGHT,
                  dpi         = FIG_DPI,
                  bg          = "white")
cat("  Wrote supp PNG:", supp_png, "\n")

cat("\n=== Figure 3 outputs ===\n")
cat("  Main PNG :", out_png, "\n")
if (HAS_SVG) cat("  Main SVG :", out_svg, "\n")
cat("  Supp PNG :", supp_png, "\n")
cat("  Panel A table:", file.path(OUT_DIR, "figure3_panelA_signal_burden.tsv"), "\n")
cat("  Lambda GC table:", file.path(OUT_DIR, "figure3_panelB_lambda_gc.tsv"), "\n")
cat("  GWS loci table:", file.path(OUT_DIR, "figure3_gws_loci_deduplicated.tsv"), "\n")
cat("Done.\n")
