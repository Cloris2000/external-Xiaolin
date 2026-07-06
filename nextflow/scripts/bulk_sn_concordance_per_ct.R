#!/usr/bin/env Rscript
# For one bulk cell type:
#   1. Read bulk meta-analysis tbl, filter to P < 1e-5
#   2. Greedy distance-based clumping (500 kb window) to get independent lead SNPs
#   3. For each lead SNP, look it up in the matched sn meta-analysis tbl
#   4. Classify: same_direction / opposite_direction / not_found (sn P < 1e-5 threshold)
#   5. Write per-locus TSV to results/meta_analysis_15cohorts/concordance/<ct>_concordance.tsv

suppressPackageStartupMessages({
  library(data.table)
})

# ---------- args -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
idx  <- which(args == "--celltype")
if (length(idx) == 0 || idx == length(args))
  stop("Usage: Rscript bulk_sn_concordance_per_ct.R --celltype <CellType>")
ct <- args[idx + 1L]

# ---------- cell type → sn prefix map ----------------------------------
ct_map <- c(
  Astrocyte     = "ASTRO",
  Endothelial   = "ENDO",
  "L5.6.NP"     = "L56NP",
  "L5.ET"       = "L5ET",
  "L6.CT"       = "L6CT",
  L6b           = "L6B",
  LAMP5         = "LAMP5LHX6",
  OPC           = "OPC",
  Oligodendrocyte = "OLIGO",
  PVALB         = "PVALB",
  SST           = "SST",
  VIP           = "VIP",
  VLMC          = "VLMC"
  # IT, L4.IT, L5.6.IT.Car3, Microglia, PAX6, Pericyte have no sn match
)

BULK_DIR <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts"
SN_DIR   <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_sn_rosmap_msbb_hbcc_alias15"
OUT_DIR  <- file.path(BULK_DIR, "concordance")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BULK_THRESH  <- 1e-5
CLUMP_WINDOW <- 500e3   # 500 kb

# ---------- find bulk tbl ----------------------------------------------
bulk_files <- list.files(BULK_DIR, pattern = paste0("^", ct, "_meta_analysis.*1\\.tbl$"),
                         full.names = TRUE)
if (length(bulk_files) == 0) stop("No bulk tbl found for cell type: ", ct)
bulk_file <- bulk_files[1]
message("[", ct, "] Reading bulk: ", basename(bulk_file))

bulk <- fread(bulk_file,
              select  = c("MarkerName", "Effect", "StdErr", "P-value"),
              nThread = 4,
              verbose = FALSE)
setnames(bulk, "P-value", "Pvalue")

# filter to suggestive threshold
bulk <- bulk[Pvalue < BULK_THRESH]
message("[", ct, "] Suggestive SNPs (P<1e-5): ", nrow(bulk))

# ---------- parse chromosome and position ------------------------------
# MarkerName format: chr7:12252119:G:GT
bulk[, c("CHR", "POS") := {
  parts <- tstrsplit(MarkerName, ":", fixed = TRUE)
  list(parts[[1]], as.integer(parts[[2]]))
}]

# ---------- greedy distance-based clumping ----------------------------
setorder(bulk, Pvalue)  # most significant first

leads <- character(0)
used  <- logical(nrow(bulk))

for (i in seq_len(nrow(bulk))) {
  if (used[i]) next
  leads <- c(leads, bulk$MarkerName[i])
  # mark all SNPs within 500kb on same chromosome as used
  same_chr <- bulk$CHR == bulk$CHR[i]
  near     <- abs(bulk$POS - bulk$POS[i]) <= CLUMP_WINDOW
  used[same_chr & near] <- TRUE
}

message("[", ct, "] Independent loci after clumping: ", length(leads))
lead_dt <- bulk[MarkerName %in% leads, .(MarkerName, CHR, POS, bulk_Effect = Effect, bulk_P = Pvalue)]

# ---------- find sn tbl ------------------------------------------------
sn_prefix <- ct_map[ct]
has_sn    <- !is.na(sn_prefix)

if (has_sn) {
  sn_files <- list.files(SN_DIR,
                         pattern = paste0("^", sn_prefix, "_meta_analysis.*1\\.tbl$"),
                         full.names = TRUE)
  has_sn <- length(sn_files) > 0
}

if (has_sn) {
  sn_file <- sn_files[1]
  message("[", ct, "] Reading sn: ", basename(sn_file))
  sn <- fread(sn_file,
              select  = c("MarkerName", "Effect", "P-value"),
              nThread = 4,
              verbose = FALSE)
  setnames(sn, "P-value", "sn_Pvalue")
  setnames(sn, "Effect",  "sn_Effect")

  # Parse sn positions
  sn[, c("sn_CHR", "sn_POS") := {
    parts <- tstrsplit(MarkerName, ":", fixed = TRUE)
    list(parts[[1]], as.integer(parts[[2]]))
  }]

  # Regional matching: for each bulk lead locus, find the best (lowest P)
  # sn SNP within ±500kb on the same chromosome.
  MATCH_WINDOW <- 500e3

  lead_dt[, c("sn_Effect", "sn_Pvalue", "sn_MarkerName") := {
    res <- mapply(function(chr, pos) {
      candidates <- sn[sn_CHR == chr & abs(sn_POS - pos) <= MATCH_WINDOW]
      if (nrow(candidates) == 0) return(list(NA_real_, NA_real_, NA_character_))
      best <- candidates[which.min(sn_Pvalue)]
      list(best$sn_Effect, best$sn_Pvalue, best$MarkerName)
    }, CHR, POS, SIMPLIFY = FALSE)
    list(
      vapply(res, `[[`, numeric(1), 1),
      vapply(res, `[[`, numeric(1), 2),
      vapply(res, `[[`, character(1), 3)
    )
  }]
} else {
  message("[", ct, "] No sn counterpart — all loci will be 'not_found'")
  lead_dt[, `:=`(sn_Effect = NA_real_, sn_Pvalue = NA_real_, sn_MarkerName = NA_character_)]
}

# ---------- classify ---------------------------------------------------
# No p-value threshold on sn (small sample size = weak power).
# Classify purely by direction of effect; "not_found" only when no
# sn proxy exists within ±500kb.
lead_dt[, category := fifelse(
  is.na(sn_Effect),
  "not_found",
  fifelse(
    sign(bulk_Effect) == sign(sn_Effect),
    "same_direction",
    "opposite_direction"
  )
)]

lead_dt[, bulk_sig := ifelse(bulk_P < 5e-8, "GW_significant", "suggestive")]
lead_dt[, celltype := ct]

# ---------- write ------------------------------------------------------
out_file <- file.path(OUT_DIR, paste0(ct, "_concordance.tsv"))
fwrite(lead_dt[, .(celltype, MarkerName, CHR, POS,
                   bulk_Effect, bulk_P, bulk_sig,
                   sn_MarkerName, sn_Effect, sn_Pvalue, category)],
       out_file, sep = "\t")
message("[", ct, "] Written: ", out_file)
