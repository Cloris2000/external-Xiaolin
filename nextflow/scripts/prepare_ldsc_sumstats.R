#!/usr/bin/env Rscript
# Prepare LDSC-ready summary statistics from a REGENIE .raw_p file.
#
# Maps chr:pos variants to rsIDs using the HapMap3 LD score reference,
# then outputs a tab-delimited file ready for munge_sumstats.py.
#
# Usage:
#   Rscript prepare_ldsc_sumstats.R \
#     --input  path/to/cohort_celltype_step2.regenie.raw_p \
#     --output path/to/output/cohort_celltype.ldsc_input.txt \
#     --ld_dir /external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr \
#     [--rsid_cache path/to/rsid_lookup.rds]
#
# Output columns: SNP A1 A2 BETA SE N P
#   SNP  = rsID (from HapMap3 LD score reference, matched by CHR+BP)
#   A1   = effect allele (ALLELE1 in Regenie)
#   A2   = other allele  (ALLELE0 in Regenie)
#   BETA = effect size
#   SE   = standard error
#   N    = sample size
#   P    = p-value

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--input",      type="character", help="Input .regenie.raw_p file"),
  make_option("--output",     type="character", help="Output LDSC-ready .txt file"),
  make_option("--ld_dir",     type="character",
              default="/external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr",
              help="Directory containing {chr}.l2.ldscore.gz LD score files"),
  make_option("--rsid_cache", type="character", default=NULL,
              help="Optional: path to cached rsID lookup RDS (speeds up repeated runs)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input))  stop("--input is required")
if (is.null(opt$output)) stop("--output is required")
if (!file.exists(opt$input)) stop(paste("Input not found:", opt$input))

# ── Build or load rsID lookup table (CHR + BP → rsID) ────────────────────────
cache_path <- opt$rsid_cache
if (!is.null(cache_path) && file.exists(cache_path)) {
  cat("Loading cached rsID lookup from:", cache_path, "\n")
  rsid_lookup <- readRDS(cache_path)
} else {
  cat("Building rsID lookup from LD score files in:", opt$ld_dir, "\n")
  chrs <- 1:22
  rsid_list <- lapply(chrs, function(chr) {
    f <- file.path(opt$ld_dir, paste0(chr, ".l2.ldscore.gz"))
    if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
    fread(f, select=c("CHR","SNP","BP"), showProgress=FALSE) %>%
      rename(CHROM=CHR, GENPOS=BP, rsID=SNP) %>%
      mutate(CHROM=as.integer(CHROM), GENPOS=as.integer(GENPOS))
  })
  rsid_lookup <- bind_rows(rsid_list)
  cat("  ->", nrow(rsid_lookup), "HapMap3 SNPs loaded\n")
  if (!is.null(cache_path)) {
    saveRDS(rsid_lookup, cache_path)
    cat("  Cached to:", cache_path, "\n")
  }
}

# ── Read raw_p file ───────────────────────────────────────────────────────────
cat("Reading:", opt$input, "\n")
gwas <- fread(opt$input, header=TRUE) %>%
  as.data.frame() %>%
  mutate(
    CHROM  = as.integer(gsub("chr", "", as.character(CHROM))),
    GENPOS = as.integer(GENPOS),
    BETA   = as.numeric(BETA),
    SE     = as.numeric(SE),
    N      = as.numeric(N),
    P      = as.numeric(P)
  ) %>%
  filter(!is.na(CHROM), !is.na(GENPOS), !is.na(BETA), !is.na(SE),
         !is.na(P), P > 0, P <= 1, SE > 0)

cat("  ->", nrow(gwas), "variants before rsID mapping\n")

# ── Map to rsIDs ──────────────────────────────────────────────────────────────
gwas_rsid <- gwas %>%
  inner_join(rsid_lookup, by=c("CHROM","GENPOS")) %>%
  filter(!is.na(rsID), rsID != ".") %>%
  # Keep one entry per rsID (deduplicate multi-allelic)
  group_by(rsID) %>%
  slice_min(order_by=P, n=1, with_ties=FALSE) %>%
  ungroup()

cat("  ->", nrow(gwas_rsid), "variants retained after rsID mapping\n")

if (nrow(gwas_rsid) < 100000) {
  warning("Fewer than 100k variants retained — LDSC estimates may be unreliable")
}

# ── Write output ──────────────────────────────────────────────────────────────
out <- gwas_rsid %>%
  dplyr::select(
    SNP  = rsID,
    A1   = ALLELE1,   # effect allele in Regenie
    A2   = ALLELE0,   # other allele
    BETA = BETA,
    SE   = SE,
    N    = N,
    P    = P
  )

dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
fwrite(out, opt$output, sep="\t", quote=FALSE)
cat("Written:", opt$output, "(", nrow(out), "SNPs )\n")
