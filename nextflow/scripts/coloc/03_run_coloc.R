#!/usr/bin/env Rscript
# Run colocalization analysis between cell-type GWAS loci and disease GWAS.
# Uses coloc.abf (primary) and coloc.susie (for top loci if susieR available).
#
# For each locus in the cell-type GWAS and each disease GWAS:
#   1. Subset disease GWAS to the locus window
#   2. Align alleles between datasets (flip beta if REF/ALT swapped)
#   3. Run coloc.abf
#   4. Optionally run coloc.susie for loci with PP.H4 > 0.5
#   5. Save results
#
# Usage:
#   Rscript 03_run_coloc.R \
#     --cell_type VIP \
#     --loci_file results/coloc/loci/VIP_loci.tsv \
#     --locus_data results/coloc/loci/VIP_locus_data.tsv.gz \
#     --disease_dir results/coloc/disease_gwas \
#     --output_dir results/coloc/coloc_results \
#     --cell_type_N 2812 \
#     [--run_susie TRUE]

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(coloc)
})

option_list <- list(
  make_option("--cell_type",    type = "character", help = "Cell type name (e.g. VIP)"),
  make_option("--loci_file",    type = "character", help = "Lead loci TSV from 01_extract_loci.R"),
  make_option("--locus_data",   type = "character", help = "Locus data TSV.gz from 01_extract_loci.R"),
  make_option("--disease_dir",  type = "character", default = "results/coloc/disease_gwas",
              help = "Directory with standardized disease GWAS files"),
  make_option("--output_dir",   type = "character", default = "results/coloc/coloc_results",
              help = "Output directory for coloc results"),
  make_option("--cell_type_N",  type = "integer",   default = 2812,
              help = "Total effective N for cell-type meta-analysis [default: 2812]"),
  make_option("--run_susie",    type = "logical",   default = TRUE,
              help = "Run coloc.susie for loci with PP.H4 > 0.5 [default: TRUE]"),
  make_option("--pp_h4_susie",  type = "double",    default = 0.5,
              help = "PP.H4 threshold to trigger coloc.susie [default: 0.5]"),
  make_option("--min_snps",     type = "integer",   default = 10,
              help = "Minimum SNPs per locus/window to run coloc [default: 10]"),
  make_option("--disease_include", type = "character", default = NULL,
              help = "Comma-separated disease labels to include (basename without _hg19.tsv); default all")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  !is.null(opt$cell_type),
  !is.null(opt$loci_file),
  !is.null(opt$locus_data)
)

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Check susieR availability
susie_available <- opt$run_susie && requireNamespace("susieR", quietly = TRUE)
if (opt$run_susie && !susie_available) {
  cat("NOTE: susieR not installed — coloc.susie will be skipped.\n")
}

# ---------------------------------------------------------------------------
# Load cell-type loci
# ---------------------------------------------------------------------------
cat(sprintf("[%s] Loading loci from %s\n", opt$cell_type, opt$loci_file))
leads <- fread(opt$loci_file, header = TRUE, sep = "\t", data.table = FALSE)

if (nrow(leads) == 0) {
  cat(sprintf("[%s] No loci found. Exiting.\n", opt$cell_type))
  quit(status = 0)
}
cat(sprintf("[%s] %d independent loci to test\n", opt$cell_type, nrow(leads)))

cat(sprintf("[%s] Loading locus data from %s\n", opt$cell_type, opt$locus_data))
locus_data <- fread(opt$locus_data, header = TRUE, sep = "\t", data.table = FALSE)

# ---------------------------------------------------------------------------
# Discover disease GWAS files
# Expects: {disease_dir}/{Disease}/{Disease}_hg19.tsv
# ---------------------------------------------------------------------------
disease_files <- list.files(
  opt$disease_dir,
  pattern = "_hg19\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(disease_files) == 0) {
  stop(paste("No *_hg19.tsv disease files found in", opt$disease_dir,
             "\nRun 02_download_disease_gwas.sh first."))
}
if (!is.null(opt$disease_include) && nzchar(opt$disease_include)) {
  include_labels <- trimws(strsplit(opt$disease_include, ",", fixed = TRUE)[[1]])
  disease_labels <- sub("_hg19\\.tsv$", "", basename(disease_files))
  keep <- disease_labels %in% include_labels
  disease_files <- disease_files[keep]
  if (length(disease_files) == 0) {
    stop(paste("No disease files matched --disease_include:", opt$disease_include))
  }
  cat(sprintf("Filtered to %d disease(s) via --disease_include\n", length(disease_files)))
}
cat(sprintf("Found %d disease GWAS files:\n", length(disease_files)))
for (f in disease_files) cat("  ", basename(f), "\n")

# ---------------------------------------------------------------------------
# Allele harmonization helper
# Returns the dataset1 beta and the shared SNP list, accounting for
# REF/ALT flips between the two datasets (flip beta if needed).
# Drops strand-ambiguous SNPs (A/T and C/G).
# ---------------------------------------------------------------------------
STRAND_AMB <- list(c("A","T"), c("T","A"), c("C","G"), c("G","C"))
is_strand_ambiguous <- function(ref, alt) {
  paste(ref, alt) %in% c("A T", "T A", "C G", "G C")
}

harmonize_datasets <- function(d1, d2) {
  # d1, d2: data frames with columns snp, beta, se, ref, alt
  # Returns merged data frame with harmonized beta_d1, beta_d2, varbeta_d1, varbeta_d2
  m <- merge(d1, d2, by = "snp", suffixes = c("_ct", "_dis"))

  # Remove strand-ambiguous SNPs from both sides
  strand_amb <- is_strand_ambiguous(m$ref_ct, m$alt_ct)
  m <- m[!strand_amb, ]

  if (nrow(m) == 0) return(NULL)

  # Detect allele flips: when effect alleles are swapped, flip beta_dis
  flipped <- (toupper(m$ref_ct) == toupper(m$alt_dis)) &
             (toupper(m$alt_ct) == toupper(m$ref_dis))
  m$beta_dis[flipped] <- -m$beta_dis[flipped]

  # Remove mismatched alleles (neither matching nor flipped)
  matched <- (toupper(m$ref_ct) == toupper(m$ref_dis) &
              toupper(m$alt_ct) == toupper(m$alt_dis)) | flipped
  m <- m[matched, ]

  m
}

# ---------------------------------------------------------------------------
# Run coloc.abf for one (locus, disease) pair
# ---------------------------------------------------------------------------
run_coloc_abf_pair <- function(ct_window, dis_window, cell_type_N, disease_label,
                                locus_id, sdY = 1.0) {
  # ct_window:  data frame with snp, beta, se, ref, alt, varbeta (cell-type)
  # dis_window: data frame with snp, beta, se, ref, alt, varbeta, Ncase, Ncont

  m <- harmonize_datasets(ct_window, dis_window)
  if (is.null(m) || nrow(m) < 10) return(NULL)

  # Remove entries with missing or zero varbeta
  m <- m[is.finite(m$beta_ct) & is.finite(m$varbeta_ct) & m$varbeta_ct > 0 &
         is.finite(m$beta_dis) & is.finite(m$varbeta_dis) & m$varbeta_dis > 0, ]
  if (nrow(m) < 10) return(NULL)

  ds1 <- list(
    beta    = m$beta_ct,
    varbeta = m$varbeta_ct,
    N       = cell_type_N,
    sdY     = sdY,
    type    = "quant",
    snp     = m$snp
  )

  Ncase <- m$Ncase_dis[1]
  Ncont <- m$Ncont_dis[1]
  N_dis <- m$N_dis[1]
  s_frac <- if (!is.na(Ncase) && !is.na(Ncont) && Ncont > 0)
    Ncase / (Ncase + Ncont) else 0.5

  ds2 <- list(
    beta    = m$beta_dis,
    varbeta = m$varbeta_dis,
    N       = ifelse(!is.na(N_dis), N_dis, Ncase + Ncont),
    s       = s_frac,
    type    = "cc",
    snp     = m$snp
  )

  res <- tryCatch(
    coloc.abf(ds1, ds2),
    error = function(e) { cat("  coloc.abf error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res)) return(NULL)

  summ <- res$summary
  data.frame(
    locus_id       = locus_id,
    disease        = disease_label,
    n_snps         = nrow(m),
    PP.H0          = summ["PP.H0.abf"],
    PP.H1          = summ["PP.H1.abf"],
    PP.H2          = summ["PP.H2.abf"],
    PP.H3          = summ["PP.H3.abf"],
    PP.H4          = summ["PP.H4.abf"],
    coloc_method   = "coloc.abf",
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Run coloc.susie for one (locus, disease) pair
# ---------------------------------------------------------------------------
run_coloc_susie_pair <- function(ct_window, dis_window, cell_type_N,
                                  disease_label, locus_id, sdY = 1.0) {
  if (!requireNamespace("susieR", quietly = TRUE)) return(NULL)

  m <- harmonize_datasets(ct_window, dis_window)
  if (is.null(m) || nrow(m) < 50) return(NULL)
  m <- m[is.finite(m$beta_ct) & is.finite(m$varbeta_ct) & m$varbeta_ct > 0 &
         is.finite(m$beta_dis) & is.finite(m$varbeta_dis) & m$varbeta_dis > 0, ]
  if (nrow(m) < 50) return(NULL)

  Ncase <- m$Ncase_dis[1]
  Ncont <- m$Ncont_dis[1]
  N_dis <- m$N_dis[1]
  s_frac <- if (!is.na(Ncase) && !is.na(Ncont) && Ncont > 0)
    Ncase / (Ncase + Ncont) else 0.5

  ds1 <- list(beta = m$beta_ct,  varbeta = m$varbeta_ct,
               N = cell_type_N, sdY = sdY, type = "quant", snp = m$snp)
  ds2 <- list(beta = m$beta_dis, varbeta = m$varbeta_dis,
               N = ifelse(!is.na(N_dis), N_dis, Ncase + Ncont),
               s = s_frac, type = "cc", snp = m$snp)

  # Run SuSiE on each dataset separately, then coloc.susie
  s1 <- tryCatch(
    susieR::susie_rss(bhat = ds1$beta, shat = sqrt(ds1$varbeta),
                      n = ds1$N, var_y = ds1$sdY^2, L = 10,
                      estimate_residual_variance = TRUE),
    error = function(e) NULL
  )
  s2 <- tryCatch(
    susieR::susie_rss(bhat = ds2$beta, shat = sqrt(ds2$varbeta),
                      n = ds2$N, L = 10,
                      estimate_residual_variance = TRUE),
    error = function(e) NULL
  )
  if (is.null(s1) || is.null(s2)) return(NULL)

  res <- tryCatch(
    coloc.susie(ds1, ds2, susie.args = list(L = 10)),
    error = function(e) { cat("  coloc.susie error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(res) || is.null(res$summary) || nrow(res$summary) == 0) return(NULL)

  # Return the credible set pair with highest PP.H4
  best <- res$summary[which.max(res$summary$PP.H4.abf), , drop = FALSE]
  data.frame(
    locus_id       = locus_id,
    disease        = disease_label,
    n_snps         = nrow(m),
    PP.H0          = best$PP.H0.abf,
    PP.H1          = best$PP.H1.abf,
    PP.H2          = best$PP.H2.abf,
    PP.H3          = best$PP.H3.abf,
    PP.H4          = best$PP.H4.abf,
    coloc_method   = "coloc.susie",
    ncs1           = best$nsnps.cs1,
    ncs2           = best$nsnps.cs2,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Main loop: iterate over disease files × loci
# ---------------------------------------------------------------------------
all_results <- list()

for (dis_file in disease_files) {
  disease_label <- sub("_hg19\\.tsv$", "", basename(dis_file))
  cat(sprintf("\n[%s] Processing disease: %s\n", opt$cell_type, disease_label))

  dis_raw <- tryCatch(
    fread(dis_file, header = TRUE, sep = "\t", data.table = FALSE),
    error = function(e) { cat("  ERROR reading", dis_file, "\n"); NULL }
  )
  if (is.null(dis_raw) || nrow(dis_raw) == 0) next

  # Standardize disease data
  dis_raw <- dis_raw %>%
    mutate(
      beta    = as.numeric(beta),
      se      = as.numeric(se),
      p       = as.numeric(p),
      N       = suppressWarnings(as.numeric(N)),
      Ncase   = suppressWarnings(as.numeric(Ncase)),
      Ncont   = suppressWarnings(as.numeric(Ncont)),
      varbeta = se^2,
      chr     = as.character(chr),
      pos     = as.integer(pos)
    ) %>%
    filter(!is.na(beta), !is.na(se), se > 0, !is.na(pos))

  # Index disease data by chr for fast subset
  dis_by_chr <- split(dis_raw, dis_raw$chr)

  for (i in seq_len(nrow(leads))) {
    locus   <- leads[i, ]
    locus_id <- locus$locus_id
    locus_chr <- as.character(locus$chr)

    cat(sprintf("  Locus %d/%d: %s\n", i, nrow(leads), locus_id))

    # Cell-type data for this locus
    ct_window <- locus_data %>%
      filter(locus_id == locus$locus_id) %>%
      mutate(
        varbeta = as.numeric(se)^2,
        ref     = toupper(ref),
        alt     = toupper(alt)
      ) %>%
      filter(!is.na(beta), !is.na(varbeta), varbeta > 0) %>%
      rename(beta_ct = beta, se_ct = se, varbeta_ct = varbeta,
             ref_ct = ref, alt_ct = alt) %>%
      select(snp, beta_ct, se_ct, varbeta_ct, ref_ct, alt_ct)

    if (nrow(ct_window) < opt$min_snps) {
      cat(sprintf("    Skip: only %d CT SNPs (< %d)\n", nrow(ct_window), opt$min_snps))
      next
    }

    # Disease data for this locus window
    dis_chr_data <- dis_by_chr[[locus_chr]]
    if (is.null(dis_chr_data)) {
      cat(sprintf("    Skip: no disease data on chr%s\n", locus_chr))
      next
    }
    dis_window <- dis_chr_data %>%
      filter(pos >= locus$window_start, pos <= locus$window_end) %>%
      rename(beta_dis = beta, se_dis = se, varbeta_dis = varbeta,
             ref_dis = ref, alt_dis = alt,
             N_dis = N, Ncase_dis = Ncase, Ncont_dis = Ncont) %>%
      select(snp, beta_dis, se_dis, varbeta_dis, ref_dis, alt_dis,
             N_dis, Ncase_dis, Ncont_dis)

    if (nrow(dis_window) < opt$min_snps) {
      cat(sprintf("    Skip: only %d disease SNPs (< %d)\n",
                  nrow(dis_window), opt$min_snps))
      next
    }

    # --- coloc.abf ---
    abf_res <- run_coloc_abf_pair(
      ct_window   = ct_window,
      dis_window  = dis_window,
      cell_type_N = opt$cell_type_N,
      disease_label = disease_label,
      locus_id    = locus_id
    )

    if (!is.null(abf_res)) {
      abf_res$cell_type <- opt$cell_type
      all_results[[length(all_results) + 1]] <- abf_res
      cat(sprintf("    coloc.abf PP.H4 = %.3f\n", abf_res$PP.H4))

      # --- coloc.susie for strong hits ---
      if (susie_available && !is.na(abf_res$PP.H4) &&
          abf_res$PP.H4 >= opt$pp_h4_susie) {
        cat(sprintf("    Running coloc.susie (PP.H4 threshold met)...\n"))
        susie_res <- run_coloc_susie_pair(
          ct_window     = ct_window,
          dis_window    = dis_window,
          cell_type_N   = opt$cell_type_N,
          disease_label = disease_label,
          locus_id      = locus_id
        )
        if (!is.null(susie_res)) {
          susie_res$cell_type <- opt$cell_type
          all_results[[length(all_results) + 1]] <- susie_res
          cat(sprintf("    coloc.susie PP.H4 = %.3f\n", susie_res$PP.H4))
        }
      }
    }
  }  # end loci loop
}  # end disease loop

# ---------------------------------------------------------------------------
# Write results
# ---------------------------------------------------------------------------
if (length(all_results) == 0) {
  cat(sprintf("[%s] No coloc results produced.\n", opt$cell_type))
  quit(status = 0)
}

results_df <- bind_rows(all_results) %>%
  arrange(desc(PP.H4))

outfile <- file.path(opt$output_dir,
                     paste0(opt$cell_type, "_coloc_results.tsv"))
fwrite(results_df, outfile, sep = "\t", quote = FALSE)
cat(sprintf("\n[%s] Results written to %s\n", opt$cell_type, outfile))
cat(sprintf("[%s] Total pairs tested: %d\n", opt$cell_type, nrow(results_df)))
cat(sprintf("[%s] PP.H4 >= 0.8 (strong coloc): %d\n",
            opt$cell_type, sum(results_df$PP.H4 >= 0.8, na.rm = TRUE)))
cat(sprintf("[%s] PP.H4 >= 0.5 (moderate coloc): %d\n",
            opt$cell_type, sum(results_df$PP.H4 >= 0.5, na.rm = TRUE)))
