#!/usr/bin/env Rscript
# Compare sn-derived vs bulk MGP-derived GWAS results (ROSMAP or HBCC).
# Uses sn GWAS as ground truth to evaluate bulk concordance per cell type.
#
# Metrics per cell type:
#   - Spearman correlation of -log10(P) genome-wide
#   - Concordance at top 100 hits
#   - Beta sign agreement at sn suggestive SNPs (P < 1e-5)
#   - Lambda (genomic inflation) for each
#
# Outputs:
#   - sn_vs_bulk_summary.tsv
#   - {cell_type}_sn_vs_bulk_scatter.png  (one per cell type)
#   - sn_vs_bulk_correlation_barplot.png

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option("--sn_dir",      type="character", help="Dir with sn .regenie.raw_p files"),
  make_option("--bulk_dir",    type="character", help="Dir with bulk .regenie.raw_p files"),
  make_option("--output_dir",  type="character", default=".", help="Output directory"),
  make_option("--cohort",      type="character", default="ROSMAP",
              help="Legacy: sets both sn_cohort and bulk_cohort if they are not specified"),
  make_option("--sn_cohort",   type="character", default=NULL,
              help="Cohort prefix in sn filenames (e.g. NIMH_HBCC). Overrides --cohort."),
  make_option("--bulk_cohort", type="character", default=NULL,
              help="Cohort prefix in bulk filenames (e.g. NIMH_HBCC_1M). Overrides --cohort.")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$sn_dir) || is.null(opt$bulk_dir)) stop("--sn_dir and --bulk_dir are required")

# Resolve cohort prefixes
sn_cohort   <- if (!is.null(opt$sn_cohort))   opt$sn_cohort   else opt$cohort
bulk_cohort <- if (!is.null(opt$bulk_cohort)) opt$bulk_cohort else opt$cohort

cat("sn cohort prefix  :", sn_cohort, "\n")
cat("bulk cohort prefix:", bulk_cohort, "\n\n")

dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

# ── Cell type crosswalks per cohort ──────────────────────────────────────────
# ROSMAP: sn (Green_sn subclasses) → bulk (MGP 19 types)
rosmap_crosswalk <- list(
  "Astrocyte"       = "Astrocyte",
  "Endothelial"     = "Endothelial",
  "OPC"             = "OPC",
  "Oligodendrocyte" = "Oligodendrocyte",
  "VLMC"            = "VLMC",
  "L4.IT"           = "L4.IT",
  "L5.ET"           = "L5.ET",
  "L5.6.NP"         = "L5.6.NP",
  "L6.CT"           = "L6.CT",
  "L6b"             = "L6b",
  "Microglia.PVM"   = "Microglia",
  "Pvalb"           = "PVALB",
  "Sst"             = "SST",
  "Vip"             = "VIP",
  "Lamp5"           = "LAMP5",
  "Pax6"            = "PAX6"
  # No direct bulk match: Chandelier, L2.3.IT, L5.IT, L6.IT, L6.IT.Car3,
  #                       Lamp5.Lhx6, Sncg, Sst.Chodl
)

# HBCC: sn (psychAD subclasses) → bulk (MGP 19 types from NIMH_HBCC_1M)
hbcc_crosswalk <- list(
  "Astro"        = "Astrocyte",
  "Endo"         = "Endothelial",
  "OPC"          = "OPC",
  "Oligo"        = "Oligodendrocyte",
  "VLMC"         = "VLMC",
  "Micro"        = "Microglia",
  "PVM"          = "Microglia",      # PVM grouped with Microglia in bulk
  "IN_PVALB"     = "PVALB",
  "IN_SST"       = "SST",
  "IN_VIP"       = "VIP",
  "IN_LAMP5_RELN"= "LAMP5",
  "EN_L5_ET"     = "L5.ET",
  "EN_L5_6_NP"   = "L5.6.NP",
  "EN_L6_CT"     = "EN_L6_CT",      # no bulk match → will be skipped (file won't exist)
  "EN_L6B"       = "L6b"
  # No direct bulk match: Adaptive, EN_L2_3_IT, EN_L3_5_IT_1/2/3,
  #                       EN_L6_IT_1/2, IN_ADARB2, IN_LAMP5_LHX6,
  #                       IN_PVALB_CHC, PC, SMC
)

# Select crosswalk based on sn cohort prefix
if (grepl("HBCC", sn_cohort, ignore.case=TRUE)) {
  crosswalk <- hbcc_crosswalk
  cat("Using HBCC crosswalk (", length(crosswalk), "pairs)\n\n")
} else {
  crosswalk <- rosmap_crosswalk
  cat("Using ROSMAP crosswalk (", length(crosswalk), "pairs)\n\n")
}

# ── Helper: load raw_p file ───────────────────────────────────────────────────
load_gwas <- function(path, label) {
  if (!file.exists(path)) { cat("  MISSING:", path, "\n"); return(NULL) }
  cat("  Loading", label, ":", path, "\n")
  dt <- fread(path, header=TRUE) %>%
    as.data.frame() %>%
    rename(POS = GENPOS) %>%
    mutate(
      CHROM  = as.integer(gsub("chr", "", as.character(CHROM))),
      POS    = as.integer(POS),
      P      = as.numeric(P),
      BETA   = as.numeric(BETA),
      LOG10P = -log10(P)
    ) %>%
    filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1) %>%
    dplyr::select(CHROM, POS, ID, P, BETA, LOG10P)
  n_before <- nrow(dt)
  # Deduplicate multi-allelic sites: keep lowest P per position to avoid cartesian join explosion
  dt <- dt %>%
    group_by(CHROM, POS) %>%
    slice_min(order_by = P, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.frame()
  if (nrow(dt) < n_before)
    cat("    Deduped multi-allelic sites:", n_before - nrow(dt), "removed\n")
  cat("    ->", nrow(dt), "variants\n")
  dt
}

# ── Helper: genomic lambda ────────────────────────────────────────────────────
calc_lambda <- function(p) {
  chi2 <- qchisq(1 - p[!is.na(p) & p > 0 & p <= 1], 1)
  median(chi2, na.rm=TRUE) / qchisq(0.5, 1)
}

# ── Process each cell type pair ──────────────────────────────────────────────
summary_rows  <- list()
scatter_plots <- list()

# Pre-load existing summary rows so completed cell types can be skipped
existing_summary_file <- file.path(opt$output_dir, "sn_vs_bulk_summary.tsv")
if (file.exists(existing_summary_file)) {
  existing_df <- read.table(existing_summary_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  for (i in seq_len(nrow(existing_df))) {
    ct <- existing_df$sn_cell_type[i]
    summary_rows[[ct]] <- existing_df[i, , drop=FALSE]
  }
  cat("Loaded", nrow(existing_df), "existing results from", existing_summary_file, "\n")
}

for (sn_ct in names(crosswalk)) {
  tryCatch({
    bulk_ct <- crosswalk[[sn_ct]]

    # Skip if scatter plot already exists (cell type already processed)
    scatter_file <- file.path(opt$output_dir, paste0(sn_ct, "_sn_vs_bulk_scatter.png"))
    if (file.exists(scatter_file)) {
      cat("\n=== Skipping (already done):", sn_ct, "===\n")
      next
    }

    cat("\n=== Processing:", sn_ct, "(sn) vs", bulk_ct, "(bulk) ===\n")

    sn_file   <- file.path(opt$sn_dir,   paste0(sn_cohort,   "_", sn_ct,   "_step2.regenie.raw_p"))
    bulk_file <- file.path(opt$bulk_dir, paste0(bulk_cohort, "_", bulk_ct, "_step2.regenie.raw_p"))

    sn   <- load_gwas(sn_file,   "sn")
    bulk <- load_gwas(bulk_file, "bulk")
    if (is.null(sn) || is.null(bulk)) {
      cat("  Skipping - file(s) missing\n"); next
    }

    merged <- inner_join(sn, bulk, by=c("CHROM","POS"), suffix=c(".sn",".bulk")) %>%
      filter(!is.na(LOG10P.sn), !is.na(LOG10P.bulk))
    cat("  Joined SNPs:", nrow(merged), "\n")
    if (nrow(merged) < 100) { cat("  Too few shared SNPs, skipping\n"); next }

    spearman_r         <- cor(merged$LOG10P.sn, merged$LOG10P.bulk, method="spearman")
    lambda_sn          <- calc_lambda(merged$P.sn)
    lambda_bulk        <- calc_lambda(merged$P.bulk)
    top100_sn          <- merged %>% arrange(P.sn)   %>% head(100) %>% pull(POS)
    top100_bulk        <- merged %>% arrange(P.bulk) %>% head(100) %>% pull(POS)
    concordance_top100 <- mean(top100_sn %in% top100_bulk)

    suggestive   <- merged %>% filter(P.sn < 1e-5)
    n_suggestive <- nrow(suggestive)
    sign_agree   <- if (n_suggestive > 0)
      mean(sign(suggestive$BETA.sn) == sign(suggestive$BETA.bulk), na.rm=TRUE)
    else NA_real_

    cat(sprintf("  Spearman r=%.4f | lambda_sn=%.3f | lambda_bulk=%.3f | top100=%.2f | sign_agree=%s (n_sugg=%d)\n",
                spearman_r, lambda_sn, lambda_bulk, concordance_top100,
                ifelse(is.na(sign_agree), "NA", sprintf("%.2f", sign_agree)), n_suggestive))

    summary_rows[[sn_ct]] <- data.frame(
      sn_cell_type       = sn_ct,
      bulk_cell_type     = bulk_ct,
      n_shared_snps      = nrow(merged),
      spearman_r         = round(spearman_r, 4),
      lambda_sn          = round(lambda_sn,  4),
      lambda_bulk        = round(lambda_bulk, 4),
      concordance_top100 = round(concordance_top100, 3),
      n_suggestive_sn    = n_suggestive,
      beta_sign_agree    = round(sign_agree, 3),
      stringsAsFactors   = FALSE
    )

    plot_data <- merged %>%
      mutate(category = case_when(
        P.sn < 5e-8 ~ "Genome-wide sig (sn)",
        P.sn < 1e-5 ~ "Suggestive (sn)",
        TRUE        ~ "Other"
      )) %>%
      arrange(desc(category))

    p_scatter <- ggplot(plot_data, aes(x=LOG10P.sn, y=LOG10P.bulk, color=category)) +
      geom_point(size=0.4, alpha=0.5) +
      scale_color_manual(values=c(
        "Other"                = "grey70",
        "Suggestive (sn)"      = "#2196F3",
        "Genome-wide sig (sn)" = "#F44336"
      )) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="black", linewidth=0.4) +
      labs(
        title    = paste0(sn_ct, " (sn) vs ", bulk_ct, " (bulk)"),
        subtitle = sprintf("Spearman r = %.3f | top-100 concordance = %.0f%%",
                           spearman_r, concordance_top100 * 100),
        x     = expression(-log[10](P)~"sn"),
        y     = expression(-log[10](P)~"bulk"),
        color = NULL
      ) +
      theme_cowplot(font_size=10) +
      theme(legend.position="bottom")

    scatter_file <- file.path(opt$output_dir, paste0(sn_ct, "_sn_vs_bulk_scatter.png"))
    ggsave(scatter_file, plot=p_scatter, width=6, height=5, dpi=200)
    cat("  Saved scatter:", scatter_file, "\n")
    scatter_plots[[sn_ct]] <- p_scatter

  }, error = function(e) {
    cat(paste0("\nERROR processing ", sn_ct, ": ", e$message, "\n"))
  })
}

# ── Summary table ─────────────────────────────────────────────────────────────
if (length(summary_rows) == 0) stop("No cell types were processed successfully.")

summary_df <- bind_rows(summary_rows) %>% arrange(desc(spearman_r))

summary_file <- file.path(opt$output_dir, "sn_vs_bulk_summary.tsv")
write.table(summary_df, summary_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("\nSummary table saved:", summary_file, "\n")
print(summary_df[, c("sn_cell_type","spearman_r","concordance_top100","beta_sign_agree","lambda_sn","lambda_bulk")])

# ── Correlation barplot ───────────────────────────────────────────────────────
p_bar <- ggplot(summary_df,
                aes(x=reorder(sn_cell_type, spearman_r), y=spearman_r, fill=spearman_r)) +
  geom_col() +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_fill_gradient2(low="#d73027", mid="#fee090", high="#1a9850", midpoint=0.1,
                       name="Spearman r") +
  coord_flip() +
  labs(
    title    = paste0("sn vs Bulk GWAS Concordance (", sn_cohort, ")"),
    subtitle = "Spearman correlation of -log10(P), genome-wide",
    x = "Cell type (sn name)",
    y = "Spearman r"
  ) +
  theme_cowplot(font_size=11)

bar_file <- file.path(opt$output_dir, "sn_vs_bulk_correlation_barplot.png")
ggsave(bar_file, plot=p_bar, width=7, height=6, dpi=200)
cat("Barplot saved:", bar_file, "\n")

cat("\nDone. All outputs in:", opt$output_dir, "\n")
