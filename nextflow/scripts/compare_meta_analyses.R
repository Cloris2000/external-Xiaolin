#!/usr/bin/env Rscript
# compare_meta_analyses.R
#
# Compares two meta-analysis runs:
#   A) meta_analysis_13cohorts_design_matrix   — tech-covariate-corrected, design= preserved bio signal
#   B) meta_analysis_no_tech_except_rosmap_libprep — no tech-covariate regression
#
# Outputs written to results/meta_analysis_comparison/:
#   lambda_comparison.png          — genomic inflation factor per cell type
#   n_hits_comparison.png          — GW + suggestive hit counts per cell type
#   ldsc_h2_comparison.png         — SNP h² per cell type (both runs)
#   ldsc_intercept_comparison.png  — LDSC intercept (confounding signal)
#   per_celltype_pvalue_scatter/   — -log10(P) scatter for each cell type
#   top_hits_comparison.tsv        — top GW hits unique to each run or shared

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(optparse)
})

option_list <- list(
  make_option("--dir_design", type="character",
    default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_13cohorts_design_matrix"),
  make_option("--dir_notech", type="character",
    default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_no_tech_except_rosmap_libprep"),
  make_option("--outdir", type="character",
    default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_comparison"),
  make_option("--gw_thresh",   type="double",  default=5e-8),
  make_option("--sugg_thresh", type="double",  default=1e-5)
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir, "per_celltype_scatter"), showWarnings=FALSE)

label_design <- "With tech-cov correction\n(design= protects bio signal)"
label_notech <- "No tech-cov correction"

# ── Helper ────────────────────────────────────────────────────────────────────
read_tbl <- function(dir, ct) {
  f <- list.files(dir, pattern=paste0("^", ct, "_meta_analysis_.*\\.tbl$"),
                  full.names=TRUE)
  f <- f[!grepl("1\\.tbl$", f)]
  if (length(f) == 0) return(NULL)
  dt <- fread(f[1], showProgress=FALSE)
  # Standardise column names (METAL output)
  setnames(dt, old=c("MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","P-value","Direction"),
           new=c("SNP","A1","A2","Freq1","Beta","SE","P","Direction"), skip_absent=TRUE)
  if (!"P" %in% colnames(dt) && "Pvalue" %in% colnames(dt)) setnames(dt, "Pvalue", "P")
  dt[, cell_type := ct]
  dt
}

cell_types <- c("Astrocyte","Endothelial","IT","L4.IT","L5.6.IT.Car3","L5.6.NP",
                "L5.ET","L6.CT","L6b","LAMP5","Microglia","OPC","Oligodendrocyte",
                "PAX6","PVALB","Pericyte","SST","VIP","VLMC")

# ── 1. Load meta-analysis summary tables ──────────────────────────────────────
read_summary <- function(dir) {
  f <- file.path(dir, "plots", "meta_analysis_summary.tsv")
  if (!file.exists(f)) return(NULL)
  read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}
summ_d <- read_summary(opt$dir_design)
summ_n <- read_summary(opt$dir_notech)

# If summary is incomplete (only 1 row), rebuild from .tbl files
rebuild_summary <- function(dir, cts, gw_thresh, sugg_thresh) {
  cat("  Rebuilding summary from .tbl files for", basename(dir), "...\n")
  rows <- lapply(cts, function(ct) {
    dt <- read_tbl(dir, ct)
    if (is.null(dt) || !"P" %in% colnames(dt)) return(NULL)
    p <- as.numeric(dt$P)
    chisq <- qchisq(p, df=1, lower.tail=FALSE)
    lambda <- round(median(chisq, na.rm=TRUE) / qchisq(0.5, 1), 4)
    data.frame(cell_type=ct, n_variants=sum(!is.na(p)),
               lambda=lambda,
               n_gw=sum(p < gw_thresh, na.rm=TRUE),
               n_sugg=sum(p < sugg_thresh, na.rm=TRUE))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

if (is.null(summ_d) || nrow(summ_d) < 5)
  summ_d <- rebuild_summary(opt$dir_design, cell_types, opt$gw_thresh, opt$sugg_thresh)
if (is.null(summ_n) || nrow(summ_n) < 5)
  summ_n <- rebuild_summary(opt$dir_notech, cell_types, opt$gw_thresh, opt$sugg_thresh)

summ_both <- bind_rows(
  summ_d %>% mutate(run=label_design),
  summ_n %>% mutate(run=label_notech)
)
ct_order_lambda <- summ_d %>% arrange(lambda) %>% pull(cell_type)
summ_both$cell_type <- factor(summ_both$cell_type, levels=ct_order_lambda)
summ_both$run <- factor(summ_both$run, levels=c(label_notech, label_design))

# ── 2. Lambda plot ────────────────────────────────────────────────────────────
cat("Plotting lambda comparison...\n")
p_lambda <- ggplot(summ_both, aes(x=cell_type, y=lambda, fill=run)) +
  geom_col(position=position_dodge(0.75), width=0.7) +
  geom_hline(yintercept=1, linetype="dashed", colour="black", linewidth=0.4) +
  scale_fill_manual(values=c("#ef6548","#2171b5")) +
  labs(title="Genomic Inflation Factor (λ) per Cell Type",
       subtitle="λ=1 ideal; <1 conservative; >1 inflated",
       x=NULL, y="λ (genomic inflation)", fill=NULL) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=40, hjust=1),
        legend.position="top", panel.grid.major.x=element_blank())
ggsave(file.path(opt$outdir, "lambda_comparison.png"), p_lambda,
       width=13, height=6, dpi=300)
cat("  Saved: lambda_comparison.png\n")

# ── 3. Hit counts ─────────────────────────────────────────────────────────────
cat("Plotting hit count comparison...\n")
ct_order_hits <- summ_d %>% arrange(desc(n_gw)) %>% pull(cell_type)
summ_both$cell_type <- factor(summ_both$cell_type, levels=ct_order_hits)

hits_long <- summ_both %>%
  select(cell_type, run, n_gw, n_sugg) %>%
  pivot_longer(cols=c(n_gw, n_sugg), names_to="threshold",
               values_to="n_hits") %>%
  mutate(threshold=recode(threshold,
    n_gw   = paste0("Genome-wide (p<", opt$gw_thresh, ")"),
    n_sugg = paste0("Suggestive (p<", opt$sugg_thresh, ")")))

p_hits <- ggplot(hits_long, aes(x=cell_type, y=n_hits, fill=run)) +
  geom_col(position=position_dodge(0.75), width=0.7) +
  facet_wrap(~threshold, scales="free_y", ncol=1) +
  scale_fill_manual(values=c("#ef6548","#2171b5")) +
  labs(title="Number of Significant Hits per Cell Type",
       x=NULL, y="Number of variants", fill=NULL) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=40, hjust=1),
        legend.position="top", panel.grid.major.x=element_blank())
ggsave(file.path(opt$outdir, "n_hits_comparison.png"), p_hits,
       width=13, height=9, dpi=300)
cat("  Saved: n_hits_comparison.png\n")

# ── 4. LDSC h² comparison ─────────────────────────────────────────────────────
cat("Plotting LDSC h² comparison...\n")
read_ldsc <- function(dir) {
  f <- file.path(dir, "ldsc", "summary", "ldsc_h2_summary.tsv")
  if (!file.exists(f)) return(NULL)
  read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}
ldsc_d <- read_ldsc(opt$dir_design)
ldsc_n <- read_ldsc(opt$dir_notech)

if (!is.null(ldsc_d) && !is.null(ldsc_n)) {
  ldsc_both <- bind_rows(
    ldsc_d %>% mutate(run=label_design),
    ldsc_n %>% mutate(run=label_notech)
  )
  ct_order_h2 <- ldsc_d %>% arrange(desc(observed_h2)) %>% pull(cell_type)
  ldsc_both$cell_type <- factor(ldsc_both$cell_type, levels=ct_order_h2)
  ldsc_both$run <- factor(ldsc_both$run, levels=c(label_notech, label_design))

  p_h2 <- ggplot(ldsc_both, aes(x=cell_type, y=observed_h2, fill=run)) +
    geom_col(position=position_dodge(0.75), width=0.7) +
    geom_errorbar(aes(ymin=observed_h2 - observed_h2_se,
                      ymax=observed_h2 + observed_h2_se),
                  position=position_dodge(0.75), width=0.3, linewidth=0.4) +
    geom_hline(yintercept=0, linewidth=0.4) +
    scale_fill_manual(values=c("#ef6548","#2171b5")) +
    labs(title="SNP Heritability (h²) per Cell Type — LDSC",
         subtitle="Error bars = ±1 SE",
         x=NULL, y="Observed h²", fill=NULL) +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=40, hjust=1),
          legend.position="top", panel.grid.major.x=element_blank())
  ggsave(file.path(opt$outdir, "ldsc_h2_comparison.png"), p_h2,
         width=13, height=6, dpi=300)
  cat("  Saved: ldsc_h2_comparison.png\n")

  # LDSC intercept — elevated intercept signals confounding/stratification
  p_int <- ggplot(ldsc_both, aes(x=cell_type, y=intercept, fill=run)) +
    geom_col(position=position_dodge(0.75), width=0.7) +
    geom_errorbar(aes(ymin=intercept - intercept_se,
                      ymax=intercept + intercept_se),
                  position=position_dodge(0.75), width=0.3, linewidth=0.4) +
    geom_hline(yintercept=1, linetype="dashed", linewidth=0.4) +
    scale_fill_manual(values=c("#ef6548","#2171b5")) +
    labs(title="LDSC Intercept per Cell Type",
         subtitle="Intercept=1 means no confounding; >1 suggests inflation from structure/confounders",
         x=NULL, y="LDSC intercept", fill=NULL) +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=40, hjust=1),
          legend.position="top", panel.grid.major.x=element_blank())
  ggsave(file.path(opt$outdir, "ldsc_intercept_comparison.png"), p_int,
         width=13, height=6, dpi=300)
  cat("  Saved: ldsc_intercept_comparison.png\n")

  # Save LDSC comparison table
  write.table(ldsc_both, file.path(opt$outdir, "ldsc_comparison.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("  Saved: ldsc_comparison.tsv\n")
}

# ── 5. Top hits with heterogeneity scores ────────────────────────────────────
# Read only significant hits from each .tbl (fast — uses fread with row filtering).
# Extracts GW (p<5e-8) and suggestive (p<1e-5) hits with full heterogeneity columns.
cat("Extracting top hits with heterogeneity scores...\n")

het_cols <- c("MarkerName","Allele1","Allele2","Freq1","Effect","StdErr",
              "P-value","Direction","HetISq","HetChiSq","HetDf","HetPVal")

read_hits <- function(dir, ct, gw_thresh, sugg_thresh) {
  f <- list.files(dir, pattern=paste0("^", ct, "_meta_analysis_.*\\.tbl$"),
                  full.names=TRUE)
  f <- f[!grepl("1\\.tbl$", f)]
  if (length(f) == 0) return(NULL)
  dt <- fread(f[1], showProgress=FALSE, select=intersect(het_cols, fread(f[1], nrows=0) |> colnames()))
  # standardise P column name
  if ("P-value" %in% colnames(dt)) setnames(dt, "P-value", "P")
  dt[as.numeric(P) < sugg_thresh]
}

hits_design_list <- list()
hits_notech_list <- list()

for (ct in cell_types) {
  cat(" ", ct, "...")
  h_d <- read_hits(opt$dir_design, ct, opt$gw_thresh, opt$sugg_thresh)
  h_n <- read_hits(opt$dir_notech, ct, opt$gw_thresh, opt$sugg_thresh)

  if (!is.null(h_d) && nrow(h_d) > 0) {
    h_d[, `:=`(cell_type=ct, run="design",
               threshold=ifelse(as.numeric(P) < opt$gw_thresh, "genome-wide", "suggestive"))]
    hits_design_list[[ct]] <- h_d
  }
  if (!is.null(h_n) && nrow(h_n) > 0) {
    h_n[, `:=`(cell_type=ct, run="no_tech",
               threshold=ifelse(as.numeric(P) < opt$gw_thresh, "genome-wide", "suggestive"))]
    hits_notech_list[[ct]] <- h_n
  }
  cat(sprintf(" design=%d  no-tech=%d\n",
              if (!is.null(h_d)) nrow(h_d) else 0,
              if (!is.null(h_n)) nrow(h_n) else 0))
}

hits_design <- rbindlist(hits_design_list, fill=TRUE)
hits_notech <- rbindlist(hits_notech_list, fill=TRUE)
hits_all    <- rbindlist(list(hits_design, hits_notech), fill=TRUE)

# Save full tables
write.table(hits_all[run=="design"],
            file.path(opt$outdir, "top_hits_design_with_het.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(hits_all[run=="no_tech"],
            file.path(opt$outdir, "top_hits_notech_with_het.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("  Saved: top_hits_design_with_het.tsv\n")
cat("  Saved: top_hits_notech_with_het.tsv\n")

# Shared vs unique hits (match on MarkerName + cell_type)
if (nrow(hits_design) > 0 && nrow(hits_notech) > 0) {
  gw_d <- hits_design[threshold=="genome-wide", .(MarkerName, cell_type)]
  gw_n <- hits_notech[threshold=="genome-wide", .(MarkerName, cell_type)]
  shared   <- merge(gw_d, gw_n, by=c("MarkerName","cell_type"))
  only_d   <- gw_d[!gw_n, on=c("MarkerName","cell_type")]
  only_n   <- gw_n[!gw_d, on=c("MarkerName","cell_type")]
  cat(sprintf("  GW hits — shared: %d  design-only: %d  no-tech-only: %d\n",
              nrow(shared), nrow(only_d), nrow(only_n)))
}

# ── Het score summary plot for GW + suggestive hits ──────────────────────────
cat("Plotting heterogeneity at top hits...\n")

if (nrow(hits_all) > 0) {
  hits_plot <- hits_all[!is.na(HetISq)]
  hits_plot[, P_num := as.numeric(P)]
  hits_plot[, logP  := -log10(P_num)]
  hits_plot[, het_flag := ifelse(HetISq > 75, "High het (I²>75%)",
                          ifelse(HetISq > 50, "Moderate het (I²>50%)", "Low het (I²≤50%)"))]
  hits_plot[, het_flag := factor(het_flag,
              levels=c("Low het (I²≤50%)","Moderate het (I²>50%)","High het (I²>75%)"))]
  hits_plot[, run_label := ifelse(run=="design", label_design, label_notech)]

  # I² vs -log10(P) scatter coloured by run, faceted by threshold
  p_het_scatter <- ggplot(hits_plot, aes(x=logP, y=HetISq, colour=run_label,
                                          shape=het_flag)) +
    geom_point(alpha=0.6, size=1.5) +
    geom_hline(yintercept=50, linetype="dashed", colour="#fd8d3c", linewidth=0.4) +
    geom_hline(yintercept=75, linetype="dashed", colour="#e31a1c", linewidth=0.4) +
    facet_wrap(~threshold, scales="free_x") +
    scale_colour_manual(values=c("#2171b5","#ef6548")) +
    scale_shape_manual(values=c(16, 17, 4)) +
    labs(title="Heterogeneity (I²) at Top Hits — design= vs no-tech correction",
         subtitle="Dashed lines: I²=50% (orange), I²=75% (red)",
         x="-log10(P)", y="I² (%)", colour=NULL, shape=NULL) +
    theme_bw(base_size=12) +
    theme(legend.position="top")
  ggsave(file.path(opt$outdir, "het_at_top_hits_scatter.png"),
         p_het_scatter, width=13, height=6, dpi=300)
  cat("  Saved: het_at_top_hits_scatter.png\n")

  # I² distribution barplot by cell type and run
  het_summary_ct <- hits_plot[, .(
    n_hits      = .N,
    n_high_het  = sum(HetISq > 75, na.rm=TRUE),
    n_mod_het   = sum(HetISq > 50 & HetISq <= 75, na.rm=TRUE),
    pct_high_het= round(mean(HetISq > 75, na.rm=TRUE)*100, 1),
    median_I2   = round(median(HetISq, na.rm=TRUE), 1)
  ), by=.(cell_type, threshold, run_label)]

  ct_order_het <- het_summary_ct[threshold=="genome-wide" & run_label==label_design][
    order(desc(pct_high_het))]$cell_type
  if (length(ct_order_het) == 0)
    ct_order_het <- het_summary_ct[order(desc(pct_high_het))]$cell_type |> unique()
  het_summary_ct[, cell_type := factor(cell_type,
                   levels=unique(c(ct_order_het, cell_types)))]

  p_het_bar <- ggplot(het_summary_ct,
                      aes(x=cell_type, y=pct_high_het, fill=run_label)) +
    geom_col(position=position_dodge(0.75), width=0.7) +
    facet_wrap(~threshold, ncol=1) +
    scale_fill_manual(values=c("#2171b5","#ef6548")) +
    labs(title="% Hits with High Heterogeneity (I²>75%) per Cell Type",
         subtitle="High I² at a hit = effect inconsistent across cohorts",
         x=NULL, y="% hits with I²>75%", fill=NULL) +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=40, hjust=1),
          legend.position="top", panel.grid.major.x=element_blank())
  ggsave(file.path(opt$outdir, "het_pct_high_per_celltype.png"),
         p_het_bar, width=13, height=9, dpi=300)
  cat("  Saved: het_pct_high_per_celltype.png\n")

  write.table(het_summary_ct, file.path(opt$outdir, "het_summary_at_top_hits.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("  Saved: het_summary_at_top_hits.tsv\n")
}

# ── 6. Summary delta table ────────────────────────────────────────────────────
if (!is.null(summ_d) && !is.null(summ_n)) {
  delta <- merge(summ_d, summ_n, by="cell_type", suffixes=c("_design","_notech"))
  delta$delta_lambda  <- delta$lambda_design  - delta$lambda_notech
  delta$delta_n_gw    <- delta$n_gw_design    - delta$n_gw_notech
  delta$delta_n_sugg  <- delta$n_sugg_design  - delta$n_sugg_notech
  write.table(delta, file.path(opt$outdir, "summary_delta.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saved: summary_delta.tsv\n")

  cat("\n=== Quick summary: design vs no-tech ===\n")
  cat(sprintf("  %-20s  λ_design  λ_notech  Δλ      GW_design  GW_notech  ΔGW\n", "cell_type"))
  for (i in seq_len(nrow(delta))) {
    cat(sprintf("  %-20s  %6.3f    %6.3f    %+.3f   %4d       %4d       %+d\n",
      delta$cell_type[i],
      delta$lambda_design[i], delta$lambda_notech[i], delta$delta_lambda[i],
      delta$n_gw_design[i],   delta$n_gw_notech[i],  delta$delta_n_gw[i]))
  }
}

cat("\n========================================\n")
cat("All outputs written to:", opt$outdir, "\n")
cat("  lambda_comparison.png\n")
cat("  n_hits_comparison.png\n")
cat("  ldsc_h2_comparison.png\n")
cat("  ldsc_intercept_comparison.png\n")
cat("  ldsc_comparison.tsv\n")
cat("  top_hits_design_with_het.tsv    — GW+suggestive hits (design run) + HetISq/HetPVal\n")
cat("  top_hits_notech_with_het.tsv    — GW+suggestive hits (no-tech run) + HetISq/HetPVal\n")
cat("  het_at_top_hits_scatter.png     — I² vs -log10(P) scatter coloured by run\n")
cat("  het_pct_high_per_celltype.png   — % hits with I²>75% per cell type\n")
cat("  het_summary_at_top_hits.tsv\n")
cat("  summary_delta.tsv\n")
cat("========================================\n")
