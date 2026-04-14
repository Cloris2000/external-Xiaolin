#!/usr/bin/env Rscript
# compare_mgp_snrnaseq.R
#
# Compares three MGP estimates for ROSMAP DLPFC:
#   1. Pipeline MGP with design= (biological signal protected)
#   2. Pipeline MGP without design= (re-computed here from no-design corrected expression)
#   3. snRNA-seq cell type proportions (ground truth reference)
#
# Outputs (written to results/ROSMAP/comparison_design_matrix/):
#   scatter_<celltype>.png     — per-cell-type scatter: MGP_design vs MGP_nodesign vs snRNAseq
#   correlation_summary.png   — heatmap of Pearson r between each MGP vs snRNAseq per cell type
#   correlation_summary.tsv   — tabular version

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(optparse)
})

option_list <- list(
  make_option("--results_dir", type="character",
              default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP"),
  make_option("--sn_proportions", type="character",
              default="/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv"),
  make_option("--meta_file", type="character",
              default="/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv"),
  make_option("--marker_file", type="character",
              default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/new_MTGnCgG_lfct2.5_Publication.csv"),
  make_option("--batch_covariates", type="character", default="sequencingBatch,libraryPrep")
)
opt <- parse_args(OptionParser(option_list=option_list))

outdir <- file.path(opt$results_dir, "comparison_design_matrix")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ── Cell type name mapping: pipeline MGP → snRNA-seq ─────────────────────────
# Only map cell types that have a clear 1-to-1 equivalent.
# snRNA-seq has finer subtypes (Lamp5.Lhx6 vs Lamp5, L5.IT vs IT, etc.)
# We use the broadest reasonable match.
ct_map <- c(
  "Astrocyte"    = "Astro",
  "Endothelial"  = "Endo",
  "PVALB"        = "Pvalb",
  "VIP"          = "Vip",
  "LAMP5"        = "Lamp5",
  "SST"          = "Sst",
  "PAX6"         = "Pax6",
  "L5.6.NP"      = "L5.6.NP",
  "L6b"          = "L6b",
  "L6.CT"        = "L6.CT",
  "L5.ET"        = "L5.ET",
  "L4.IT"        = "L4.IT",
  "Microglia"    = "Micro",
  "OPC"          = "OPC",
  "Oligodendrocyte" = "Oligo"
  # IT, L5.6.IT.Car3, Pericyte, VLMC: no direct 1-to-1 in snRNA-seq
)

# ── 1. Load pipeline MGP (with design=) ──────────────────────────────────────
cat("Loading pipeline MGP estimates (with design=)...\n")
mgp_design <- read_csv(file.path(opt$results_dir, "cell_proportions.csv"),
                       show_col_types=FALSE)
# specimenID is sample identifier
cat("  Samples:", nrow(mgp_design), "| Cell types:", ncol(mgp_design)-1, "\n")

# ── 2. Compute MGP without design= ───────────────────────────────────────────
cat("Loading expression data for no-design MGP computation...\n")
zscore_mat    <- readRDS(file.path(opt$results_dir, "zscore_data.RData"))
meta_raw      <- read_csv(file.path(opt$results_dir, "metadata_cleaned.csv"),
                           show_col_types=FALSE)

# Find column that aligns with expression sample IDs
best_col <- NULL; best_overlap <- 0
for (i in seq_along(meta_raw)) {
  ov <- length(intersect(colnames(zscore_mat), as.character(meta_raw[[i]])))
  if (ov > best_overlap) { best_overlap <- ov; best_col <- i }
}
cat("  Using metadata column '", colnames(meta_raw)[best_col],
    "' as sample ID (", best_overlap, " samples)\n", sep="")
id_col <- as.character(meta_raw[[best_col]])
meta_raw <- meta_raw[!duplicated(id_col), ]
meta_df  <- as.data.frame(meta_raw)           # convert to data.frame BEFORE setting rownames
rownames(meta_df) <- as.character(meta_df[[best_col]])

# Load top tech covariates
tech_covs_raw <- readLines(file.path(opt$results_dir, "top_tech_covariates.txt"), warn=FALSE)
tech_covs <- gsub('^"|"$', '', trimws(tech_covs_raw[tech_covs_raw != "" & tech_covs_raw != "\"\""]))
tech_covs <- tech_covs[!startsWith(tech_covs, "#")]

# Align samples
common_samples <- intersect(colnames(zscore_mat), rownames(meta_df))
zscore_df      <- as.data.frame(zscore_mat[, common_samples, drop=FALSE])
meta_aligned   <- meta_df[common_samples, , drop=FALSE]

# ── Replicate the pipeline's exact removeBatchEffect() logic ─────────────────
# Pipeline uses: batch1=sequencingBatch, batch2=libraryPrep,
#                covariates=metadata_nobatch (top_tech_cov columns, numeric),
#                design=design_for_batch (or matrix(1,n,1) if no bio vars)
# The no-design version replaces design_for_batch with matrix(1, n, 1) — intercept only.
# Everything else (NA handling, batch vectors, covariate matrix) is identical.

batch_cols    <- trimws(strsplit(opt$batch_covariates, ",")[[1]])
batch_present <- batch_cols[batch_cols %in% colnames(meta_aligned)]

# Build tech covariate matrix — pipeline uses top_tech_cov columns only
tech_present <- intersect(tech_covs, colnames(meta_aligned))
cat("  Tech covariates (", length(tech_present), "):", paste(head(tech_present, 5), collapse=", "),
    if (length(tech_present) > 5) "..." else "", "\n")
metadata_nobatch <- meta_aligned[, tech_present, drop=FALSE]

# Coerce non-numeric columns to numeric (same as pipeline)
for (col in colnames(metadata_nobatch)) {
  if (!is.numeric(metadata_nobatch[[col]])) {
    metadata_nobatch[[col]] <- as.numeric(as.factor(metadata_nobatch[[col]]))
  }
}

# NA handling — mirror pipeline exactly:
# columns with >2 NAs: drop the column
# columns with <=2 NAs: remove those specific samples
if (any(is.na(metadata_nobatch))) {
  na_per_col <- sapply(metadata_nobatch, function(col) sum(is.na(col)))
  cols_too_many_na <- names(na_per_col[na_per_col > 2])
  if (length(cols_too_many_na) > 0) {
    cat("  Dropping tech covariate(s) with >2 NA samples:", paste(cols_too_many_na, collapse=", "), "\n")
    metadata_nobatch <- metadata_nobatch[, !colnames(metadata_nobatch) %in% cols_too_many_na, drop=FALSE]
  }
  na_samples <- rownames(metadata_nobatch)[!complete.cases(metadata_nobatch)]
  if (length(na_samples) > 0) {
    cat("  Removing", length(na_samples), "sample(s) with NA in tech covariates:", paste(na_samples, collapse=", "), "\n")
    keep_samples     <- setdiff(colnames(zscore_df), na_samples)
    zscore_df        <- zscore_df[, keep_samples, drop=FALSE]
    metadata_nobatch <- metadata_nobatch[keep_samples, , drop=FALSE]
    meta_aligned     <- meta_aligned[keep_samples, , drop=FALSE]
  }
}

# Batch vectors (same samples as zscore_df after NA removal)
batch1 <- if (length(batch_present) >= 1) as.vector(meta_aligned[[batch_present[1]]]) else NULL
batch2 <- if (length(batch_present) >= 2) as.vector(meta_aligned[[batch_present[2]]]) else NULL

has_covariates <- ncol(metadata_nobatch) > 0
has_batch1     <- !is.null(batch1)
has_batch2     <- !is.null(batch2)
n_nd           <- ncol(zscore_df)

# design= intercept-only (no biological variable protection)
design_intercept <- matrix(1, n_nd, 1)

cat("  Applying removeBatchEffect WITHOUT design= (",
    if (has_batch1) batch_present[1] else "", "+",
    if (has_batch2) batch_present[2] else "", "+",
    ncol(metadata_nobatch), "tech covs,", n_nd, "samples)...\n")

nodesign_mat <- tryCatch({
  if (has_batch1 && has_batch2 && has_covariates) {
    limma::removeBatchEffect(zscore_df, batch=batch1, batch2=batch2,
                             covariates=data.matrix(metadata_nobatch),
                             design=design_intercept)
  } else if (has_batch1 && has_covariates) {
    limma::removeBatchEffect(zscore_df, batch=batch1,
                             covariates=data.matrix(metadata_nobatch),
                             design=design_intercept)
  } else if (has_batch1 && has_batch2) {
    limma::removeBatchEffect(zscore_df, batch=batch1, batch2=batch2,
                             design=design_intercept)
  } else if (has_batch1) {
    limma::removeBatchEffect(zscore_df, batch=batch1, design=design_intercept)
  } else if (has_covariates) {
    limma::removeBatchEffect(zscore_df, covariates=data.matrix(metadata_nobatch),
                             design=design_intercept)
  } else {
    cat("  NOTE: no batch or tech covariates — no-design result = uncorrected\n")
    as.matrix(zscore_df)
  }
}, error = function(e) { stop("removeBatchEffect (no-design) failed: ", e$message) })

# Run MGP on no-design corrected matrix
cat("  Running MGP on no-design corrected matrix...\n")
markers <- read_csv(opt$marker_file, show_col_types=FALSE)
cat("  Marker file loaded:", nrow(markers), "markers,", ncol(markers)-1, "cell types\n")

# MGP estimation using marker gene means (z-score of mean z-score across markers)
run_mgp <- function(expr_mat, marker_df) {
  cell_types <- colnames(marker_df)[-1]  # first col is gene name
  gene_col   <- colnames(marker_df)[1]
  results <- lapply(cell_types, function(ct) {
    markers_ct <- marker_df[[gene_col]][marker_df[[ct]] == 1]
    genes_present <- intersect(markers_ct, rownames(expr_mat))
    if (length(genes_present) < 3) return(rep(NA, ncol(expr_mat)))
    colMeans(expr_mat[genes_present, , drop=FALSE], na.rm=TRUE)
  })
  names(results) <- cell_types
  as.data.frame(do.call(cbind, results))
}

# Check marker file structure
cat("  Marker file columns (first 5):", paste(colnames(markers)[1:min(5,ncol(markers))], collapse=", "), "\n")
cat("  Marker gene column:", colnames(markers)[1], "\n")

mgp_nodesign_raw <- run_mgp(nodesign_mat, as.data.frame(markers))
mgp_nodesign_raw$specimenID <- colnames(nodesign_mat)

# Rename columns to match pipeline naming
colnames(mgp_nodesign_raw) <- gsub(" ", ".", gsub("/", ".", colnames(mgp_nodesign_raw)))
cat("  No-design MGP computed:", nrow(mgp_nodesign_raw), "samples,",
    ncol(mgp_nodesign_raw)-1, "cell types\n")

# ── 3. Load and harmonize snRNA-seq proportions ───────────────────────────────
cat("Loading snRNA-seq proportions...\n")
sn_raw <- read_csv(opt$sn_proportions, show_col_types=FALSE)

# Your harmonization code
meta_rnaseq <- read_csv(opt$meta_file, show_col_types=FALSE)
matcher <- meta_rnaseq[!duplicated(meta_rnaseq$specimenID), c("specimenID","projid","Study","msex","age_death","tissue")]
matcher <- matcher[!duplicated(matcher$projid), ]
sn_props <- merge(sn_raw, matcher, by.x="ID", by.y="projid", all.x=TRUE, all.y=FALSE)
sn_props <- sn_props[!is.na(sn_props$specimenID), ]
cat("  snRNA-seq samples with matched specimenID:", nrow(sn_props), "\n")
# Rename snRNA-seq columns to match ct_map values
colnames(sn_props) <- gsub(" ", ".", gsub("/", ".", colnames(sn_props)))

# ── 4. Merge all three ────────────────────────────────────────────────────────
cat("Merging datasets...\n")
# pipeline MGP already uses specimenID
mgp_d  <- as.data.frame(mgp_design)
mgp_nd <- mgp_nodesign_raw

# Merge design vs no-design
merged <- merge(mgp_d,  mgp_nd,  by="specimenID", suffixes=c("_design","_nodesign"))
# Merge with snRNA-seq (by specimenID)
merged <- merge(merged, sn_props, by="specimenID")
cat("  Final merged samples:", nrow(merged), "\n")

# ── 5. Correlation analysis ───────────────────────────────────────────────────
cat("Computing correlations with snRNA-seq...\n")
cor_results <- lapply(names(ct_map), function(mgp_ct) {
  sn_ct <- ct_map[mgp_ct]

  col_design   <- mgp_ct  # already named by pipeline
  col_nodesign <- paste0(mgp_ct, "_nodesign")

  # Fallback: some cell types may have dot-substituted names
  col_design_alt   <- gsub("\\.", "\\.", mgp_ct)

  if (!col_design %in% colnames(merged) || !col_nodesign %in% colnames(merged) || !sn_ct %in% colnames(merged)) {
    cat("  SKIP", mgp_ct, "— column(s) missing\n")
    return(NULL)
  }
  x_d  <- merged[[col_design]]
  x_nd <- merged[[col_nodesign]]
  y    <- merged[[sn_ct]]

  ok <- complete.cases(x_d, x_nd, y)
  if (sum(ok) < 10) { cat("  SKIP", mgp_ct, "— too few complete cases\n"); return(NULL) }

  r_design   <- cor(x_d[ok],  y[ok], method="pearson")
  r_nodesign <- cor(x_nd[ok], y[ok], method="pearson")
  cat(sprintf("  %-20s  r_design=%6.3f  r_nodesign=%6.3f  delta=%+.3f\n",
              mgp_ct, r_design, r_nodesign, r_design - r_nodesign))
  data.frame(cell_type=mgp_ct, sn_cell_type=sn_ct,
             n=sum(ok),
             r_design=r_design, r_nodesign=r_nodesign,
             delta_r=r_design - r_nodesign)
})
cor_df <- do.call(rbind, Filter(Negate(is.null), cor_results))
write.table(cor_df, file.path(outdir, "correlation_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("  Saved: correlation_summary.tsv\n")

# ── 6. Summary bar plot: r_design vs r_nodesign per cell type ─────────────────
cat("Plotting correlation summary...\n")
cor_long <- cor_df %>%
  select(cell_type, r_design, r_nodesign) %>%
  pivot_longer(cols=c(r_design, r_nodesign),
               names_to="method", values_to="r") %>%
  mutate(method=recode(method,
    r_design   = "With design= (bio-signal protected)",
    r_nodesign = "Without design="))

# Order cell types by r_design descending
ct_order <- cor_df %>% arrange(desc(r_design)) %>% pull(cell_type)
cor_long$cell_type <- factor(cor_long$cell_type, levels=ct_order)

p_bar <- ggplot(cor_long, aes(x=cell_type, y=r, fill=method)) +
  geom_col(position=position_dodge(0.75), width=0.7) +
  geom_hline(yintercept=0, linewidth=0.4) +
  scale_fill_manual(values=c(
    "With design= (bio-signal protected)" = "#2171b5",
    "Without design="                     = "#ef6548")) +
  labs(title="ROSMAP: MGP vs snRNA-seq proportions — with vs without design= in removeBatchEffect()",
       subtitle=paste0("n=", unique(cor_df$n)[1], " samples with matched snRNA-seq"),
       x=NULL, y="Pearson r (MGP vs snRNA-seq)",
       fill=NULL) +
  theme_bw(base_size=13) +
  theme(axis.text.x=element_text(angle=40, hjust=1),
        legend.position="top",
        panel.grid.major.x=element_blank())
ggsave(file.path(outdir, "correlation_summary_barplot.png"),
       p_bar, width=12, height=6, dpi=300)
cat("  Saved: correlation_summary_barplot.png\n")

# Delta plot: positive = design= improved correlation with snRNA-seq
p_delta <- ggplot(cor_df, aes(x=reorder(cell_type, delta_r), y=delta_r,
                               fill=delta_r > 0)) +
  geom_col(width=0.7) +
  geom_hline(yintercept=0, linewidth=0.4) +
  scale_fill_manual(values=c("TRUE"="#2171b5", "FALSE"="#ef6548"),
                    labels=c("TRUE"="design= improves r", "FALSE"="design= reduces r")) +
  labs(title="Delta r = r(design=) − r(no design=)",
       subtitle="Positive = design= makes MGP more similar to snRNA-seq",
       x=NULL, y="Δ Pearson r", fill=NULL) +
  theme_bw(base_size=13) +
  theme(axis.text.x=element_text(angle=40, hjust=1),
        legend.position="top",
        panel.grid.major.x=element_blank())
ggsave(file.path(outdir, "correlation_delta_barplot.png"),
       p_delta, width=10, height=5, dpi=300)
cat("  Saved: correlation_delta_barplot.png\n")

# ── 7. Per-cell-type scatter plots ────────────────────────────────────────────
cat("Plotting per-cell-type scatters...\n")
for (mgp_ct in names(ct_map)) {
  sn_ct        <- ct_map[mgp_ct]
  col_design   <- mgp_ct
  col_nodesign <- paste0(mgp_ct, "_nodesign")
  if (!all(c(col_design, col_nodesign, sn_ct) %in% colnames(merged))) next

  plot_df <- merged[, c(col_design, col_nodesign, sn_ct)] %>%
    setNames(c("design", "nodesign", "snRNAseq")) %>%
    filter(complete.cases(.)) %>%
    pivot_longer(cols=c(design, nodesign), names_to="method", values_to="MGP") %>%
    mutate(method=recode(method,
      design   = "With design=",
      nodesign = "Without design="))

  r_vals <- plot_df %>%
    group_by(method) %>%
    summarise(r=round(cor(MGP, snRNAseq, method="pearson"), 3), .groups="drop") %>%
    mutate(label=paste0(method, "  r=", r))

  p_sc <- ggplot(plot_df, aes(x=MGP, y=snRNAseq, colour=method)) +
    geom_point(alpha=0.35, size=0.8) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.9) +
    scale_colour_manual(values=c("With design="="#2171b5","Without design="="#ef6548")) +
    annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.4, size=3.5,
             label=paste(r_vals$label, collapse="\n")) +
    labs(title=paste0(mgp_ct, " — MGP estimate vs snRNA-seq proportion"),
         x="MGP estimate (z-score)", y="snRNA-seq proportion",
         colour=NULL) +
    theme_bw(base_size=12) +
    theme(legend.position="top")
  ggsave(file.path(outdir, paste0("scatter_", gsub("\\.", "_", mgp_ct), ".png")),
         p_sc, width=6, height=5, dpi=250)
}
cat("  Scatter plots saved to:", outdir, "\n")

cat("\n========================================\n")
cat("All outputs written to:", outdir, "\n")
cat("  correlation_summary.tsv\n")
cat("  correlation_summary_barplot.png\n")
cat("  correlation_delta_barplot.png\n")
cat("  scatter_<celltype>.png (", length(names(ct_map)), "cell types)\n")
cat("========================================\n")
