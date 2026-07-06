#!/usr/bin/env Rscript
# regen_eigencor_plots.R
#
# Regenerates before/after batch-correction eigencor plots for all cohorts
# using already-saved intermediate files (zscore_data.RData, corrected_data.RData,
# metadata_cleaned.csv, top_tech_covariates.txt).
#
# Outputs per cohort (written into results/<cohort>/):
#   eigencor_plot_before_correction.png        — raw z-score data
#   eigencor_plot_nodesign_correction.png      — corrected WITHOUT design= (no bio-variable protection)
#   eigencor_plot_removedBatchEff_cov.png      — corrected WITH design= (bio-variable protected, overwritten)
#
# All three plots show ALL candidate tech covariates so differences are directly visible.
#
# Usage:
#   Rscript scripts/regen_eigencor_plots.R
#   Rscript scripts/regen_eigencor_plots.R --cohorts ROSMAP,Mayo

suppressPackageStartupMessages({
  library(PCAtools)
  library(limma)
  library(optparse)
})

option_list <- list(
  make_option("--results_dir", type="character",
              default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results",
              help="Base results directory [default: %default]"),
  make_option("--cohorts", type="character", default=NULL,
              help="Comma-separated list of cohorts to process (default: all 13)"),
  make_option("--batch_covariates", type="character",
              default="sequencingBatch,libraryPrep",
              help="Comma-separated batch column names used in removeBatchEffect [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

all_cohorts <- c("ROSMAP", "ROSMAP_array", "Mayo", "MSBB",
                 "CMC_MSSM", "CMC_PENN", "CMC_PITT",
                 "GTEx", "NABEC", "GVEX",
                 "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650")

cohorts <- if (!is.null(opt$cohorts)) {
  strsplit(opt$cohorts, ",")[[1]]
} else {
  all_cohorts
}

make_eigencor <- function(pca_obj, tech_covs, title_str) {
  eigencorplot(pca_obj,
               metavars = tech_covs,
               cexLabY = 0.5,
               rotLabY = 0.8,
               corFUN = "pearson",
               main = title_str,
               titleX = "PCs",
               titleY = "technical covariates",
               corMultipleTestCorrection = "hochberg",
               signifSymbols = c('***', '**', '*', ''),
               signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
               scale = FALSE)
}

save_png <- function(plot_obj, path, width=15, height=10, res=300) {
  png(path, width=width, height=height, units='in', res=res)
  print(plot_obj)
  dev.off()
  cat("    Saved:", path, "\n")
}

for (cohort in cohorts) {
  cat("\n========================================\n")
  cat("Processing:", cohort, "\n")
  cat("========================================\n")

  base <- file.path(opt$results_dir, cohort)

  zscore_path    <- file.path(base, "zscore_data.RData")
  corrected_path <- file.path(base, "corrected_data.RData")
  meta_path      <- file.path(base, "metadata_cleaned.csv")
  techcov_path   <- file.path(base, "top_tech_covariates.txt")

  # Verify files
  missing <- c()
  for (f in c(zscore_path, corrected_path, meta_path, techcov_path)) {
    if (!file.exists(f)) missing <- c(missing, basename(f))
  }
  if (length(missing) > 0) {
    cat("  SKIP — missing files:", paste(missing, collapse=", "), "\n")
    next
  }

  # Load inputs — files are saved with saveRDS() so use readRDS(), not load()
  cat("  Loading zscore_data.RData...\n")
  zscore_mat <- readRDS(zscore_path)

  cat("  Loading corrected_data.RData...\n")
  corrected_mat <- readRDS(corrected_path)

  cat("  Loading metadata_cleaned.csv...\n")
  meta_raw <- read.csv(meta_path, check.names=FALSE)

  # Find which metadata column best aligns with expression matrix sample IDs.
  # The script loads the expression matrix first, so we can check overlap here.
  best_col <- NULL
  best_overlap <- 0
  for (i in seq_along(meta_raw)) {
    ov <- length(intersect(colnames(zscore_mat), as.character(meta_raw[[i]])))
    if (ov > best_overlap) { best_overlap <- ov; best_col <- i }
  }
  if (is.null(best_col) || best_overlap == 0) {
    cat("  SKIP — no metadata column overlaps with expression sample IDs\n")
    next
  }
  cat("  Using metadata column '", colnames(meta_raw)[best_col],
      "' as sample ID (", best_overlap, " overlapping samples)\n", sep="")
  id_col <- as.character(meta_raw[[best_col]])
  if (anyDuplicated(id_col)) {
    cat("  NOTE: duplicate sample IDs — keeping first occurrence of each\n")
    meta_raw <- meta_raw[!duplicated(id_col), ]
    id_col   <- as.character(meta_raw[[best_col]])
  }
  rownames(meta_raw) <- id_col
  meta <- meta_raw[, -best_col, drop=FALSE]

  # Load the *selected* (top-N) tech covariates — used for labelling only
  cat("  Loading top_tech_covariates.txt (selected covariates)...\n")
  selected_covs <- readLines(techcov_path, warn=FALSE)
  selected_covs <- trimws(selected_covs[selected_covs != "" & !startsWith(selected_covs, "#") & selected_covs != "\"\""])
  selected_covs <- gsub('^"|"$', '', selected_covs)
  cat("  Selected covariates (", length(selected_covs), "):", paste(head(selected_covs, 5), collapse=", "),
      if (length(selected_covs) > 5) "..." else "", "\n")

  if (length(selected_covs) == 0) {
    cat("  SKIP — no technical covariates found\n")
    next
  }

  # Align samples between expression and metadata
  common_samples <- intersect(colnames(zscore_mat), rownames(meta))
  zscore_aligned    <- zscore_mat[, common_samples, drop=FALSE]
  corrected_aligned <- corrected_mat[, common_samples, drop=FALSE]
  meta_aligned      <- meta[common_samples, , drop=FALSE]

  # Build the FULL candidate tech-covariate list from all numeric-ish metadata columns.
  # We exclude obvious non-tech columns (IDs, biological variables, clinical prefixes)
  # using the same logic as remove_tech_covar.R, then keep only those present in metadata.
  non_tech_cols <- c(
    "X", "specimenID", "individualID", "sampleID", "id", "synapseID",
    "subjectID", "Subject_ID", "dbGaP_Subject_ID", "projid", "wgs_id",
    "subject_id", "biospecimen_repository_sample_id",
    "msex", "sex", "Sex", "ageDeath", "age_death", "age", "AOD",
    "cogdx", "diagnosis", "Diagnosis", "primaryDiagnosis",
    "PMI", "pmi", "RIN", "rin",
    "tissue", "Tissue", "Region", "BrainRegion",
    "notes", "batch", "Batch"
  )
  non_tech_pats <- c("^MH", "^DTH", "^LB", "^TRISCH", "^TRCL", "^TRORG",
                     "^TRAMP", "^TRCRTMP", "^TRTPT", "^TRVNT", "^TRDNI",
                     "^DTHPRNINT", "^DTHWT", "^DTHCLS", "^DTHTYP", "^DTHCAT",
                     "^DTHICD", "^DTHFUC", "^DTHLUC", "^DTHLU", "^DTHCOD",
                     "^DTHSEASON", "^DTHTIME", "^DTHPLCE", "^executed\\.")

  all_meta_cols <- colnames(meta_aligned)
  pattern_flagged <- sapply(all_meta_cols, function(col) {
    any(sapply(non_tech_pats, function(pat) grepl(pat, col)))
  })
  all_tech_covs <- all_meta_cols[
    !all_meta_cols %in% non_tech_cols & !pattern_flagged
  ]
  # Keep only columns with > 1 unique non-NA value (useless otherwise)
  all_tech_covs <- all_tech_covs[sapply(all_tech_covs, function(col) {
    length(unique(na.omit(meta_aligned[[col]]))) > 1
  })]

  cat("  All candidate tech covariates for plot (", length(all_tech_covs), "):",
      paste(head(all_tech_covs, 5), collapse=", "),
      if (length(all_tech_covs) > 5) "..." else "", "\n")

  if (length(all_tech_covs) == 0) {
    cat("  SKIP — no candidate tech covariates found in metadata\n")
    next
  }

  # Plot height scales with number of covariates
  plot_height <- max(10, ceiling(length(all_tech_covs) * 0.35))

  cat("  Running PCA on uncorrected data (", ncol(zscore_aligned), "samples x",
      nrow(zscore_aligned), "genes)...\n")
  p_before <- tryCatch(
    PCAtools::pca(zscore_aligned, metadata=meta_aligned),
    error = function(e) { cat("  ERROR running PCA (before):", e$message, "\n"); NULL }
  )

  # ── Compute no-design correction (batch + tech covariates, no design= argument) ──
  # This shows what removeBatchEffect does WITHOUT protecting biological signal.
  cat("  Computing removeBatchEffect WITHOUT design= ...\n")
  batch_cov_names <- trimws(strsplit(opt$batch_covariates, ",")[[1]])
  batch1_present  <- batch_cov_names[batch_cov_names %in% colnames(meta_aligned)]
  tech_for_rbe    <- intersect(selected_covs, colnames(meta_aligned))

  nodesign_mat <- tryCatch({
    mat_nd <- zscore_aligned
    # Build covariates matrix from tech covariates (numeric only)
    if (length(tech_for_rbe) > 0) {
      cov_df <- meta_aligned[colnames(mat_nd), tech_for_rbe, drop=FALSE]
      # coerce factors/characters to numeric
      for (col in colnames(cov_df)) {
        if (!is.numeric(cov_df[[col]])) cov_df[[col]] <- as.numeric(as.factor(cov_df[[col]]))
      }
      # drop NA samples
      ok_samples <- rownames(cov_df)[complete.cases(cov_df)]
      mat_nd  <- mat_nd[, ok_samples, drop=FALSE]
      cov_df  <- cov_df[ok_samples, , drop=FALSE]
      cov_mat <- data.matrix(cov_df)
    } else {
      cov_mat <- NULL
    }
    # Apply removeBatchEffect without design=
    if (length(batch1_present) >= 2) {
      b1 <- meta_aligned[colnames(mat_nd), batch1_present[1]]
      b2 <- meta_aligned[colnames(mat_nd), batch1_present[2]]
      if (!is.null(cov_mat)) {
        limma::removeBatchEffect(mat_nd, batch=b1, batch2=b2, covariates=cov_mat)
      } else {
        limma::removeBatchEffect(mat_nd, batch=b1, batch2=b2)
      }
    } else if (length(batch1_present) == 1) {
      b1 <- meta_aligned[colnames(mat_nd), batch1_present[1]]
      if (!is.null(cov_mat)) {
        limma::removeBatchEffect(mat_nd, batch=b1, covariates=cov_mat)
      } else {
        limma::removeBatchEffect(mat_nd, batch=b1)
      }
    } else {
      if (!is.null(cov_mat)) {
        limma::removeBatchEffect(mat_nd, covariates=cov_mat)
      } else {
        cat("  NOTE: no batch or tech covariates to apply — no-design result = uncorrected\n")
        mat_nd
      }
    }
  }, error = function(e) {
    cat("  WARNING: no-design removeBatchEffect failed:", e$message, "\n")
    NULL
  })

  cat("  Running PCA on no-design corrected data...\n")
  p_nodesign <- if (!is.null(nodesign_mat)) {
    # align metadata to possibly reduced sample set
    nd_samples <- intersect(colnames(nodesign_mat), rownames(meta_aligned))
    tryCatch(
      PCAtools::pca(nodesign_mat[, nd_samples, drop=FALSE],
                    metadata=meta_aligned[nd_samples, , drop=FALSE]),
      error = function(e) { cat("  ERROR running PCA (no-design):", e$message, "\n"); NULL }
    )
  } else NULL

  cat("  Running PCA on corrected data (with design=)...\n")
  p_after <- tryCatch(
    PCAtools::pca(corrected_aligned, metadata=meta_aligned),
    error = function(e) { cat("  ERROR running PCA (after):", e$message, "\n"); NULL }
  )

  # ── Before plot — ALL candidate covariates ──
  if (!is.null(p_before)) {
    tryCatch({
      plt <- make_eigencor(p_before, all_tech_covs,
                           paste0(cohort, " — BEFORE correction (all tech covariates; selected: ",
                                  length(intersect(selected_covs, all_tech_covs)), " of ", length(all_tech_covs), ")"))
      save_png(plt, file.path(base, "eigencor_plot_before_correction.png"),
               height=plot_height)
    }, error = function(e) {
      cat("  WARNING: before-correction plot failed:", e$message, "\n")
    })
  }

  # ── No-design plot — ALL candidate covariates ──
  if (!is.null(p_nodesign)) {
    tryCatch({
      plt <- make_eigencor(p_nodesign, all_tech_covs,
                           paste0(cohort, " — AFTER correction WITHOUT design= (batch+tech only; selected: ",
                                  length(intersect(selected_covs, all_tech_covs)), " of ", length(all_tech_covs), ")"))
      save_png(plt, file.path(base, "eigencor_plot_nodesign_correction.png"),
               height=plot_height)
    }, error = function(e) {
      cat("  WARNING: no-design plot failed:", e$message, "\n")
    })
  }

  # ── After (with design=) plot — ALL candidate covariates ──
  if (!is.null(p_after)) {
    tryCatch({
      plt <- make_eigencor(p_after, all_tech_covs,
                           paste0(cohort, " — AFTER correction (all tech covariates; selected: ",
                                  length(intersect(selected_covs, all_tech_covs)), " of ", length(all_tech_covs), ")"))
      save_png(plt, file.path(base, "eigencor_plot_removedBatchEff_cov.png"),
               height=plot_height)
    }, error = function(e) {
      cat("  WARNING: after-correction plot failed:", e$message, "\n")
    })
  }

  cat("  Done:", cohort, "\n")
  gc()
}

cat("\n========================================\n")
cat("All cohorts processed.\n")
cat("========================================\n")
