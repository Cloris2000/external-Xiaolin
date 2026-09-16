#!/usr/bin/env Rscript
# Recover April 7 / June 9 RNA library (cell_proportions specimenID) per GWAS FID.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(RNOmni)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--scc_results", type = "character",
              default = "/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow/results"),
  make_option("--out_dir", type = "character",
              default = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/validation/june9_fid_specimen_maps"),
  make_option("--cohorts", type = "character", default = "ROSMAP,ROSMAP_array,MSBB")
)))

biospec_path <- list(
  ROSMAP = "/project/rrg-shreejoy/ROSMAP/Metadata/ROSMAP_biospecimen_metadata.csv",
  ROSMAP_array = "/project/rrg-shreejoy/ROSMAP/Metadata/ROSMAP_biospecimen_metadata.csv",
  MSBB = "/project/rrg-shreejoy/MSBB/Metadata/MSBB_biospecimen_metadata.csv"
)

rint_from_raw <- function(combined_df, cell_types) {
  bind_rows(lapply(cell_types, function(cell_type) {
    lm_formula <- paste0(
      "scale(", cell_type, ") ~ msex + scale(age_death) + scale(age_death_sex) + ",
      "scale(age_death2) + scale(age_death2_sex)"
    )
    results <- residuals(lm(lm_formula, data = combined_df))
    keep <- !is.na(results)
    out <- as.data.frame(RNOmni::RankNorm(results[keep]))
    names(out)[1] <- "transformed_residuals"
    out$cell_type <- cell_type
    out$FID <- combined_df$FID[keep]
    out
  })) %>% pivot_wider(names_from = "cell_type", values_from = "transformed_residuals")
}

score_vs_scc <- function(new_wide, scc_pheno, cell_types) {
  m <- merge(new_wide, scc_pheno, by = "FID", suffixes = c("_n", "_o"))
  mean(sapply(cell_types, function(ct) {
    suppressWarnings(cor(m[[paste0(ct, "_n")]], m[[paste0(ct, "_o")]], use = "complete.obs"))
  }), na.rm = TRUE)
}

recover_one <- function(cohort) {
  scc <- opt$scc_results
  prop <- read.csv(file.path(scc, cohort, "cell_proportions.csv"), stringsAsFactors = FALSE)
  meta <- read.csv(file.path(scc, cohort, "combined_metrics.csv"), stringsAsFactors = FALSE)
  clin <- read.table(file.path(scc, cohort, "clinical_covariates.txt"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  pheno <- read.table(file.path(scc, cohort, "phenotypes_RINT.txt"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  id_col <- if ("specimenID" %in% names(prop)) "specimenID" else names(prop)[1]
  syn_col <- if ("synapseID" %in% names(meta)) "synapseID" else id_col
  cell_types <- intersect(setdiff(names(prop), c(id_col, "FID", "IID", "individualID")),
                          setdiff(names(pheno), c("FID", "IID")))
  prop$specimenID <- as.character(prop[[id_col]])
  meta$._syn <- as.character(meta[[syn_col]])
  meta$._ind <- as.character(meta$individualID)
  prop_ind <- merge(prop, unique(meta[, c("._syn", "._ind")]),
                    by.x = "specimenID", by.y = "._syn", all.x = TRUE)

  clin$FID <- as.character(clin$FID)
  fid_to_ind <- setNames(clin$FID, clin$FID)
  bp <- biospec_path[[cohort]]
  if (!is.null(bp) && file.exists(bp)) {
    bio <- read.csv(bp, stringsAsFactors = FALSE)
    if ("assay" %in% names(bio) && cohort %in% c("ROSMAP", "ROSMAP_array")) {
      bio <- bio[bio$assay == "wholeGenomeSeq", , drop = FALSE]
    }
    spec_col <- if ("specimenID" %in% names(bio)) "specimenID" else names(bio)[1]
    fid_to_ind[as.character(bio[[spec_col]])] <- as.character(bio$individualID)
  }

  rows <- list()
  for (i in seq_len(nrow(clin))) {
    fid <- clin$FID[i]
    ind <- fid_to_ind[[fid]]
    if (is.null(ind) || is.na(ind)) ind <- fid
    hits <- prop_ind[prop_ind$._ind == ind | prop_ind$specimenID == fid, , drop = FALSE]
    if (nrow(hits) == 0) next
    hits$FID <- fid
    hits$msex <- clin$msex[i]
    hits$age_death <- clin$age_death[i]
    hits$age_death_sex <- clin$age_death_sex[i]
    hits$age_death2 <- clin$age_death2[i]
    hits$age_death2_sex <- clin$age_death2_sex[i]
    rows[[length(rows) + 1]] <- hits
  }
  merged <- bind_rows(rows)
  if (nrow(merged) == 0) {
    cat(cohort, ": no FID-specimen join, skip\n")
    return(invisible(NULL))
  }
  n_per <- merged %>% count(FID)
  multi <- n_per$FID[n_per$n > 1]
  cat(cohort, ": joined", nrow(merged), "rows /", length(unique(merged$FID)),
      "FIDs /", length(multi), "multi-library\n")

  pick <- merged %>%
    group_by(FID) %>%
    arrange(specimenID, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()

  for (fid in multi) {
    cands <- unique(merged$specimenID[merged$FID == fid])
    best_score <- -Inf
    best_spec <- cands[1]
    for (spec in cands) {
      trial <- bind_rows(pick[pick$FID != fid, ],
                         merged[merged$FID == fid & merged$specimenID == spec, ][1, ])
      wide <- tryCatch(rint_from_raw(trial, cell_types), error = function(e) NULL)
      if (is.null(wide)) next
      sc <- score_vs_scc(wide, pheno, cell_types)
      if (is.finite(sc) && sc > best_score) {
        best_score <- sc
        best_spec <- spec
      }
    }
    pick <- bind_rows(pick[pick$FID != fid, ],
                      merged[merged$FID == fid & merged$specimenID == best_spec, ][1, ])
    cat("  ", fid, "->", best_spec, "r", round(best_score, 6), "\n")
  }

  dest <- file.path(opt$out_dir, paste0(cohort, ".tsv"))
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(pick[, c("FID", "specimenID")], dest, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat(cohort, ": wrote", dest, "n=", nrow(pick), "\n")
}

for (co in trimws(strsplit(opt$cohorts, ",", fixed = TRUE)[[1]])) recover_one(co)
