#!/usr/bin/env Rscript
# Compare a newly written cell_proportions.csv / phenotypes_RINT.txt to the SCC copy.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--new_prop",  type = "character"),
  make_option("--old_prop",  type = "character"),
  make_option("--new_pheno", type = "character", default = NULL),
  make_option("--old_pheno", type = "character", default = NULL),
  make_option("--cohort",    type = "character"),
  make_option("--out",       type = "character")
)))

id_col <- function(dt) {
  hits <- intersect(c("specimenID", "sample_id", "IID", "individualID", "V1"), names(dt))
  if (!length(hits)) names(dt)[1] else hits[1]
}

compare_tables <- function(new_path, old_path, label) {
  if (is.null(new_path) || is.null(old_path) || !file.exists(new_path) || !file.exists(old_path)) {
    return(data.table(table = label, status = "MISSING", n_new = NA, n_old = NA,
                      n_shared = NA, median_r = NA, min_r = NA, mean_abs_diff = NA))
  }
  a <- fread(new_path)
  b <- fread(old_path)
  ia <- id_col(a); ib <- id_col(b)
  setnames(a, ia, "id"); setnames(b, ib, "id")
  a[, id := as.character(id)]; b[, id := as.character(id)]
  nums <- intersect(names(a)[vapply(a, is.numeric, TRUE)],
                    names(b)[vapply(b, is.numeric, TRUE)])
  nums <- setdiff(nums, "id")
  shared <- intersect(a$id, b$id)
  if (!length(shared) || !length(nums)) {
    return(data.table(table = label, status = "NO_OVERLAP", n_new = nrow(a), n_old = nrow(b),
                      n_shared = length(shared), median_r = NA, min_r = NA, mean_abs_diff = NA))
  }
  a2 <- a[id %in% shared]; b2 <- b[id %in% shared]
  setkey(a2, id); setkey(b2, id)
  rs <- sapply(nums, function(col) suppressWarnings(cor(a2[[col]], b2[[col]], use = "complete.obs")))
  diffs <- sapply(nums, function(col) mean(abs(a2[[col]] - b2[[col]]), na.rm = TRUE))
  data.table(table = label, status = "OK", n_new = nrow(a), n_old = nrow(b),
             n_shared = length(shared), median_r = median(rs, na.rm = TRUE),
             min_r = min(rs, na.rm = TRUE), mean_abs_diff = mean(diffs, na.rm = TRUE))
}

res <- rbind(
  compare_tables(opt$new_prop,  opt$old_prop,  "cell_proportions"),
  compare_tables(opt$new_pheno, opt$old_pheno, "phenotypes_RINT")
)
res[, cohort := opt$cohort]
setcolorder(res, "cohort")
fwrite(res, opt$out, sep = "\t")
print(res)
