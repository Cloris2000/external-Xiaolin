#!/usr/bin/env Rscript
# Compare a new REGENIE .raw_p to the SCC copy (BETA/SE/P correlation + lead SNP).

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--new", type = "character"),
  make_option("--old", type = "character"),
  make_option("--cohort", type = "character"),
  make_option("--cell_type", type = "character"),
  make_option("--out", type = "character")
)))

read_one <- function(path) {
  dt <- fread(path)
  nms <- names(dt)
  id <- intersect(c("ID", "SNP", "MarkerName"), nms)[1]
  beta <- intersect(c("BETA", "Effect"), nms)[1]
  se <- intersect(c("SE", "StdErr"), nms)[1]
  p <- intersect(c("P", "Pvalue", "P-value"), nms)[1]
  setnames(dt, c(id, beta, se, p), c("id", "beta", "se", "p"), skip_absent = TRUE)
  dt[, `:=`(id = as.character(id), beta = as.numeric(beta), se = as.numeric(se), p = as.numeric(p))]
  dt[!is.na(id) & !is.na(p)]
}

a <- read_one(opt$new)
b <- read_one(opt$old)
setkey(a, id); setkey(b, id)
m <- a[b, nomatch = 0]
lead_a <- a[which.min(p)]$id
lead_b <- b[which.min(p)]$id
out <- data.table(
  cohort = opt$cohort, cell_type = opt$cell_type,
  n_new = nrow(a), n_old = nrow(b), n_shared = nrow(m),
  r_beta = if (nrow(m)) cor(m$beta, m$i.beta, use = "complete.obs") else NA_real_,
  r_se   = if (nrow(m)) cor(m$se,   m$i.se,   use = "complete.obs") else NA_real_,
  r_logp = if (nrow(m)) cor(-log10(pmax(m$p, 1e-300)), -log10(pmax(m$i.p, 1e-300)), use = "complete.obs") else NA_real_,
  lead_new = lead_a, lead_old = lead_b, lead_same = identical(lead_a, lead_b)
)
fwrite(out, opt$out, sep = "\t")
print(out)
