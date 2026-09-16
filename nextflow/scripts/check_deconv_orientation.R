#!/usr/bin/env Rscript
# Orientation gate for MGP cell-type proportions.
#
# MGP estimates each cell type's relative proportion as PC1 of its marker genes.
# PC1's sign is arbitrary and MGP orients it from the marker loadings; for cell
# types with few markers (PAX6: 7) that orientation can flip between runs, which
# would invert every GWAS BETA for that cell type in that cohort.  This checks,
# per cell type, that the estimated proportion correlates POSITIVELY with the mean
# (per-gene z-scored) expression of that cell type's markers in the same
# batch-corrected matrix MGP consumed.
#
# Exit 0 if every cell type is positively oriented, 2 if any is negative.
# Cell types whose markers are absent from the matrix are reported as NA (not a failure).

suppressPackageStartupMessages({ library(data.table); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cohort",      type = "character"),
  make_option("--proportions", type = "character", help = "cell_proportions.csv (sample id in column 1)"),
  make_option("--corrected",   type = "character", help = "corrected_data.RData (saveRDS matrix genes x samples)"),
  make_option("--marker_file", type = "character", help = "Sonny marker CSV (Ensembl gene ID, Subclass)"),
  make_option("--out",         type = "character"),
  make_option("--weak_r",      type = "double", default = 0.3, help = "flag |r| below this as WEAK")
)))
stopifnot(!is.null(opt$proportions), !is.null(opt$corrected), !is.null(opt$marker_file), !is.null(opt$out))

prop <- fread(opt$proportions); setnames(prop, 1, "id"); prop[, id := as.character(id)]
expr <- readRDS(opt$corrected)
if (is.data.frame(expr)) expr <- as.matrix(expr)
rownames(expr) <- sub("\\.[0-9]+$", "", rownames(expr))
markers <- fread(opt$marker_file, check.names = FALSE)
gcol <- grep("^Ensembl", names(markers), value = TRUE)[1]
ccol <- grep("^Subclass$", names(markers), value = TRUE)[1]
stopifnot(!is.na(gcol), !is.na(ccol))
markers <- markers[!is.na(get(gcol)) & get(gcol) != ""]
markers[, ct := make.names(get(ccol))]

common <- intersect(prop$id, colnames(expr))
if (length(common) < 10) stop("fewer than 10 samples shared between proportions and expression matrix")
prop <- prop[match(common, id)]
expr <- expr[, common, drop = FALSE]

# Plain data.frame for the marker lookup: inside a data.table `[`, a bare `ct`
# would resolve to the column, not the loop variable, and select every marker.
mk <- as.data.frame(markers)
rows <- lapply(setdiff(names(prop), "id"), function(ct) {
  genes <- intersect(mk[mk$ct == ct, gcol], rownames(expr))
  if (length(genes) == 0)
    return(data.table(cohort = opt$cohort, cell_type = ct, n_markers = 0L, r = NA_real_, verdict = "NA_NO_MARKERS"))
  z <- expr[genes, , drop = FALSE]
  z <- sweep(z, 1, rowMeans(z)); sds <- apply(z, 1, sd); sds[sds == 0] <- 1; z <- sweep(z, 1, sds, "/")
  score <- colMeans(z)
  r <- suppressWarnings(cor(score, prop[[ct]]))
  verdict <- if (is.na(r)) "NA" else if (r < 0) "NEGATIVE" else if (r < opt$weak_r) "WEAK" else "OK"
  data.table(cohort = opt$cohort, cell_type = ct, n_markers = length(genes), r = round(r, 4), verdict = verdict)
})
res <- rbindlist(rows)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(res, opt$out, sep = "\t")
print(res[order(r)])

neg <- res[verdict == "NEGATIVE", cell_type]
if (length(neg)) {
  cat(sprintf("\nORIENTATION FAIL [%s]: negative marker correlation for: %s\n", opt$cohort, paste(neg, collapse = ", ")))
  quit(status = 2)
}
cat(sprintf("\nORIENTATION OK [%s]: %d cell types positively oriented (%d weak)\n",
            opt$cohort, sum(res$verdict %in% c("OK", "WEAK")), sum(res$verdict == "WEAK")))
