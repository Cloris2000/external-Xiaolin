#!/usr/bin/env Rscript
# Compare the corrected meta-analysis (v2) with the SCC manuscript meta, per cell type,
# and draw a Miami plot (new above the axis, SCC mirrored below).
#
#   Rscript compare_meta_v2_vs_scc.R --cell_type Astrocyte \
#     --new_dir  /scratch/zhoux156/results/meta_analysis_15cohorts_hg19_v2 \
#     --scc_dir  /project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow/results/meta_analysis_15cohorts \
#     --out_dir  /scratch/zhoux156/results/meta_analysis_15cohorts_hg19_v2/compare_manuscript
#
# SCC side uses METAL's <prefix>1.tbl (the table behind the manuscript figures);
# the plain <prefix>.tbl in that tree is a larger pre-filter artifact.
#
# Outputs, per cell type:
#   <ct>_miami.png          Miami plot, gw/suggestive lines, lead SNP labelled on each side
#   <ct>_leads.tsv          top loci (new and SCC) side by side with the other run's stats
#   <ct>_summary.tsv        one row: counts, correlations, lead concordance, locus overlap
#
# A "locus" is a 1 Mb window around a lead; loci are called independently in each
# run and matched by position, so a lead that merely shifts within a locus counts
# as shared rather than as a discordance.

suppressPackageStartupMessages({
  library(data.table); library(optparse); library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cell_type", type = "character"),
  make_option("--new_dir",   type = "character"),
  make_option("--scc_dir",   type = "character"),
  make_option("--out_dir",   type = "character"),
  make_option("--gw",        type = "double", default = 5e-8),
  make_option("--sugg",      type = "double", default = 1e-5),
  make_option("--window",    type = "double", default = 1e6)
)))
ct <- opt$cell_type
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

pick <- function(dir, ct, scc) {
  f <- list.files(dir, pattern = sprintf("^%s_meta_analysis_.*%s\\.tbl$",
                                         gsub("\\.", "\\\\.", ct), if (scc) "1" else ""),
                  full.names = TRUE)
  if (scc) f <- f[grepl("1\\.tbl$", f)] else f <- f[!grepl("1\\.tbl$", f)]
  f <- f[file.size(f) > 0]
  if (!length(f)) stop("no .tbl for ", ct, " in ", dir)
  f[which.max(file.size(f))]
}

read_meta <- function(path) {
  dt <- fread(path, select = c("MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","P-value","Direction"),
              showProgress = FALSE)
  setnames(dt, c("snp","a1","a2","freq","beta","se","p","dir"))
  dt[, c("chr","pos") := tstrsplit(snp, ":", keep = 1:2)]
  dt[, chr := as.integer(sub("^chr", "", chr))][, pos := as.numeric(pos)][, p := as.numeric(p)]
  dt[!is.na(chr) & !is.na(pos) & !is.na(p) & p > 0 & p <= 1]
}

# Greedy 1 Mb clumping on p-value order
call_loci <- function(dt, thresh, window) {
  s <- dt[p < thresh][order(p)]
  if (!nrow(s)) return(s[0])
  keep <- logical(nrow(s)); taken <- list()
  for (i in seq_len(nrow(s))) {
    c_i <- s$chr[i]; p_i <- s$pos[i]; ok <- TRUE
    for (t in taken) if (t[1] == c_i && abs(t[2] - p_i) < window) { ok <- FALSE; break }
    if (ok) { keep[i] <- TRUE; taken[[length(taken) + 1]] <- c(c_i, p_i) }
  }
  s[keep]
}

new <- read_meta(pick(opt$new_dir, ct, FALSE))
scc <- read_meta(pick(opt$scc_dir, ct, TRUE))
setkey(new, snp); setkey(scc, snp)
m <- merge(new, scc, by = "snp", suffixes = c(".new", ".scc"))

ln <- call_loci(new, opt$sugg, opt$window); ls_ <- call_loci(scc, opt$sugg, opt$window)
matched <- function(a, b) if (!nrow(a) || !nrow(b)) rep(FALSE, nrow(a)) else
  sapply(seq_len(nrow(a)), function(i) any(b$chr == a$chr[i] & abs(b$pos - a$pos[i]) < opt$window))
ln[, in_scc := matched(ln, ls_)]; ls_[, in_new := matched(ls_, ln)]

# per-locus table: each run's lead with the other run's stats at that SNP
side <- function(leads, other, label) {
  if (!nrow(leads)) return(NULL)
  o <- other[match(leads$snp, other$snp)]
  data.table(cell_type = ct, source = label, snp = leads$snp, chr = leads$chr, pos = leads$pos,
             effect_allele = leads$a1, beta = leads$beta, p = leads$p, dir = leads$dir,
             other_beta = o$beta, other_p = o$p,
             shared_locus = if (label == "new") leads$in_scc else leads$in_new)
}
leads <- rbindlist(list(side(ln, scc, "new"), side(ls_, new, "scc")), fill = TRUE)
if (!is.null(leads) && nrow(leads)) fwrite(leads[order(source, p)], file.path(opt$out_dir, paste0(ct, "_leads.tsv")), sep = "\t")

top_new <- if (nrow(new)) new[which.min(p)] else NULL
top_scc <- if (nrow(scc)) scc[which.min(p)] else NULL
summ <- data.table(
  cell_type = ct,
  n_new = nrow(new), n_scc = nrow(scc), n_shared = nrow(m),
  n_new_only = nrow(new) - nrow(m), n_scc_only = nrow(scc) - nrow(m),
  gw_new = new[p < opt$gw, .N], gw_scc = scc[p < opt$gw, .N],
  sugg_new = new[p < opt$sugg, .N], sugg_scc = scc[p < opt$sugg, .N],
  loci_new = nrow(ln), loci_scc = nrow(ls_),
  loci_shared = sum(ln$in_scc), loci_new_only = sum(!ln$in_scc), loci_scc_only = sum(!ls_$in_new),
  r_beta = if (nrow(m)) cor(m$beta.new, m$beta.scc, use = "complete.obs") else NA_real_,
  r_logp = if (nrow(m)) cor(-log10(pmax(m$p.new, 1e-300)), -log10(pmax(m$p.scc, 1e-300)), use = "complete.obs") else NA_real_,
  top_new = top_new$snp, top_new_p = top_new$p,
  top_scc = top_scc$snp, top_scc_p = top_scc$p,
  top_same = identical(top_new$snp, top_scc$snp)
)
fwrite(summ, file.path(opt$out_dir, paste0(ct, "_summary.tsv")), sep = "\t")
print(summ)

# ---- Miami plot: new above, SCC mirrored below -----------------------------
thin <- function(dt, keep_all_below = 1e-4, frac = 0.02) {
  sig <- dt[p < keep_all_below]
  rest <- dt[p >= keep_all_below]
  if (nrow(rest) > 0) rest <- rest[sample(.N, max(1, floor(.N * frac)))]
  rbind(sig, rest)
}
pn <- thin(new)[, .(chr, pos, p, side = "new")]
ps <- thin(scc)[, .(chr, pos, p, side = "scc")]
pl <- rbind(pn, ps)[chr %in% 1:22]
off <- pl[, .(mx = max(pos, na.rm = TRUE)), by = chr][order(chr)]
off[, cum := cumsum(as.numeric(mx)) - mx]
pl <- merge(pl, off[, .(chr, cum)], by = "chr")
pl[, x := pos + cum]
pl[, y := ifelse(side == "new", -log10(p), log10(p))]
axis_df <- pl[, .(center = mean(range(x))), by = chr][order(chr)]
ymax <- max(abs(pl$y), na.rm = TRUE)

p <- ggplot(pl, aes(x, y, colour = interaction(chr %% 2 == 0, side))) +
  geom_point(size = 0.35, alpha = 0.7) +
  geom_hline(yintercept = c(-log10(opt$gw), log10(opt$gw)), colour = "#D94701", linewidth = 0.3) +
  geom_hline(yintercept = c(-log10(opt$sugg), log10(opt$sugg)), colour = "#2171B5", linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  scale_colour_manual(values = c("FALSE.new"="#2171B5","TRUE.new"="#6BAED6",
                                 "FALSE.scc"="#A63603","TRUE.scc"="#FD8D3C"), guide = "none") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr, expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(-ymax*1.12, ymax*1.12),
                     labels = function(v) format(abs(v), digits = 2)) +
  labs(title = sprintf("%s — corrected meta (v2, top) vs SCC manuscript meta (bottom)", ct),
       subtitle = sprintf("new: %s variants, %d GW, %d loci  |  SCC: %s variants, %d GW, %d loci  |  r(beta)=%.3f",
                          format(nrow(new), big.mark=","), summ$gw_new, summ$loci_new,
                          format(nrow(scc), big.mark=","), summ$gw_scc, summ$loci_scc, summ$r_beta),
       x = "Chromosome", y = expression(-log[10](P)~~"     new  /  SCC     ")) +
  theme_classic(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x = element_text(size = 7))

lab <- rbind(
  if (nrow(ln)) data.table(x = ln$pos + off$cum[match(ln$chr, off$chr)], y = -log10(ln$p), snp = ln$snp)[1:min(3,.N)] else NULL,
  if (nrow(ls_)) data.table(x = ls_$pos + off$cum[match(ls_$chr, off$chr)], y = log10(ls_$p), snp = ls_$snp)[1:min(3,.N)] else NULL)
if (!is.null(lab) && nrow(lab)) {
  if (requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(data = lab, aes(x, y, label = snp), inherit.aes = FALSE,
                                      size = 2.3, min.segment.length = 0, max.overlaps = 20)
}
ggsave(file.path(opt$out_dir, paste0(ct, "_miami.png")), p, width = 13, height = 7, dpi = 200)
cat(sprintf("  wrote %s_miami.png / _leads.tsv / _summary.tsv\n", ct))
