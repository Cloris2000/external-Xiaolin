#!/usr/bin/env Rscript
# Figure 6 (v2): cross-cohort consistency of top CTP-GWAS loci -- data-driven.
#
# The original figure6_cohort_heterogeneity.R has the per-cohort betas typed in by
# hand for three loci from the old run (11 cohorts).  This version reads them:
#
#   per-cohort effects : the harmonized, hg19-lifted per-cohort files METAL actually
#                        consumed, <meta_dir>/harmonized/<ct>/harmonized/<cohort>_<ct>_harmonized.raw_p
#                        (so GTEx_v10 / AMP-AD / ROSMAP_array are included)
#   meta effect        : <meta_dir>/<ct>_meta_analysis_*.annotated.tsv (BETA_ALT, ALT-oriented)
#   ancestry           : docs/meta_cohort_metadata.tsv ancestry_class
#                        EUR_homogeneous -> "European"; AFR_enriched / mixed_ancestry -> "Mixed ancestries"
#
# Loci: by default the three manuscript loci; for each, the lead (min P) variant of
# the cell type within +/- window_kb in the v2 meta is used, so a lead that moved
# within the locus is followed.  Override with --loci <tsv: panel,cell_type,chr,pos,label>.
# All effects are per copy of the ALT allele (REGENIE ALLELE1; METAL BETA_ALT).
#
#   Rscript figure6_cohort_heterogeneity_v2.R --meta_dir ... --metadata docs/meta_cohort_metadata.tsv --out_dir ...
#
# Outputs: figure6_cohort_heterogeneity_v2.{png,pdf}, figure6_v2_loci.tsv, figure6_v2_per_cohort_effects.tsv

suppressPackageStartupMessages({
  library(data.table); library(optparse); library(ggplot2); library(dplyr); library(cowplot)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--meta_dir",  type = "character"),
  make_option("--metadata",  type = "character"),
  make_option("--out_dir",   type = "character"),
  make_option("--loci",      type = "character", default = NULL),
  make_option("--window_kb", type = "integer",   default = 500),
  make_option("--min_cohorts", type = "integer", default = 10,
              help = "Prefer the most significant window variant carried by at least this many cohorts (METAL Direction != '?'); falls back to the raw lead if none. 0 = raw lead. [default 10]")
)))
stopifnot(!is.null(opt$meta_dir), !is.null(opt$metadata), !is.null(opt$out_dir))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

CLR_EUR <- "#2166AC"; CLR_AFR <- "#D95F02"; CLR_META <- "#B2182B"; FIG_DPI <- 200

loci <- if (!is.null(opt$loci)) fread(opt$loci) else data.table(
  panel = c("A", "B", "C"),
  cell_type = c("VIP", "L5.6.IT.Car3", "Microglia"),
  chr = c(7L, 12L, 6L),
  pos = c(12284378L, 2324042L, 164862615L),
  label = c("TMEM106B", "CACNA1C", "PRKN region"))

meta_cohorts <- fread(opt$metadata)
anc_map <- setNames(ifelse(meta_cohorts$ancestry_class == "EUR_homogeneous", "EUR", "AFR"), meta_cohorts$cohort)

annotated_for <- function(ct) {
  f <- list.files(opt$meta_dir, pattern = sprintf("^%s_meta_analysis_.*\\.annotated\\.tsv$", gsub("\\.", "\\\\.", ct)), full.names = TRUE)
  if (!length(f)) stop("no annotated.tsv for ", ct); f[1]
}

# Lead within the window.  The annotated table is ~1.2 GB, so pre-filter with a
# fast grep on the MarkerName prefix (chr:<1 Mb bucket>) and do the exact window
# test in R; an awk field-parse of the whole file took ~10 min per locus.
find_lead <- function(ct, chr, pos, win) {
  f <- annotated_for(ct)
  lo <- max(1, pos - win); hi <- pos + win
  buckets <- seq(lo %/% 1e6, hi %/% 1e6)
  alts <- ifelse(buckets == 0, "[0-9]{1,6}", paste0(buckets, "[0-9]{6}"))
  rx <- sprintf("^chr%s:(%s):", chr, paste(alts, collapse = "|"))
  # `|| true`: grep exits 1 on no match, which fread(cmd=) would report as a failure.
  dt <- fread(cmd = sprintf("{ head -1 %s; grep -E %s %s || true; }", shQuote(f), shQuote(rx), shQuote(f)))
  dt <- dt[as.numeric(POS) >= lo & as.numeric(POS) <= hi]
  if (!nrow(dt)) stop(sprintf("%s: no variants within %d kb of chr%s:%d", ct, win / 1000, chr, pos))
  dt[, p := as.numeric(`P-value`)]
  dt[, n_coh := nchar(gsub("\\?", "", Direction))]
  raw_lead <- dt[which.min(p)]
  # A forest plot of cross-cohort consistency is most informative for a variant most
  # cohorts actually carry; the raw lead can be an indel present in only a few.
  if (opt$min_cohorts > 0 && any(dt$n_coh >= opt$min_cohorts)) {
    lead <- dt[n_coh >= opt$min_cohorts][which.min(p)]
    if (lead$MarkerName != raw_lead$MarkerName)
      cat(sprintf("    raw lead %s (P=%.3g, %d cohorts) -> using %s (P=%.3g, %d cohorts) [--min_cohorts %d]\n",
                  raw_lead$MarkerName, raw_lead$p, raw_lead$n_coh, lead$MarkerName, lead$p, lead$n_coh, opt$min_cohorts))
    lead[, raw_lead := raw_lead$MarkerName][, raw_lead_P := raw_lead$p]
    return(lead)
  }
  raw_lead[, raw_lead := MarkerName][, raw_lead_P := p]
}

# Per-cohort effect at one marker from the harmonized (hg19) file METAL consumed.
# Files are space-delimited with ID in column 3, so the marker is always " <id> ".
cohort_effect <- function(ct, cohort, marker) {
  # HARMONIZE_META_SUMSTATS publishes to <meta_dir>/harmonized/<ct>/harmonized/
  f <- file.path(opt$meta_dir, "harmonized", ct, "harmonized", sprintf("%s_%s_harmonized.raw_p", cohort, ct))
  if (!file.exists(f)) return(NULL)
  # A cohort may simply not carry this variant: grep then exits 1, so `|| true`.
  dt <- fread(cmd = sprintf("{ head -1 %s; grep -F -m1 -- %s %s || true; }", shQuote(f), shQuote(paste0(" ", marker, " ")), shQuote(f)))
  if (nrow(dt) < 1) return(NULL)
  data.table(label = cohort, beta = as.numeric(dt$BETA[1]), se = as.numeric(dt$SE[1]),
             n = as.integer(dt$N[1]), pval = as.numeric(dt$P[1]),
             effect_allele = toupper(dt$ALLELE1[1]), ancestry = anc_map[[cohort]], is_meta = FALSE)
}

compute_heterogeneity <- function(beta, se) {
  wi <- 1 / se^2; bhat <- sum(wi * beta) / sum(wi); Q <- sum(wi * (beta - bhat)^2)
  df <- length(beta) - 1L
  list(Q = Q, df = df, I2 = if (Q > 0) max(0, (Q - df) / Q) * 100 else 0, Qp = pchisq(Q, df, lower.tail = FALSE))
}

make_forest <- function(cohorts_df, meta_df, panel_letter, title_str, x_label) {
  # EUR first (by N desc), then mixed-ancestry cohorts, then the meta diamond
  ord <- cohorts_df %>% arrange(ancestry != "EUR", desc(n)) %>% pull(label)
  cohort_order <- c(ord, "Meta-analysis")
  n_eur <- sum(cohorts_df$ancestry == "EUR")
  df <- bind_rows(cohorts_df, meta_df) %>%
    mutate(label = factor(label, levels = rev(cohort_order)),
           ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se,
           plabel = ifelse(pval < 0.001, formatC(pval, format = "e", digits = 2), formatC(pval, format = "f", digits = 3)),
           nlabel = ifelse(is.na(n), "", paste0("N=", formatC(n, format = "d", big.mark = ","))),
           anc_label = factor(case_when(ancestry == "EUR" ~ "European", ancestry == "AFR" ~ "Mixed ancestries", TRUE ~ "Meta-analysis"),
                              levels = c("European", "Mixed ancestries", "Meta-analysis")))
  het <- compute_heterogeneity(cohorts_df$beta, cohorts_df$se)
  het_label <- sprintf("I² = %.1f%%  (Q = %.1f, df = %d, P = %.2g)", het$I2, het$Q, het$df, het$Qp)
  xr <- range(c(df$ci_lo, df$ci_hi), na.rm = TRUE); pad <- diff(xr) * 0.05
  x_left <- xr[1] - pad; x_right <- xr[2] + pad
  # dotted separators: between EUR / mixed blocks and above the meta row
  seps <- c(1.5, if (n_eur < nrow(cohorts_df)) nrow(cohorts_df) - n_eur + 1.5)
  p <- ggplot(df, aes(y = label, x = beta, colour = anc_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = seps, linetype = "dotted", colour = "grey65", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, linewidth = is_meta), height = 0.28) +
    geom_point(aes(shape = is_meta, size = is_meta)) +
    geom_text(aes(x = x_left - 0.02, label = nlabel), hjust = 1, size = 3.0, colour = "grey40") +
    geom_text(aes(x = x_right + 0.02, label = paste0("P=", plabel)), hjust = 0, size = 3.0, colour = "grey20") +
    scale_colour_manual(name = "Ancestry", values = c("European" = CLR_EUR, "Mixed ancestries" = CLR_AFR, "Meta-analysis" = CLR_META),
                        guide = guide_legend(override.aes = list(shape = 16, size = 3, linewidth = 0.6))) +
    scale_linewidth_manual(values = c("FALSE" = 0.7, "TRUE" = 1.4), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 4.2), guide = "none") +
    scale_x_continuous(expand = expansion(add = c(0.30, 0.44))) +
    # heterogeneity goes in the subtitle: inside the panel it collides with the top P label
    labs(title = paste0(panel_letter, ". ", title_str), subtitle = het_label, x = x_label, y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
          plot.subtitle = element_text(size = 9.5, colour = "grey30", face = "italic", hjust = 0),
          axis.text.y = element_text(size = 10, colour = "black"), axis.text.x = element_text(size = 9, colour = "black"),
          axis.title.x = element_text(size = 10), panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          legend.position = "none", plot.margin = margin(8, 8, 8, 8))
  list(plot = p, het = het)
}

panels <- list(); lead_rows <- list(); effect_rows <- list()
for (i in seq_len(nrow(loci))) {
  L <- loci[i]
  lead <- find_lead(L$cell_type, L$chr, L$pos, opt$window_kb * 1000L)
  marker <- lead$MarkerName
  cat(sprintf("Panel %s  %s x %s: lead %s  P=%.3g  (%s)\n", L$panel, L$cell_type, L$label, marker, lead$p,
              if (lead$p < 1e-5) "suggestive" else "NOT suggestive in v2"))
  eff <- rbindlist(lapply(meta_cohorts$cohort, function(co) cohort_effect(L$cell_type, co, marker)))
  eff <- eff[!is.na(beta) & !is.na(se) & se > 0]
  meta_row <- data.table(label = "Meta-analysis", beta = as.numeric(lead$BETA_ALT), se = as.numeric(lead$StdErr),
                         n = NA_integer_, pval = lead$p, effect_allele = lead$ALT, ancestry = "META", is_meta = TRUE)
  title <- sprintf("%s × %s (%s, effect allele = %s)", L$cell_type, L$label, marker, lead$ALT)
  fp <- make_forest(eff, meta_row, L$panel, title, sprintf("Effect on %s proportion (β per ALT allele, 95%% CI)", L$cell_type))
  panels[[L$panel]] <- fp$plot
  lead_rows[[i]] <- data.table(panel = L$panel, cell_type = L$cell_type, requested_chr = L$chr, requested_pos = L$pos,
                               label = L$label, lead = marker, lead_P = lead$p, ALT = lead$ALT, BETA_ALT = lead$BETA_ALT,
                               n_cohorts = nrow(eff), I2 = round(fp$het$I2, 1), Q = round(fp$het$Q, 2), Q_df = fp$het$df, Q_P = fp$het$Qp,
                               suggestive_v2 = lead$p < 1e-5,
                               raw_window_lead = lead$raw_lead, raw_window_lead_P = lead$raw_lead_P)
  effect_rows[[i]] <- cbind(panel = L$panel, cell_type = L$cell_type, marker = marker, eff)
}
fwrite(rbindlist(lead_rows), file.path(opt$out_dir, "figure6_v2_loci.tsv"), sep = "\t")
fwrite(rbindlist(effect_rows), file.path(opt$out_dir, "figure6_v2_per_cohort_effects.tsv"), sep = "\t")

legend_df <- data.frame(x = 1:3, y = 1, g = factor(c("European", "Mixed ancestries", "Meta-analysis"),
                                                    levels = c("European", "Mixed ancestries", "Meta-analysis")))
leg <- get_legend(ggplot(legend_df, aes(x, y, colour = g)) + geom_point(size = 3) +
  scale_colour_manual(name = "Ancestry", values = c("European" = CLR_EUR, "Mixed ancestries" = CLR_AFR, "Meta-analysis" = CLR_META)) +
  theme_classic() + theme(legend.position = "bottom"))
fig <- plot_grid(plot_grid(plotlist = panels, ncol = 1, align = "v"), leg, ncol = 1, rel_heights = c(1, 0.05))
h <- 3.6 * length(panels) + 0.6
ggsave(file.path(opt$out_dir, "figure6_cohort_heterogeneity_v2.png"), fig, width = 11, height = h, dpi = FIG_DPI, bg = "white")
ggsave(file.path(opt$out_dir, "figure6_cohort_heterogeneity_v2.pdf"), fig, width = 11, height = h, device = cairo_pdf)
cat("Saved figure6_cohort_heterogeneity_v2.{png,pdf} to", opt$out_dir, "\n")
