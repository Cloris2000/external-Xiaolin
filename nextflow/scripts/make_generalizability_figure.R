#!/usr/bin/env Rscript
# =============================================================================
# Generalizability of lead CTP associations across ancestries and cohorts
# -----------------------------------------------------------------------------
# Reads the plot-ready TSVs produced by scripts/select_ld_based_leads.py and
# builds a multi-panel figure:
#   Panel A  Ancestry forest   (Pooled / EUR / AFR / AMR effect + 95% CI per lead)
#   Panel B  Cohort matrix     (per-cohort harmonized beta; diverging fill)
#   Panel C  Two locus examples (cohort forest + ancestry subtotals + FE/RE)
#
# Outputs (results/.../generalizability/figures/ and copied to manuscript_figure/):
#   generalizability_compact.{pdf,svg,png}        (A + B)
#   generalizability_with_examples.{pdf,svg,png}  (A + B + C)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(scales)
  library(grid)
})

ROOT    <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
GEN_DIR <- file.path(ROOT, "results/meta_sensitivity/generalizability")
FIG_DIR <- file.path(GEN_DIR, "figures")
MS_DIR  <- file.path(ROOT, "manuscript_figure")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── palette (shared with other manuscript figures) ───────────────────────────
CLR_POOLED <- "#000000"
CLR_EUR    <- "#2166AC"
CLR_AFR    <- "#D95F02"
CLR_AMR    <- "#1B9E77"
CLR_LOW    <- "#2166AC"   # negative beta (blue)
CLR_MID    <- "#F7F7F7"
CLR_HIGH   <- "#B2182B"   # positive beta (red)
NA_GREY    <- "grey88"

BASE_SIZE  <- 8
theme_gen <- function() {
  theme_minimal(base_size = BASE_SIZE, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey90"),
      axis.text = element_text(colour = "grey20"),
      plot.title = element_text(face = "bold", size = BASE_SIZE + 2),
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = BASE_SIZE - 1),
      legend.title = element_text(size = BASE_SIZE - 1)
    )
}

# ── cohort order / labels ────────────────────────────────────────────────────
COHORT_ORDER <- c(
  "ROSMAP", "ROSMAP_array", "Mayo", "MSBB", "CMC_MSSM", "CMC_PENN", "CMC_PITT",
  "GTEx_v10", "NABEC", "GVEX", "NIMH_HBCC_Omni5M",
  "NIMH_HBCC_1M", "NIMH_HBCC_h650", "AMP_AD_Rush", "AMP_AD_Mayo"
)
COHORT_LABEL <- c(
  ROSMAP = "ROSMAP (WGS)", ROSMAP_array = "ROSMAP (array)", Mayo = "Mayo",
  MSBB = "MSBB", CMC_MSSM = "CMC MSSM", CMC_PENN = "CMC PENN",
  CMC_PITT = "CMC PITT", GTEx_v10 = "GTEx", NABEC = "NABEC", GVEX = "GVEX",
  NIMH_HBCC_Omni5M = "HBCC (Omni5M)", NIMH_HBCC_1M = "HBCC (1M)",
  NIMH_HBCC_h650 = "HBCC (h650)", AMP_AD_Rush = "AMP-AD Rush",
  AMP_AD_Mayo = "AMP-AD Mayo"
)

# Default detailed examples for Panel C (stable + heterogeneous). Fall back to
# data-driven picks if these leads are absent after LD clumping.
EXAMPLES_DEFAULT <- list(
  stable        = list(cell_type = "Microglia", lead = "chr16:31298939:T:G"),
  heterogeneous = list(cell_type = "L5.ET",     lead = "chr2:10862188:G:A")
)

# ── load data ────────────────────────────────────────────────────────────────
panelA <- fread(file.path(GEN_DIR, "generalizability_panelA.tsv"))
panelB <- fread(file.path(GEN_DIR, "generalizability_panelB_matrix.tsv"))
annot  <- fread(file.path(GEN_DIR, "lead_heterogeneity_sensitivity.tsv"))
long   <- fread(file.path(GEN_DIR, "lead_effects_long.tsv"))

# unique row id + display label
mk_rowlab <- function(cell_type, rsID, lead) {
  id <- ifelse(rsID != "" & !is.na(rsID), rsID, lead)
  paste0(cell_type, "  ", id)
}
panelA[, rowlab := mk_rowlab(cell_type, rsID, lead_variant)]
panelA <- panelA[order(row_order)]
row_levels <- panelA$rowlab              # top-to-bottom order (row_order 0 = top)
panelA[, rowlab := factor(rowlab, levels = rev(row_levels))]

annot[, rowlab := mk_rowlab(cell_type, rsID, lead_variant)]
setkey(annot, cell_type, lead_variant)

# =============================================================================
# Panel A : ancestry forest
# =============================================================================
build_panelA <- function() {
  strata <- c("pooled", "EUR", "AFR", "AMR")
  d <- rbindlist(lapply(strata, function(s) {
    data.table(
      rowlab = panelA$rowlab,
      row_order = panelA$row_order,
      stratum = s,
      beta = panelA[[paste0(ifelse(s == "pooled", "pooled", s), "_beta")]],
      se   = panelA[[paste0(ifelse(s == "pooled", "pooled", s), "_se")]]
    )
  }))
  d <- d[!is.na(beta)]
  d[, ci_lo := beta - 1.96 * se]
  d[, ci_hi := beta + 1.96 * se]
  d[, stratum := factor(stratum, levels = c("pooled", "EUR", "AFR", "AMR"),
                        labels = c("Pooled", "EUR", "AFR", "AMR"))]

  # ancestry-heterogeneity side annotation
  ann <- panelA[, .(rowlab, ancestry_i2, n_ancestry_groups,
                    cross_ancestry_assessable)]
  ann[, lab := ifelse(cross_ancestry_assessable == "yes",
                      sprintf("I²=%.0f%%", ancestry_i2),
                      "n/a")]
  xr <- range(c(d$ci_lo, d$ci_hi), na.rm = TRUE)
  xpad <- diff(xr) * 0.04

  ggplot(d, aes(x = beta, y = rowlab, colour = stratum, shape = stratum)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   position = position_dodge(width = 0.7), height = 0,
                   linewidth = 0.45) +
    geom_point(aes(size = stratum),
               position = position_dodge(width = 0.7)) +
    geom_text(data = ann, inherit.aes = FALSE,
              aes(x = xr[2] + xpad * 3, y = rowlab, label = lab),
              hjust = 0, size = 2.3, colour = "grey35") +
    scale_colour_manual(values = c(Pooled = CLR_POOLED, EUR = CLR_EUR,
                                   AFR = CLR_AFR, AMR = CLR_AMR), name = NULL) +
    scale_shape_manual(values = c(Pooled = 18, EUR = 16, AFR = 17, AMR = 15),
                       name = NULL) +
    scale_size_manual(values = c(Pooled = 2.4, EUR = 1.7, AFR = 1.7, AMR = 1.7),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    labs(x = "Effect size (beta, effect allele)", y = NULL,
         title = "A  Effect size across ancestries") +
    theme_gen() +
    theme(legend.position = "bottom",
          legend.margin = margin(t = -4),
          axis.text.y = element_text(size = BASE_SIZE - 2),
          plot.margin = margin(t = 5, r = 2, b = 2, l = 34))
}

# =============================================================================
# Panel B : cohort matrix
# =============================================================================
build_panelB <- function(show_ylab = FALSE) {
  d <- copy(panelB)
  d[, rowlab := mk_rowlab(cell_type, rsID, lead_variant)]
  d <- d[rowlab %in% row_levels]
  d[, rowlab := factor(rowlab, levels = rev(row_levels))]
  d[, cohort := factor(cohort, levels = COHORT_ORDER,
                       labels = COHORT_LABEL[COHORT_ORDER])]
  d[, beta_num := suppressWarnings(as.numeric(beta))]
  lim <- max(abs(d$beta_num), na.rm = TRUE)

  # side annotation: cohort I2 / n cohorts / % concordant / LOO stability
  ann <- annot[, .(rowlab, cohort_i2, cohort_k_extracted,
                   pct_cohorts_concordant, loo_sign_stable)]
  ann <- ann[rowlab %in% row_levels]
  ann[, rowlab := factor(rowlab, levels = rev(row_levels))]
  ann[, lab := sprintf("I²=%.0f%%  %.0f%%dir  k=%d%s",
                       cohort_i2, pct_cohorts_concordant, cohort_k_extracted,
                       ifelse(loo_sign_stable %in% c("True", TRUE), "  LOO ok", ""))]

  p <- ggplot(d, aes(x = cohort, y = rowlab, fill = beta_num)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = CLR_LOW, mid = CLR_MID, high = CLR_HIGH,
                         midpoint = 0, limits = c(-lim, lim),
                         na.value = NA_GREY, name = "Beta (effect allele)",
                         guide = guide_colourbar(barheight = unit(0.3, "cm"),
                                                 barwidth = unit(2.6, "cm"),
                                                 title.position = "top")) +
    labs(x = NULL, y = NULL, title = "B  Effect size across cohorts") +
    theme_gen() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 1),
          panel.grid.major = element_blank(),
          legend.position = "bottom",
          legend.margin = margin(t = -4))
  if (!show_ylab) {
    p <- p + theme(axis.text.y = element_blank())
  }
  list(plot = p, ann = ann)
}

# =============================================================================
# Panel C : detailed locus examples (cohort forest)
# =============================================================================
resolve_examples <- function() {
  present <- unique(paste(panelA$cell_type, panelA$lead_variant))
  pick <- list()
  for (nm in names(EXAMPLES_DEFAULT)) {
    e <- EXAMPLES_DEFAULT[[nm]]
    if (paste(e$cell_type, e$lead) %in% present) pick[[nm]] <- e
  }
  # data-driven fallbacks
  if (is.null(pick$stable)) {
    cand <- panelA[cross_ancestry_assessable == "yes"][order(ancestry_i2)][1]
    if (nrow(cand)) pick$stable <- list(cell_type = cand$cell_type, lead = cand$lead_variant)
  }
  if (is.null(pick$heterogeneous)) {
    cand <- panelA[cross_ancestry_assessable == "yes"][order(-ancestry_i2)][1]
    if (nrow(cand)) pick$heterogeneous <- list(cell_type = cand$cell_type, lead = cand$lead_variant)
  }
  pick
}

build_forest_one <- function(ct, lead_id, letter, subtitle) {
  a <- annot[cell_type == ct & lead_variant == lead_id]
  rs <- if (nrow(a) >= 1 && !is.na(a$rsID[1]) && a$rsID[1] != "") a$rsID[1] else lead_id

  # cohort rows
  ce <- long[cell_type == ct & lead_variant == lead_id &
             level == "cohort" & !is.na(beta)]
  ce[, ord := match(stratum, COHORT_ORDER)]
  ce <- ce[order(ord)]
  ce[, label := COHORT_LABEL[stratum]]
  ce[, kind := "Cohort"]
  ce[, anc := "cohort"]

  # ancestry subtotals + pooled + RE
  anc <- long[cell_type == ct & lead_variant == lead_id &
              level == "ancestry" & !is.na(beta)]
  anc[, label := stratum]
  anc[, kind := ifelse(stratum == "Pooled", "Pooled", "Ancestry")]
  anc[, anc := stratum]
  re_row <- data.table(label = "RE (DL)", beta = a$re_beta, se = a$re_se,
                       ci_lo = a$re_beta - 1.96 * a$re_se,
                       ci_hi = a$re_beta + 1.96 * a$re_se,
                       kind = "Pooled", anc = "RE")
  ce[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]
  anc[, `:=`(ci_lo = beta - 1.96 * se, ci_hi = beta + 1.96 * se)]

  ord_lab <- c(rev(ce$label), "", "Pooled", "EUR", "AFR", "AMR", "RE (DL)")
  ord_lab <- ord_lab[ord_lab %in% c(ce$label, anc$label, "RE (DL)", "")]

  df <- rbindlist(list(
    ce[, .(label, beta, ci_lo, ci_hi, kind, anc)],
    anc[, .(label, beta, ci_lo, ci_hi, kind, anc)],
    re_row[, .(label, beta, ci_lo, ci_hi, kind, anc)]
  ), use.names = TRUE, fill = TRUE)
  df <- df[!is.na(beta)]
  lev <- c("RE (DL)", "AMR", "AFR", "EUR", "Pooled", rev(ce$label))
  lev <- lev[lev %in% df$label]
  df[, label := factor(label, levels = lev)]
  df[, colgrp := fifelse(anc == "cohort", "Cohort",
                  fifelse(anc == "EUR", "EUR",
                   fifelse(anc == "AFR", "AFR",
                    fifelse(anc == "AMR", "AMR", "Pooled/RE"))))]

  i2c <- if (nrow(a)) a$cohort_i2 else NA
  i2a <- if (nrow(a)) a$ancestry_i2 else NA
  het_lab <- sprintf("cohort I²=%.0f%%   ancestry I²=%s",
                     i2c, ifelse(is.na(i2a), "n/a", sprintf("%.0f%%", i2a)))

  ggplot(df, aes(x = beta, y = label, colour = colgrp)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.5) +
    geom_point(aes(shape = colgrp, size = colgrp)) +
    annotate("text", x = Inf, y = Inf, label = het_lab, hjust = 1.02,
             vjust = 1.4, size = 2.2, colour = "grey35", fontface = "italic") +
    scale_colour_manual(values = c(Cohort = "grey45", EUR = CLR_EUR,
                                   AFR = CLR_AFR, AMR = CLR_AMR,
                                   `Pooled/RE` = CLR_POOLED), guide = "none") +
    scale_shape_manual(values = c(Cohort = 16, EUR = 16, AFR = 17, AMR = 15,
                                  `Pooled/RE` = 18), guide = "none") +
    scale_size_manual(values = c(Cohort = 1.6, EUR = 2.4, AFR = 2.4, AMR = 2.4,
                                 `Pooled/RE` = 2.8), guide = "none") +
    labs(x = "Effect size (beta)", y = NULL,
         title = paste0(letter, "  ", ct, "  ", rs),
         subtitle = subtitle) +
    theme_gen() +
    theme(plot.subtitle = element_text(size = BASE_SIZE - 1, colour = "grey35"))
}

# =============================================================================
# assemble + save
# =============================================================================
save_all <- function(plot, stem, w, h) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), plot, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(stem, ".svg")), plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), plot, width = w, height = h,
         dpi = 600)
  for (ext in c("pdf", "svg", "png")) {
    file.copy(file.path(FIG_DIR, paste0(stem, ".", ext)),
              file.path(MS_DIR, paste0(stem, ".", ext)), overwrite = TRUE)
  }
}

pA <- build_panelA()
pB_list <- build_panelB(show_ylab = FALSE)
pB <- pB_list$plot

n_rows <- length(row_levels)
fig_h <- max(4.5, 0.32 * n_rows + 2.2)

# compact = A + B (cowplot for robust panel alignment)
top_row <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1.45, 1.0),
                     align = "h", axis = "tb")
save_all(top_row, "generalizability_compact", w = 12.5, h = fig_h)

# with examples = (A + B) over C
ex <- resolve_examples()
c_plots <- list()
if (!is.null(ex$stable))
  c_plots[["stable"]] <- build_forest_one(ex$stable$cell_type, ex$stable$lead,
      "C", "Concordant across ancestries & cohorts")
if (!is.null(ex$heterogeneous))
  c_plots[["het"]] <- build_forest_one(ex$heterogeneous$cell_type,
      ex$heterogeneous$lead, "D", "Consistent direction, variable magnitude")

if (length(c_plots) >= 1) {
  panelC <- plot_grid(plotlist = c_plots, nrow = 1,
                      rel_widths = rep(1, length(c_plots)))
  full <- plot_grid(top_row, panelC, ncol = 1, rel_heights = c(fig_h, 3.4))
  save_all(full, "generalizability_with_examples", w = 12.5, h = fig_h + 3.6)
} else {
  save_all(top_row, "generalizability_with_examples", w = 12.5, h = fig_h)
}

cat("Figures written to", FIG_DIR, "and copied to", MS_DIR, "\n")
cat("Panel C examples:\n")
for (nm in names(ex)) cat("  ", nm, ":", ex[[nm]]$cell_type, ex[[nm]]$lead, "\n")
