#!/usr/bin/env Rscript
# =============================================================================
# Generalizability figure (revised: Panels A + B + C) -- publication style
# -----------------------------------------------------------------------------
# A  Ancestry-stratified effects  (forest; ONLY leads with >=2 ancestry estimates)
# B  Cohort-level effects          (heatmap; ONLY leads present in >=7 cohorts)
# C  Cohort-level forest plots for 3 representative loci (1 row x 3 col)
#
# Panels A and B are DECOUPLED: each shows only the leads for which its own
# heterogeneity comparison is meaningful, and each carries its own row labels.
#
# Uses ONLY the plot-ready tables produced by scripts/select_ld_based_leads.py.
# SNP selection is NEVER recomputed here. LOCO / RE are derived (meta-analysis
# arithmetic) from the already-harmonized cohort betas/SEs (full LOCO in supp).
#
# Display orientation is controlled by ORIENT_TO_INCREASING (see config below):
#   TRUE  = flip to pooled CTP-increasing allele (color = concordance with meta);
#   FALSE = keep ALT allele (color = biological direction of the ALT allele).
# All flips are recorded in figure_effect_orientation_audit.tsv.
#
# The original scripts/make_generalizability_figure.R is preserved unchanged.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
})

MIN_ANCESTRIES <- 2   # Panel A inclusion
MIN_COHORTS    <- 7   # Panel B inclusion

# Display orientation:
#   TRUE  -> flip each SNP to the pooled CTP-increasing allele (color = concordance
#            with the meta-analysis; pooled beta always positive; heatmap ~all red).
#   FALSE -> keep the ALT allele as the displayed effect allele (color = biological
#            direction of the ALT allele: red = raises CTP, blue = lowers CTP).
# Cohorts are harmonized to the same per-SNP allele in either mode.
ORIENT_TO_INCREASING <- FALSE

# Panel C examples (VARIANT "" -> auto; candidate table printed)
EXAMPLE_CELLTYPES <- c("Microglia",          "L6b",                "L5.ET")
EXAMPLE_VARIANTS  <- c("chr16:31298939:T:G",  "chr3:142612091:C:T", "chr2:10862188:G:A")
EXAMPLE_QUALITY   <- c("stable",             "intermediate",       "heterogeneous")

ROOT    <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
GEN_DIR <- file.path(ROOT, "results/meta_sensitivity/generalizability")
OUT_DIR <- file.path(GEN_DIR, "figures")
TAB_DIR <- file.path(GEN_DIR, "figure_tables")
MS_DIR  <- file.path(ROOT, "manuscript_figure")
for (d in c(OUT_DIR, TAB_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

CLR_POOLED <- "#000000"; CLR_EUR <- "#2166AC"; CLR_AFR <- "#D95F02"
CLR_AMR <- "#1B9E77"
CLR_LOW <- "#2166AC"; CLR_MID <- "#F7F7F7"; CLR_HIGH <- "#B2182B"
NA_GREY <- "grey85"
COMP_COL <- c(`EUR-majority` = "#0072B2", `AFR-enriched` = "#D55E00",
              `Mixed/diverse` = "#009E73")

BASE <- 10
theme_base <- function() {
  theme_minimal(base_size = BASE, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text = element_text(colour = "grey20"),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
}

COHORT_ORDER <- c(
  "ROSMAP", "ROSMAP_array", "Mayo", "MSBB", "CMC_MSSM", "CMC_PENN", "CMC_PITT",
  "GTEx_v10", "NABEC", "GVEX", "NIMH_HBCC_Omni5M",
  "NIMH_HBCC_1M", "NIMH_HBCC_h650", "AMP_AD_Rush", "AMP_AD_Mayo")
COHORT_LABEL <- c(
  ROSMAP = "ROSMAP (WGS)", ROSMAP_array = "ROSMAP (array)", Mayo = "Mayo",
  MSBB = "MSBB", CMC_MSSM = "CMC MSSM", CMC_PENN = "CMC PENN",
  CMC_PITT = "CMC PITT", GTEx_v10 = "GTEx", NABEC = "NABEC", GVEX = "GVEX",
  NIMH_HBCC_Omni5M = "HBCC (Omni5M)", NIMH_HBCC_1M = "HBCC (1M)",
  NIMH_HBCC_h650 = "HBCC (h650)", AMP_AD_Rush = "AMP-AD Rush",
  AMP_AD_Mayo = "AMP-AD Mayo")
COHORT_COMP <- c(
  ROSMAP = "EUR-majority", ROSMAP_array = "EUR-majority", Mayo = "EUR-majority",
  MSBB = "EUR-majority", CMC_MSSM = "EUR-majority", CMC_PENN = "EUR-majority",
  CMC_PITT = "EUR-majority", GTEx_v10 = "EUR-majority", NABEC = "EUR-majority",
  GVEX = "EUR-majority", NIMH_HBCC_Omni5M = "EUR-majority",
  NIMH_HBCC_1M = "AFR-enriched", NIMH_HBCC_h650 = "AFR-enriched",
  AMP_AD_Rush = "AFR-enriched", AMP_AD_Mayo = "Mixed/diverse")
FAMILY_BREAKS <- c(2, 7, 10, 11, 13)

fe_meta <- function(b, s) { w <- 1 / s^2; list(b = sum(w * b) / sum(w),
                                                se = sqrt(1 / sum(w))) }
het_stats <- function(b, s) {
  if (length(b) < 2) return(list(Q = NA, df = NA, I2 = NA, p = NA))
  w <- 1 / s^2; bb <- sum(w * b) / sum(w); Q <- sum(w * (b - bb)^2)
  df <- length(b) - 1; I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  list(Q = Q, df = df, I2 = I2, p = pchisq(Q, df, lower.tail = FALSE))
}
dl_re <- function(b, s) {
  if (length(b) < 2) { m <- fe_meta(b, s); return(list(b = m$b, se = m$se)) }
  w <- 1 / s^2; bb <- sum(w * b) / sum(w); Q <- sum(w * (b - bb)^2)
  df <- length(b) - 1; C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (C > 0) max(0, (Q - df) / C) else 0
  wr <- 1 / (s^2 + tau2)
  list(b = sum(wr * b) / sum(wr), se = sqrt(1 / sum(wr)))
}

annot <- fread(file.path(GEN_DIR, "lead_heterogeneity_sensitivity.tsv"))
long  <- fread(file.path(GEN_DIR, "lead_effects_long.tsv"))
setorder(annot, row_order)
CLASS_LV <- c("Non-neuronal", "Inhibitory", "Excitatory")

sv <- tstrsplit(annot$lead_variant, ":", fixed = TRUE)
annot[, `:=`(REF = sv[[3]], ALT = sv[[4]])]

annot[, orient_sign := if (ORIENT_TO_INCREASING) ifelse(pooled_beta < 0, -1, 1) else 1]
annot[, original_effect_allele  := ALT]
annot[, displayed_effect_allele := ifelse(orient_sign < 0, REF, ALT)]
annot[, displayed_pooled_beta   := pooled_beta * orient_sign]

fwrite(annot[, .(cell_type, lead_variant, rsID, original_effect_allele,
  displayed_effect_allele, original_pooled_beta = pooled_beta, displayed_pooled_beta,
  orientation_flipped = orient_sign < 0,
  notes = ifelse(orient_sign < 0, "beta*-1; alleles swapped for display", "none"))],
  file.path(TAB_DIR, "figure_effect_orientation_audit.tsv"), sep = "\t")

# label convention: celltype_rsID (or celltype_chr:pos:REF:ALT if no rsID)
disp_label <- function(ct, rs, var) {
  id <- ifelse(!is.na(rs) & rs != "", rs, var); paste0(ct, "_", id) }
annot[, rowlab := disp_label(cell_type, rsID, lead_variant)]

sign_map <- setNames(annot$orient_sign, paste(annot$cell_type, annot$lead_variant))
long[, k := paste(cell_type, lead_variant)]
long[, os := sign_map[k]]
long[, dbeta := beta * os]
long[, dlo := ifelse(is.na(beta), NA, ifelse(os < 0, -ci_hi, ci_lo))]
long[, dhi := ifelse(is.na(beta), NA, ifelse(os < 0, -ci_lo, ci_hi))]

mk_panel <- function(a) { a <- copy(a); setorder(a, row_order)
  a[, rank := seq_len(.N)]; n <- nrow(a); a[, ypos := n - rank + 1]; a }
sep_of <- function(a) { n <- nrow(a)
  last <- sapply(CLASS_LV, function(cl) { r <- a[broad_class == cl, rank]
    if (length(r)) max(r) else NA }); last <- last[!is.na(last)]
  if (length(last) > 1) n - head(last, -1) + 0.5 else numeric(0) }
annotA <- mk_panel(annot[n_ancestry_groups >= MIN_ANCESTRIES]); nA <- nrow(annotA)
annotB <- mk_panel(annot[cohort_k_extracted >= MIN_COHORTS]);  nB <- nrow(annotB)
sepA <- sep_of(annotA); sepB <- sep_of(annotB)

# =============================================================================
# PANEL A : ancestry-stratified forest
# =============================================================================
ANC_LV  <- c("Pooled", "EUR", "AFR", "AMR")
ANC_LAB <- c(Pooled = "Meta-analysis", EUR = "EUR", AFR = "AFR", AMR = "AMR")
ANC_OFF <- c(Pooled = 0.26, EUR = 0.09, AFR = -0.09, AMR = -0.26)
ANC_COL <- c(Pooled = CLR_POOLED, EUR = CLR_EUR, AFR = CLR_AFR, AMR = CLR_AMR)
ANC_SHP <- c(Pooled = 18, EUR = 16, AFR = 17, AMR = 15)

a_dat <- long[level == "ancestry" & stratum %in% ANC_LV & !is.na(dbeta)]
a_dat <- merge(a_dat, annotA[, .(cell_type, lead_variant, ypos)],
               by = c("cell_type", "lead_variant"))
a_dat[, stratum := factor(stratum, levels = ANC_LV)]
a_dat[, yy := ypos + ANC_OFF[as.character(stratum)]]

xr <- range(c(a_dat$dlo, a_dat$dhi), na.rm = TRUE); xw <- diff(xr)
YA_HI <- nA + 0.7

pA <- ggplot(a_dat) +
  {if (length(sepA)) geom_hline(yintercept = sepA, colour = "grey75", linewidth = 0.3)} +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(y = yy, xmin = dlo, xmax = dhi, colour = stratum),
                 height = 0, linewidth = 0.45) +
  geom_point(aes(y = yy, x = dbeta, colour = stratum, shape = stratum, size = stratum)) +
  scale_colour_manual(values = ANC_COL, labels = ANC_LAB, name = NULL, drop = FALSE) +
  scale_shape_manual(values = ANC_SHP, labels = ANC_LAB, name = NULL, drop = FALSE) +
  scale_size_manual(values = c(Pooled = 2.4, EUR = 1.7, AFR = 1.7, AMR = 1.7), guide = "none") +
  scale_x_continuous(limits = c(xr[1] - xw * 0.05, xr[2] + xw * 0.05),
                     breaks = scales::pretty_breaks(4)(xr), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, YA_HI), expand = c(0, 0),
                     breaks = annotA$ypos, labels = annotA$rowlab) +
  labs(x = expression("Effect size ("*beta*")"), y = NULL) +
  theme_base() +
  theme(legend.position = "bottom", legend.margin = margin(t = -4),
        axis.text.y = element_text(size = BASE - 2),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 5, r = 6, b = 2, l = 16))

# =============================================================================
# PANEL B : cohort heatmap (own labels; composition legend under B)
# =============================================================================
b_dat <- long[level == "cohort"]
b_dat <- merge(b_dat, annotB[, .(cell_type, lead_variant, ypos)],
               by = c("cell_type", "lead_variant"))
b_dat[, cohort := factor(stratum, levels = COHORT_ORDER)]
b_dat <- b_dat[!is.na(cohort)]; b_dat[, cx := as.integer(cohort)]

allbeta <- b_dat[!is.na(dbeta), dbeta]; lim_full <- max(abs(allbeta))
lim_use <- ceiling(lim_full / 0.05) * 0.05   # full symmetric range (no clipping)
clipped <- b_dat[!is.na(dbeta) & abs(dbeta) > lim_use]
b_dat[, fill_beta := pmax(pmin(dbeta, lim_use), -lim_use)]

nC_coh <- length(COHORT_ORDER)
COMP_Y <- nB + 1.0; YB_HI <- nB + 1.8
comp_dt <- data.table(cx = seq_len(nC_coh),
                      comp = factor(COHORT_COMP[COHORT_ORDER], levels = names(COMP_COL)))

pB <- ggplot() +
  geom_tile(data = b_dat, aes(x = cx, y = ypos, fill = fill_beta),
            colour = "white", linewidth = 0.35) +
  geom_tile(data = comp_dt, aes(x = cx, y = COMP_Y),
            fill = COMP_COL[as.character(comp_dt$comp)], height = 0.7, width = 0.9) +
  # invisible layer to build the composition legend near Panel B
  geom_point(data = comp_dt, aes(x = cx, y = COMP_Y, colour = comp), alpha = 0) +
  {if (length(sepB)) geom_hline(yintercept = sepB, colour = "grey75", linewidth = 0.3)} +
  geom_vline(xintercept = FAMILY_BREAKS + 0.5, colour = "grey80", linewidth = 0.3) +
  scale_fill_gradient2(low = CLR_LOW, mid = CLR_MID, high = CLR_HIGH, midpoint = 0,
                       limits = c(-lim_use, lim_use), na.value = NA_GREY,
                       name = expression("Standardized effect size, "*beta),
                       guide = guide_colourbar(barheight = unit(0.3, "cm"),
                                               barwidth = unit(2.8, "cm"),
                                               title.position = "top", order = 1)) +
  scale_colour_manual(values = COMP_COL, name = NULL,
                      guide = guide_legend(override.aes = list(alpha = 1, size = 3,
                                           shape = 15), order = 2)) +
  scale_x_continuous(breaks = seq_len(nC_coh), labels = COHORT_LABEL[COHORT_ORDER],
                     limits = c(0.4, nC_coh + 0.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, YB_HI), expand = c(0, 0),
                     breaks = annotB$ypos, labels = annotB$rowlab) +
  labs(x = NULL, y = NULL) +
  theme_base() +
  theme(axis.text.y = element_text(size = BASE - 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = BASE - 1),
        panel.grid.major = element_blank(),
        legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(t = -2), legend.spacing.x = unit(0.8, "cm"),
        legend.box.spacing = unit(0.4, "cm"),
        plot.margin = margin(t = 5, r = 6, b = 2, l = 6)) +
  coord_cartesian(clip = "off")

# =============================================================================
# PANEL C : 3 representative cohort forest plots (1 row x 3 col)
# =============================================================================
build_forest <- function(ct, var) {
  ar <- annot[cell_type == ct & lead_variant == var]
  rs <- if (nrow(ar) && !is.na(ar$rsID) && ar$rsID != "") ar$rsID else var
  coh <- long[cell_type == ct & lead_variant == var & level == "cohort" & !is.na(dbeta)]
  coh[, comp := factor(COHORT_COMP[as.character(stratum)], levels = names(COMP_COL))]
  coh[, ord := match(stratum, COHORT_ORDER)]; setorder(coh, ord)
  b <- coh$dbeta; s <- coh$se
  fe <- fe_meta(b, s); re <- dl_re(b, s); hs <- het_stats(b, s); k <- length(b)
  anc <- long[cell_type == ct & lead_variant == var & level == "ancestry" &
              stratum %in% c("EUR", "AFR", "AMR") & !is.na(dbeta)]
  rows <- rbindlist(list(
    coh[, .(label = COHORT_LABEL[stratum], beta = dbeta, lo = dlo, hi = dhi,
            grp = as.character(comp), kind = "cohort")],
    data.table(label = "Meta", beta = fe$b, lo = fe$b - 1.96 * fe$se,
               hi = fe$b + 1.96 * fe$se, grp = "summary", kind = "fe")),
    use.names = TRUE)
  rows[, y := .N:1]
  xrng <- range(c(rows$lo, rows$hi, 0), na.rm = TRUE); xp <- diff(xrng) * 0.04
  grp_col <- c(COMP_COL, anc = "grey40", summary = "black")
  grp_shp <- c(`EUR-majority` = 16, `AFR-enriched` = 17, `Mixed/diverse` = 15,
               anc = 18, summary = 18)
  p <- ggplot(rows, aes(x = beta, y = y)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = grp), height = 0, linewidth = 0.5) +
    geom_point(aes(colour = grp, shape = grp, size = kind)) +
    scale_colour_manual(values = grp_col, guide = "none") +
    scale_shape_manual(values = grp_shp, guide = "none") +
    scale_size_manual(values = c(cohort = 1.8, ancestry = 2.6, fe = 3, re = 3), guide = "none") +
    scale_y_continuous(breaks = rows$y, labels = rows$label, expand = c(0, 0.6)) +
    coord_cartesian(xlim = c(xrng[1] - xp, xrng[2] + xp)) +
    labs(x = expression("Effect size ("*beta*")"), y = NULL, title = rs) +
    theme_base() +
    theme(plot.title = element_text(size = BASE + 1, hjust = 0.5, face = "plain"),
          axis.text.y = element_text(size = BASE - 1),
          panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
          plot.margin = margin(t = 4, r = 8, b = 2, l = 4))
  list(plot = p, fdat = rows, rs = rs, I2 = hs$I2, k = k, fe = fe, re = re,
       cell_type = ct, lead_variant = var)
}
examples <- lapply(seq_along(EXAMPLE_CELLTYPES), function(i)
  build_forest(EXAMPLE_CELLTYPES[i], EXAMPLE_VARIANTS[i]))

# =============================================================================
# LOCO for ALL leads (supplementary only)
# =============================================================================
loco_all <- rbindlist(lapply(seq_len(nrow(annot)), function(i) {
  ct <- annot$cell_type[i]; var <- annot$lead_variant[i]
  coh <- long[cell_type == ct & lead_variant == var & level == "cohort" & !is.na(dbeta)]
  coh[, ord := match(stratum, COHORT_ORDER)]; setorder(coh, ord)
  b <- coh$dbeta; s <- coh$se; k <- length(b); if (k < 1) return(NULL)
  fe <- fe_meta(b, s)
  out <- data.table(cell_type = ct, lead_variant = var, rsID = annot$rsID[i],
                    omitted_cohort = "None (all cohorts)", loco_beta = fe$b,
                    loco_se = fe$se, is_reference = TRUE)
  if (k >= 2) for (j in seq_len(k)) { m <- fe_meta(b[-j], s[-j])
    out <- rbind(out, data.table(cell_type = ct, lead_variant = var, rsID = annot$rsID[i],
      omitted_cohort = coh$stratum[j], loco_beta = m$b, loco_se = m$se, is_reference = FALSE)) }
  out[, `:=`(loco_ci_lo = loco_beta - 1.96 * loco_se, loco_ci_hi = loco_beta + 1.96 * loco_se,
             pooled_direction_retained = all(sign(out$loco_beta) == sign(fe$b)))]; out
}))
fwrite(loco_all, file.path(TAB_DIR, "supplementary_LOCO_all_leads.tsv"), sep = "\t")
loco_all[, facet := disp_label(cell_type, rsID, lead_variant)]
loco_all[, ylab := ifelse(is_reference, "ALL", COHORT_LABEL[as.character(omitted_cohort)])]
loco_all[, yy := seq_len(.N), by = facet]
supp <- ggplot(loco_all, aes(x = loco_beta, y = reorder(ylab, -yy))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = loco_ci_lo, xmax = loco_ci_hi, colour = is_reference),
                 height = 0, linewidth = 0.4) +
  geom_point(aes(colour = is_reference, size = is_reference)) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey45"), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 1.8, `FALSE` = 1.0), guide = "none") +
  facet_wrap(~ facet, scales = "free", ncol = 4) +
  labs(x = expression("Leave-one-cohort-out effect size ("*beta*")"), y = NULL) +
  theme_base() +
  theme(strip.text = element_text(size = BASE - 2, face = "bold"),
        axis.text.y = element_text(size = BASE - 3), axis.text.x = element_text(size = BASE - 2),
        panel.spacing = unit(0.4, "lines"))
ggsave(file.path(OUT_DIR, "supplementary_LOCO_all_leads.pdf"), supp,
       width = 14, height = 2.4 * ceiling(nrow(annot) / 4), device = cairo_pdf,
       limitsize = FALSE, bg = "white")

# =============================================================================
# ASSEMBLE (panel letters only; no descriptive titles)
# =============================================================================
top_row <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1.18, 1.55),
                     labels = c("A", "B"), label_size = 15, label_fontface = "bold")
cRow <- plot_grid(plotlist = lapply(examples, `[[`, "plot"), nrow = 1)
top_h <- max(5.5, 0.32 * max(nA, nB) + 1.8)
c_h   <- 3.8
fig_ABC <- plot_grid(top_row, cRow, ncol = 1, rel_heights = c(top_h, c_h),
                     labels = c("", "C"), label_size = 15, label_fontface = "bold")

save_all <- function(p, stem, w, h) {
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, limitsize = FALSE, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(stem, ".svg")), p, width = w, height = h,
         limitsize = FALSE, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 600, limitsize = FALSE, bg = "white")
  for (ext in c("pdf", "svg", "png"))
    file.copy(file.path(OUT_DIR, paste0(stem, ".", ext)),
              file.path(MS_DIR, paste0(stem, ".", ext)), overwrite = TRUE)
}
save_all(fig_ABC, "generalizability_figure_ABC", w = 14, h = top_h + c_h)
save_all(top_row, "generalizability_figure_AB",  w = 14, h = top_h)

# =============================================================================
# PLOT-READY TABLES
# =============================================================================
fwrite(a_dat[, .(cell_type, lead_variant, stratum = as.character(stratum),
                 displayed_beta = dbeta, ci_lo = dlo, ci_hi = dhi, se, p, status)],
       file.path(TAB_DIR, "panel_A_ancestry_effects.tsv"), sep = "\t")
fwrite(b_dat[, .(cell_type, lead_variant, cohort = as.character(cohort),
                 cohort_composition = COHORT_COMP[as.character(cohort)],
                 displayed_beta = dbeta, fill_beta, se, status)],
       file.path(TAB_DIR, "panel_B_cohort_effects.tsv"), sep = "\t")
psgn <- annot[, .(cell_type, lead_variant, psgn = sign(displayed_pooled_beta))]
dir_dt <- merge(b_dat[!is.na(dbeta)], psgn, by = c("cell_type", "lead_variant"))
dir_dt <- dir_dt[, .(concordant = sum(sign(dbeta) == psgn), contributing = .N),
                 by = .(cell_type, lead_variant)]
dir_dt <- merge(annotB[, .(cell_type, lead_variant, cohort_i2)], dir_dt,
                by = c("cell_type", "lead_variant"), all.x = TRUE)
fwrite(dir_dt, file.path(TAB_DIR, "panel_B_cohort_summaries.tsv"), sep = "\t")
fwrite(rbindlist(lapply(examples, function(e)
  cbind(rsID = e$rs, cell_type = e$cell_type, lead_variant = e$lead_variant,
        e$fdat[, .(label, beta, lo, hi, grp, kind)]))),
  file.path(TAB_DIR, "panel_C_forest_data.tsv"), sep = "\t")

# =============================================================================
# QC REPORT
# =============================================================================
qc <- c(); add <- function(...) qc <<- c(qc, sprintf(...))
add("Generalizability figure (A/B/C) QC report"); add(strrep("=", 60))
add("1. Script: scripts/make_generalizability_figure_ABC.R")
add("2. Inputs: lead_heterogeneity_sensitivity.tsv, lead_effects_long.tsv")
add("3. Panels decoupled: A (>=%d ancestries)=%d rows; B (>=%d cohorts)=%d rows; total leads=%d",
    MIN_ANCESTRIES, nA, MIN_COHORTS, nB, nrow(annot))
add("4. Orient to CTP-increasing allele: %s; flips applied: %d/%d; displayed allele: %s",
    ORIENT_TO_INCREASING, sum(annot$orient_sign < 0), nrow(annot),
    if (ORIENT_TO_INCREASING) "pooled CTP-increasing" else "ALT (color = ALT-allele direction)")
add("5. Label convention: celltype_rsID (celltype_chr:pos:REF:ALT if no rsID); variants w/o rsID=%d",
    annot[is.na(rsID) | rsID == "", .N])
add("6. Style: white background, no descriptive titles, panel letters A/B/C only")
add("   Panel A: Ancestry I2 column removed; x-axis 'Effect size (beta)'; legend 'Meta-analysis'")
add("   Panel B: Cohort I2 / Direction columns removed; composition + colorbar legends under B")
add("7. Panel C: 3 cohort forest plots (1 row x 3 col), rsID-only titles; LOCO in supplement")
for (e in examples)
  add("     %s | %s : I2=%.0f%%, %d cohorts", e$cell_type, e$rs, e$I2, e$k)
add("8. Heatmap symmetric limit +/- %.3f (full +/- %.3f); clipped cells: %d",
    lim_use, lim_full, nrow(clipped))
if (nrow(clipped)) for (i in seq_len(nrow(clipped)))
  add("     %s | %s : %s beta=%.3f", clipped$cell_type[i], clipped$lead_variant[i],
      clipped$stratum[i], clipped$dbeta[i])
add("9. AMR = Admixed American (Latino; AMP-AD Mayo AMR subset)")
add("10. Outputs: generalizability_figure_{ABC,AB}.{pdf,svg,png}, supplementary_LOCO_all_leads.{tsv,pdf},")
add("    panel_*.tsv, figure_effect_orientation_audit.tsv (copied to manuscript_figure/)")
writeLines(qc, file.path(OUT_DIR, "generalizability_figure_QC_report.txt"))
cat(paste(qc, collapse = "\n"), "\n\nDONE\n")
