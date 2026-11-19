# Figure 10 supplement (DRAFT): Sst_25 in SCZ and AD.
#   a,b  forest plots: per-cohort case-vs-control beta +- 95% CI (Fig 10 primary model:
#        z-scored outcome ~ case + age + sex + PMI; SCZ adults 18-69, AD all ages), bulk
#        (production CelMod, out-of-sample donors) and snRNA (within-compartment CLR), each
#        with its DerSimonian-Laird random-effects meta-analysis (diamond).
#   c,d  snRNA Sst_25 proportion, % of all neuronal nuclei, per donor, control vs case,
#        one facet per cohort (same donors as the snRNA rows of a,b). Facet label carries the
#        model p from a,b.
#   c,d  SAME-DONOR bulk vs snRNA (fig10_sst25_matched.py): people in both a bulk and an
#        snRNA cohort, Sst_25 z-scored CLR per modality, model refitted on those donors.
#   The c/e snRNA %-of-neurons boxes are written separately (figure10_supp_sst25_pct_draft).
# Inputs: data/fig10_effects_by_cohort.tsv, fig10_effects_meta.tsv, fig10_sst25_props.tsv,
#         fig10_sst25_bulk.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig10_sst25.R
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(HERE, "..", "common", "theme_cns.R"))

TYPE <- "Sst_25"
DIS <- c(SCZ = "Schizophrenia", AD = "Alzheimer's disease")
MOD_LAB <- c(bulk = "Bulk RNA-seq (deconvolved)", snRNA = "snRNA-seq")
MOD_COL <- c(bulk = "#E69F00", snRNA = "#0072B2")   # Okabe-Ito orange / blue (colour-blind safe)
GRP_COL <- c(Control = "#B7B7B7", SCZ = "#E64B35", AD = "#3C5488")
Z <- qnorm(0.975)

pc <- fread(file.path(HERE, "data", "fig10_effects_by_cohort.tsv"))[
  cell_type == TYPE & model == "primary" & disease %in% names(DIS)]
pm <- fread(file.path(HERE, "data", "fig10_effects_meta.tsv"))[
  cell_type == TYPE & model == "primary" & disease %in% names(DIS)]
pr <- fread(file.path(HERE, "data", "fig10_sst25_props.tsv"))[disease %in% names(DIS)]
pb <- fread(file.path(HERE, "data", "fig10_sst25_bulk.tsv"))[disease %in% names(DIS)]
stopifnot(nrow(pc) > 0, nrow(pm) == 4)

fmtp <- function(p) ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p))

# one forest per disease x modality; stacked with patchwork (heights = number of rows)
forest_one <- function(d, mod, show_x) {
  # N per cohort, I² and p values are given in the legend text, not on the plot
  c1 <- pc[disease == d & modality == mod][, .(lab = cohort_label(cohort), beta, se, p, meta = FALSE)]
  m1 <- pm[disease == d & modality == mod][, .(lab = "Meta-analysis", beta, se, p, meta = TRUE)]
  x <- rbind(c1, m1)                         # cohorts by decreasing N (input order), meta last
  x[, y := rev(seq_len(.N))]
  x[, lo := beta - Z * se][, hi := beta + Z * se]
  dia <- x[meta == TRUE, .(px = c(lo, beta, hi, beta), py = y + c(0, 0.32, 0, -0.32))]
  g <- ggplot(x, aes(y = y)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50", linetype = "dashed") +
    geom_errorbarh(data = x[meta == FALSE], aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
    geom_point(data = x[meta == FALSE], aes(x = beta, size = 1 / se^2), shape = 15) +
    geom_polygon(data = dia, aes(x = px, y = py), fill = MOD_COL[[mod]], colour = "black",
                 linewidth = 0.3, inherit.aes = FALSE) +
    scale_y_continuous(breaks = x$y, labels = x$lab, expand = expansion(add = 0.6)) +
    scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5), labels = c("-1", "-0.5", "0", "0.5")) +
    scale_size_continuous(range = c(0.8, 2.4), guide = "none") +
    coord_cartesian(xlim = XLIM, clip = "off") +
    labs(x = if (show_x) "Sst_25 abundance change\n(β ± 95% CI, SD units)" else NULL,
         y = NULL, subtitle = MOD_LAB[[mod]]) +
    theme_cns +
    theme(plot.subtitle = element_text(size = fontsize, face = "bold"), plot.title.position = "plot",
          plot.margin = margin(2, 8, 2, 4))
  if (!show_x) g <- g + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                              axis.line.x = element_blank())
  g
}
XLIM <- range(c(pc$beta - Z * pc$se, pc$beta + Z * pc$se, 0))

forest <- function(d) {
  nb <- pc[disease == d & modality == "bulk", .N] + 1
  ns <- pc[disease == d & modality == "snRNA", .N] + 1
  (forest_one(d, "bulk", FALSE) + ggtitle(DIS[[d]]) +
     theme(plot.title = element_text(size = fontsize, face = "bold", margin = margin(l = 14, b = 2)))) /
    forest_one(d, "snRNA", TRUE) + plot_layout(heights = c(nb, ns))
}

boxes <- function(d, mod = "snRNA") {
  x <- if (mod == "snRNA") pr[disease == d][, yval := pct_of_neurons] else pb[disease == d][, yval := z]
  x[, group := factor(ifelse(case == 1, d, "Control"), levels = c("Control", d))]
  ps <- pc[disease == d & modality == mod, .(cohort, p, n_tot = n_case + n_ctrl)]
  x <- merge(x, ps, by = "cohort")
  x[, facet := paste0(gsub("_", "\n", cohort_label(cohort)), "\np = ", fmtp(p))]
  x[, facet := factor(facet, levels = unique(x[order(-n_tot)]$facet))]
  # display only: points above Q3 + 5 IQR of their cohort drawn as triangles at that cap
  # (e.g. a PsychAD_MSBB donor with 2 neuronal nuclei, 1 of them Sst_25 = 50%); box
  # quartiles and whiskers are unaffected; all donors stay in the models of a, b
  x[, cap := quantile(yval, 0.75) + 5 * IQR(yval), by = cohort]
  x[, capped := yval > cap][, yplot := pmin(yval, cap)]
  ggplot(x, aes(x = group, y = yplot)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.65, linewidth = 0.3, alpha = 0.7) +
    geom_jitter(data = x[capped == FALSE], width = 0.18, height = 0, size = 0.3, alpha = 0.45) +
    geom_point(data = x[capped == TRUE], shape = 2, size = 1.2, colour = "black") +
    { if (mod == "bulk") geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = "dashed") } +
    facet_wrap(~facet, nrow = 1, scales = if (mod == "snRNA") "free_y" else "fixed") +
    scale_fill_manual(values = GRP_COL, name = NULL) +
    labs(x = NULL, y = if (mod == "snRNA") "snRNA-seq Sst_25\n(% of neuronal nuclei)"
                       else "Bulk deconvolved Sst_25\n(CLR, z-scored)") +
    theme_cns +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text = element_text(size = 6.5), panel.spacing.x = unit(5, "pt"),
          legend.position = "right")
}

# ---- c, d: bulk and snRNA on ONE scale, cohorts grouped by brain bank ------------------
# Both modalities shown as the z-scored within-compartment CLR of Sst_25 (z within cohort on
# the analysis donors), i.e. exactly the outcome whose case-control beta is in a, b. Banks
# are grouped only where the two modalities share donors (fig10_donors.tsv overlap:
# HBCC 151; Mount Sinai CMC-MSSM x PsychAD MSBB 72, MSBB_FP x PsychAD MSBB 206, Diverse MSSM
# x PsychAD MSBB 31; Rush ROSMAP x CUIMC2 87, Diverse Rush x PsychAD ROSMAP 41). Green,
# Mathys and Multiome are ROSMAP people too; their donors trained CelMod and so are absent
# from bulk. Every other cohort stands alone.
BANK <- c(NIMH_HBCC = "HBCC", PsychAD_HBCC = "HBCC",
          CommonMind_MSSM = "Mount Sinai", MSBB_FP = "Mount Sinai",
          Diverse_MSSM_DLPFC = "Mount Sinai", PsychAD_MSBB = "Mount Sinai",
          ROSMAP = "Rush (ROSMAP)", Diverse_Rush2_DLPFC = "Rush (ROSMAP)", Green = "Rush (ROSMAP)",
          Mathys = "Rush (ROSMAP)", Multiome2025_DLPFC = "Rush (ROSMAP)", CUIMC2 = "Rush (ROSMAP)",
          PsychAD_ROSMAP = "Rush (ROSMAP)")
YCAP <- 3.5
compare <- function(d) {
  x <- rbind(pb[disease == d, .(cohort, modality = "bulk", case, z)],
             pr[disease == d, .(cohort, modality = "snRNA", case, z)])
  x[, group := factor(ifelse(case == 1, d, "Control"), levels = c("Control", d))]
  x[, bank := ifelse(cohort %in% names(BANK), BANK[cohort], gsub("_", " ", cohort_label(cohort)))]
  ps <- pc[disease == d, .(cohort, modality, p, n_tot = n_case + n_ctrl)]
  x <- merge(x, ps, by = c("cohort", "modality"))
  # bank order: shared banks (both modalities) by N, then bulk-only, then snRNA-only
  bk <- x[, .(both = uniqueN(modality) == 2, mod1 = modality[1], n = .N), by = bank]
  bk <- bk[order(-both, mod1, -n)]
  x[, bank := factor(bank, levels = bk$bank)]
  # within a bank: bulk cohorts first, then snRNA, each by N
  co <- unique(x[order(bank, modality, -n_tot), .(cohort, modality)])
  co[, key := paste(modality, cohort)]
  co[, lab := paste0(cohort_label(cohort), "\n", ifelse(modality == "bulk", "bulk", "snRNA"))]
  x[, key := factor(paste(modality, cohort), levels = co$key)]
  x[, capped := abs(z) > YCAP][, yplot := pmax(pmin(z, YCAP), -YCAP)]
  pl <- unique(x[, .(bank, key, modality, p)])
  ggplot(x, aes(x = key, y = yplot)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = "dashed") +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.75,
                 linewidth = 0.3, alpha = 0.75, position = position_dodge(0.85)) +
    geom_point(data = x[capped == FALSE], aes(fill = group),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85, seed = 1),
               size = 0.2, alpha = 0.35, show.legend = FALSE) +
    geom_point(data = x[capped == TRUE], aes(fill = group), shape = 2, size = 0.9,
               position = position_dodge(0.85), show.legend = FALSE) +
    geom_text(data = pl, aes(y = YCAP + 0.55, label = paste0("p=", fmtp(p))), size = 2.1) +
    facet_grid(. ~ bank, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 10)) +
    scale_x_discrete(labels = setNames(co$lab, co$key)) +
    scale_fill_manual(values = GRP_COL, name = NULL) +
    scale_y_continuous(breaks = seq(-3, 3, 1)) +
    coord_cartesian(ylim = c(-YCAP, YCAP + 0.7), clip = "off") +
    labs(x = NULL, y = "Sst_25, z-scored CLR\n(within cohort)", title = DIS[[d]]) +
    theme_cns +
    theme(axis.text.x = element_text(size = 5.5, angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(size = 6.5, face = "bold"), strip.clip = "off",
          panel.spacing.x = unit(6, "pt"),
          plot.title = element_text(size = fontsize, face = "bold"),
          legend.position = "right")
}

# ---- c, d: SAME-DONOR comparison (fig10_sst25_matched.py) ------------------------------
# People in both a bulk and an snRNA analysis set; each modality's Sst_25 CLR z-scored over
# the matched donors; p = case term of the Fig 10 model refitted on those donors; rho = per-donor
# Spearman rho between bulk and snRNA Sst_25 (rank-based, so identical on z or raw CLR).
# Pairs with < 10 matched cases or controls omitted.
md <- fread(file.path(HERE, "data", "fig10_sst25_matched_donors.tsv"), colClasses = list(character = "person_id"))
me <- fread(file.path(HERE, "data", "fig10_sst25_matched_effects.tsv"))
rs <- dcast(md, disease + bulk_cohort + snrna_cohort + person_id ~ modality, value.var = "z")
stopifnot(!anyNA(rs$bulk), !anyNA(rs$snRNA))
rs <- rs[, .(rho = cor(bulk, snRNA, method = "spearman"), n_rho = .N), by = .(disease, bulk_cohort, snrna_cohort)]
me <- merge(me, rs, by = c("disease", "bulk_cohort", "snrna_cohort"))
stopifnot(me$n_rho == me$n)
matched <- function(d) {
  x <- md[disease == d]; e <- me[disease == d]
  key <- function(b, s) paste(b, s)
  lab <- e[modality == "bulk"][order(-n)][, .(k = key(bulk_cohort, snrna_cohort),
    # n and Spearman rho per pair are given in the legend text (printed below), not on the plot
    lab = paste0(cohort_label(bulk_cohort), " bulk\n× ", cohort_label(snrna_cohort), " snRNA"))]
  x[, pair := factor(key(bulk_cohort, snrna_cohort), levels = lab$k, labels = lab$lab)]
  e[, pair := factor(key(bulk_cohort, snrna_cohort), levels = lab$k, labels = lab$lab)]
  # x = diagnosis group, fill = modality: each modality's control -> case shift side by side.
  # Both modalities are z-scored over the same matched donors (mean 0 each), so what is
  # compared is the SHIFT between groups, not the absolute level.
  x[, group := factor(ifelse(case == 1, d, "Control"), levels = c("Control", d))]
  x[, mod := factor(modality, levels = c("bulk", "snRNA"))]
  x[, capped := abs(z) > YCAP][, yplot := pmax(pmin(z, YCAP), -YCAP)]
  pl <- dcast(e, pair ~ modality, value.var = "p")
  pl[, lab := paste0("Bulk p = ", fmtp(bulk), "\nsnRNA p = ", fmtp(snRNA))]
  ggplot(x, aes(x = group, y = yplot)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = "dashed") +
    geom_boxplot(aes(fill = mod), outlier.shape = NA, width = 0.7, linewidth = 0.3, alpha = 0.7,
                 position = position_dodge(0.8)) +
    geom_point(data = x[capped == FALSE], aes(fill = mod), show.legend = FALSE,
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1),
               size = 0.3, alpha = 0.35) +
    geom_point(data = x[capped == TRUE], aes(fill = mod), shape = 2, size = 0.9,
               position = position_dodge(0.8), show.legend = FALSE) +
    geom_text(data = pl, aes(x = 1.5, y = YCAP + 0.75, label = lab), size = 2.1, lineheight = 0.9,
              inherit.aes = FALSE) +
    facet_wrap(~pair, nrow = 1) +
    scale_fill_manual(values = MOD_COL, labels = c(bulk = "Bulk RNA-seq (deconvolved)", snRNA = "snRNA-seq"),
                      name = NULL) +
    scale_y_continuous(breaks = seq(-3, 3, 1)) +
    coord_cartesian(ylim = c(-YCAP, YCAP + 1.1), clip = "off") +
    labs(x = NULL, y = "Sst_25 (z-scored CLR)", title = DIS[[d]]) +
    theme_cns +
    theme(strip.text = element_text(size = 6), strip.clip = "off",
          plot.title = element_text(size = fontsize, face = "bold"),
          legend.position = "bottom")
}

fig <- wrap_elements(full = forest("SCZ")) + wrap_elements(full = forest("AD")) +
  wrap_elements(full = matched("SCZ")) + wrap_elements(full = matched("AD")) +
  plot_layout(design = "AB\nCD", heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "a")
save_fig(fig, file.path(HERE, "figures", "figure10_supp_sst25_draft"), height_in = 6.3)

# numbers removed from the plot, for the legend text
cat("\n-- per cohort (a, b) --\n")
print(pc[, .(disease, modality, cohort = cohort_label(cohort), n_case, n_ctrl, beta = round(beta, 2), p = signif(p, 2))])
cat("\n-- meta (a, b) --\n")
print(pm[, .(disease, modality, n_cohorts, n_case, n_ctrl, beta = round(beta, 2), se = round(se, 3),
             p = signif(p, 2), I2 = round(I2))])
cat("\n-- same-donor pairs (c, d) --\n")
print(me[modality == "bulk", .(disease, bulk_cohort, snrna_cohort, n, n_case, n_ctrl, rho = round(rho, 2))])

# snRNA-only view in the unit asked for earlier (% of neuronal nuclei), kept as its own file
pct <- wrap_elements(full = boxes("SCZ", "snRNA")) / wrap_elements(full = boxes("AD", "snRNA")) +
  plot_annotation(tag_levels = "a")
save_fig(pct, file.path(HERE, "figures", "figure10_supp_sst25_pct_draft"), height_in = 4.6)
