# Figure 8. Does dissection (cortical) depth explain per-donor discordance?
# Question (user, 2026-10-19): for the same donor measured in two snRNA-seq
# datasets, is cross-dataset disagreement in cell-type composition larger
# when the two tissue blocks sampled different cortical depths? Depth is
# DISSECTION depth inferred from excitatory supertypes by the cortdepth
# package (the PI's `cortdepth discordant` run, vintage 2026-10-13); see
# fig8_data.py. The earlier sequencing-depth (UMI) version of this figure is
# superseded; its outputs stay under figures/figure8_depth_vs_discordance.*.
#
# x  wasserstein: distance between the two samples' excitatory depth profiles,
#    in depth units (0 = pia, 1 = white-matter boundary; 0.1 = the sample's
#    nuclei shifted by a tenth of the cortex on average). Direction-free.
# y  disc_inh: the donor's mean |residual CLR_A - residual CLR_B| over
#    INHIBITORY supertypes only (CLR within the inhibitory set). Excitatory
#    types are excluded because they define x; glia are excluded because
#    their proportions (oligodendrocytes in white matter) are themselves
#    largely set by sampling depth (user decision 2026-10-23). The
#    inhibitory+glial (disc_nonexc) and all-supertype (disc_all) versions are
#    run as sensitivities and printed, not plotted.
# Donors: the PI's `pass` filter (frontal, not a duplicate, no excitatory label
# collapse, sex-concordant, identity not assumed), >= 200 excitatory nuclei on
# both sides (upstream, cortdepth).
# Panels: (a) per-donor scatter per pair with Spearman rho; (b) forest plot of
# each pair's slope of y on x with 95% CI and R^2, and the pooled fixed-effects
# slope (pair as covariate; lme4 segfaults in this env, and the eight pairs of
# common/reliability_stats.py -- five before 2026-10-29 -- are still few for a
# random intercept) with its partial R^2.
# Inputs (from fig8_data.py): data/fig8_depth_pair_donor.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig8.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure8"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

# the eight pairs of common/reliability_stats.py (2026-10-29; were five)
PAIRS <- c("Green x Mathys", "Green x Multiome", "Mathys x Multiome",
          "Multiome x PsychAD ROSMAP", "BrainSCOPE CMC x PsychAD MSBB",
          "CUIMC2 x Multiome", "CUIMC2 x Mathys", "Catching HBCC x PsychAD HBCC")
# Figure 3's pair palette (cnsplots "Cell"): one colour per pair across the paper
PAIR_COL <- c("Green x Mathys" = "#2F7E8F", "Green x Multiome" = "#E1A22E",
             "Mathys x Multiome" = "#8B6FA8", "Multiome x PsychAD ROSMAP" = "#C84C3A",
             "BrainSCOPE CMC x PsychAD MSBB" = "#5F9862", "CUIMC2 x Multiome" = "#4E5A8A",
             "CUIMC2 x Mathys" = "#D07A6A", "Catching HBCC x PsychAD HBCC" = "#6B6B6B")

# Displayed dataset names are the canonical_name of common/cohort_names.tsv,
# looked up by the wp3 source_key; the pair keys above stay as data keys only.
# CUIMC2 and Catching_2026_HBCC have no row there yet (naming awaits PI review),
# so they show their registry dataset_id.
cohorts <- fread(file.path(dirname(HERE), "common", "cohort_names.tsv"))
canon <- function(key) {
  nm <- cohorts$canonical_name[cohorts$source_key == key]
  if (!length(nm)) return(key)
  stopifnot(length(nm) == 1)
  nm
}
SIDE_KEY <- c("Green" = "Green_DFC_person", "Mathys" = "Mathys_DFC_person",
              "Multiome" = "Multiome2025_DLPFC_DFC_person",
              "PsychAD ROSMAP" = "PsychAD_ROSMAP_DFC", "PsychAD MSBB" = "PsychAD_MSBB_DFC",
              "PsychAD HBCC" = "PsychAD_HBCC_DFC", "BrainSCOPE CMC" = "BrainSCOPE_CMC_DFC_person",
              "CUIMC2" = "CUIMC2", "Catching HBCC" = "Catching_2026_HBCC")
PAIR_LAB <- setNames(vapply(strsplit(PAIRS, " x "), function(s)
  paste0(canon(SIDE_KEY[[s[1]]]), "\n× ", canon(SIDE_KEY[[s[2]]])), ""), PAIRS)
# facet strips are ~45 mm wide: break a long canonical name after an underscore
# (strip labels only; panel b's row labels have room for the full line)
wrap_name <- function(nm, width = 24) {
  if (nchar(nm) <= width) return(nm)
  cut <- max(gregexpr("_", substr(nm, 1, width))[[1]])
  if (cut < 1) return(nm)
  paste0(substr(nm, 1, cut), "\n", substr(nm, cut + 1, nchar(nm)))
}
STRIP_LAB <- setNames(vapply(strsplit(PAIRS, " x "), function(s)
  paste0(wrap_name(canon(SIDE_KEY[[s[1]]])), "\n× ", wrap_name(canon(SIDE_KEY[[s[2]]]))), ""), PAIRS)
POOLED <- "Pooled (fixed-effects)"
Y_LAB <- "Donor-level cross-dataset disagreement\n(mean |CLR difference|, inhibitory supertypes)"
X_LAB <- "Dissection-depth discordance between the two samples\n(Wasserstein distance of excitatory depth profiles, cortex units)"

d0 <- fread(file.path(HERE, "data/fig8_depth_pair_donor.tsv"))
d0[, pair := factor(pair, levels = PAIRS)]
stopifnot(!anyNA(d0$pair))
cat(sprintf("rows %d, pass %d\n", nrow(d0), sum(d0$pass)))
d <- d0[pass == TRUE]
d[, abs_shift := abs(mean_shift)]

# ---- pooled fixed-effects test, main and sensitivity --------------------
pooled <- function(y, x) {
  f <- lm(reformulate(c(x, "pair"), y), data = d)
  n <- lm(reformulate("pair", y), data = d)
  cf <- summary(f)$coefficients[x, ]
  ci <- confint(f)[x, ]
  data.table(y = y, x = x, n = nrow(d), beta = unname(cf["Estimate"]), lo = ci[1], hi = ci[2],
             p = unname(cf["Pr(>|t|)"]),
             partial_r2 = 1 - sum(residuals(f)^2) / sum(residuals(n)^2))
}
sens <- rbindlist(lapply(c("disc_inh", "disc_nonexc", "disc_all"), function(y)
  rbindlist(lapply(c("wasserstein", "abs_shift", "layer_tv"), function(x) pooled(y, x)))))
cat("\npooled fixed-effects (pair as covariate), all passing donors:\n")
print(sens[, .(y, x, n, beta = round(beta, 3), lo = round(lo, 3), hi = round(hi, 3),
               p = fmt_p(p), partial_r2 = sprintf("%.1f%%", 100 * partial_r2))])
fwrite(sens, file.path(HERE, "data/fig8_depth_pooled_stats.tsv"), sep = "\t")
main <- sens[y == "disc_inh" & x == "wasserstein"]

# ---- per-pair Spearman, for panel a ---------------------------------------
stats <- d[, {
  s <- spear(wasserstein, disc_inh)
  .(n = .N, rho = unname(s$estimate), p = s$p.value)
}, by = pair]
cat("\nper pair, Spearman of depth discordance vs composition discordance (inhibitory):\n")
print(stats[order(pair)])
lab <- stats[, .(pair, lab = sprintf("n=%d\nrho=%+.2f, p=%s", n, rho, fmt_p(p)))]
lab[, `:=`(x = min(d$wasserstein), y = max(d$disc_inh))]

p_a <- ggplot(d, aes(wasserstein, disc_inh)) +
  geom_point(aes(colour = pair), size = 1.1, alpha = 0.75, stroke = 0) +
  geom_smooth(aes(colour = pair), method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = lab, aes(x, y, label = lab), hjust = 0, vjust = 1,
           size = (fontsize - 1) / ggplot2::.pt, colour = "grey30", lineheight = 0.9) +
  scale_colour_manual(values = PAIR_COL, name = NULL) +
  facet_wrap(~ pair, nrow = 2, labeller = as_labeller(STRIP_LAB)) +
  labs(x = X_LAB, y = Y_LAB,
       subtitle = "One point per donor profiled by both datasets; y uses inhibitory supertypes only (excitatory types define x; glial shares track depth)") +
  theme_cns +
  theme(legend.position = "none",
        strip.text = element_text(size = fontsize - 1.5, lineheight = 0.95),
        plot.subtitle = element_text(size = fontsize - 1, colour = "grey30"))

# ---- panel b: forest plot, per-pair slope + pooled ------------------------
per_pair <- d[, {
  m <- lm(disc_inh ~ wasserstein)
  ci <- confint(m)["wasserstein", ]
  .(estimate = coef(m)["wasserstein"], lo = ci[1], hi = ci[2],
    p = summary(m)$coefficients["wasserstein", "Pr(>|t|)"], n = .N, r2 = summary(m)$r.squared)
}, by = pair]
pooled_row <- data.table(pair = POOLED, estimate = main$beta, lo = main$lo,
                         hi = main$hi, p = main$p, n = main$n, r2 = main$partial_r2)
forest <- rbind(per_pair, pooled_row)
forest[, is_pooled := pair == POOLED]
forest[, row_lab := factor(fifelse(is_pooled, POOLED, PAIR_LAB[as.character(pair)]),
                           levels = rev(c(unname(PAIR_LAB), POOLED)))]
forest[, lab := sprintf("n=%d   beta=%+.2f [%+.2f, %+.2f]   p=%s   R²=%.1f%%",
                       n, estimate, lo, hi, fmt_p(p), 100 * r2)]

p_b <- ggplot(forest, aes(estimate, row_lab)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = pair, linewidth = is_pooled),
               width = 0, orientation = "y", show.legend = FALSE) +
  geom_point(aes(colour = pair, shape = is_pooled, size = is_pooled)) +
  geom_text(aes(label = lab), hjust = 0, vjust = -1.1,
           size = (fontsize - 1.5) / ggplot2::.pt, colour = "grey30") +
  scale_colour_manual(values = c(PAIR_COL, setNames("black", POOLED)), guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.8, "TRUE" = 3), guide = "none") +
  scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 0.9), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.75))) +
  labs(x = "Slope of composition disagreement per unit of depth discordance (Wasserstein)",
       y = NULL) +
  theme_cns +
  theme(axis.text.y = element_text(hjust = 0, lineheight = 0.95),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

# free(): panel a's facets must not be squeezed by panel b's two-line row labels
fig <- (free(p_a, type = "panel") / free(p_b, type = "panel")) +
  plot_layout(heights = c(1.5, 1.3)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, "figure8_dissection_depth_vs_discordance"),
         height_in = 6.8, width_mm = 240)
