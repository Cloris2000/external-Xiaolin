# Figure 8 (simple layer-marker version; the Figure 8 design since 2026-11-03,
# replacing the cortdepth/Xenium Wasserstein axis): sampled-depth mismatch from
# two subclass layer markers and no reference (see fig8_simple_check_data.py).
#   x  |upper_frac_A - upper_frac_B|, upper_frac = L2/3 IT / (L2/3 IT + L6 IT)
#      among each sample's excitatory nuclei
#   y  chosen by the first argument (2026-11-03, both rendered for comparison):
#        disc  disc_inh, mean |delta residual CLR| over inhibitory supertypes
#              (higher = more discordant) -> figure8_simple_depth_check.*
#        corr  agree_inh, Figure 3c's per-donor Pearson r of the two inhibitory
#              profiles (higher = more concordant) -> figure8_simple_depth_check_corr.*
#        tv    tv_inh, % of the donor's inhibitory neurons that differ between the
#              two datasets (the Figure 8 y, 2026-11-03) -> figure8_simple_depth_check_tv.*
#        icc   icc_inh, per-donor symmetric ICC of the same profiles (also scores
#              agreement in size, like the paper's metric) -> figure8_simple_depth_check_icc.*
#   Pairs: the eight of common/reliability_stats.py (2026-10-29; were five).
#   (a) per-pair scatter, same donors as the main figure (pass filter, plus
#       >= 50 marker nuclei on both sides so the fraction is not noise)
#   (b) the simple index against the Wasserstein score, pooled
# Cohort names: display_label of common/cohort_names.tsv, looked up by the
# wp3 source_key. Pair colours: the Figure 3 pair palette (cnsplots "Cell"),
# one colour per pair across the paper.
# Inputs: data/fig8_simple_depth.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig8_simple_check.R [tv|disc|corr|icc]

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure8"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")
MIN_MARKER <- 50
Y <- match.arg(commandArgs(trailingOnly = TRUE)[1], c("tv", "disc", "corr", "icc"))
YCOL <- c(tv = "tv_inh", disc = "disc_inh", corr = "agree_inh", icc = "icc_inh")[[Y]]
Y_LAB <- c(tv = "Inhibitory neurons that differ\nbetween the two datasets (%)",
           disc = "Discordance of inhibitory cell proportions\n(mean |ΔCLR|)",
           corr = "Agreement of inhibitory cell proportions\n(per-donor r between datasets)",
           icc = "Agreement of inhibitory cell proportions\n(per-donor ICC between datasets)")[[Y]]
OUT_STEM <- c(tv = "figure8_simple_depth_check_tv", disc = "figure8_simple_depth_check", corr = "figure8_simple_depth_check_corr",
              icc = "figure8_simple_depth_check_icc")[[Y]]

PAIRS <- c("Green x Mathys", "Green x Multiome", "Mathys x Multiome",
          "Multiome x PsychAD ROSMAP", "BrainSCOPE CMC x PsychAD MSBB",
          "CUIMC2 x Multiome", "CUIMC2 x Mathys", "Catching HBCC x PsychAD HBCC")
PAIR_COL <- c("Green x Mathys" = "#2F7E8F", "Green x Multiome" = "#E1A22E",
             "Mathys x Multiome" = "#8B6FA8", "Multiome x PsychAD ROSMAP" = "#C84C3A",
             "BrainSCOPE CMC x PsychAD MSBB" = "#5F9862", "CUIMC2 x Multiome" = "#4E5A8A",
             "CUIMC2 x Mathys" = "#D07A6A", "Catching HBCC x PsychAD HBCC" = "#6B6B6B")
cohorts <- fread(file.path(dirname(HERE), "common", "cohort_names.tsv"))
display <- function(key) {
  nm <- cohorts$display_label[cohorts$source_key == key]
  if (!length(nm)) return(key)   # CUIMC2, Catching_2026_HBCC: no row yet (naming awaits PI review)
  stopifnot(length(nm) == 1)
  nm
}
SIDE_KEY <- c("Green" = "Green_DFC_person", "Mathys" = "Mathys_DFC_person",
              "Multiome" = "Multiome2025_DLPFC_DFC_person",
              "PsychAD ROSMAP" = "PsychAD_ROSMAP_DFC", "PsychAD MSBB" = "PsychAD_MSBB_DFC",
              "PsychAD HBCC" = "PsychAD_HBCC_DFC", "BrainSCOPE CMC" = "BrainSCOPE_CMC_DFC_person",
              "CUIMC2" = "CUIMC2", "Catching HBCC" = "Catching_2026_HBCC")
PAIR_LAB <- setNames(vapply(strsplit(PAIRS, " x "), function(s)
  paste0(display(SIDE_KEY[[s[1]]]), " × ", display(SIDE_KEY[[s[2]]])), ""), PAIRS)
STRIP_LAB <- setNames(sub(" × ", "\n× ", PAIR_LAB), PAIRS)

d <- fread(file.path(HERE, "data/fig8_simple_depth.tsv"))
d[, pair := factor(pair, levels = PAIRS)]
d <- d[pass == TRUE & pmin(n_l23_a + n_l6it_a, n_l23_b + n_l6it_b) >= MIN_MARKER]
d[, yv := get(YCOL)]
d <- d[is.finite(yv)]
cat(sprintf("donors used: %d\n", nrow(d)))

stats <- d[, {
  s <- spear(abs_d_upper, yv)
  .(n = .N, rho = unname(s$estimate), p = s$p.value)
}, by = pair]
cat(sprintf("\nper pair, Spearman of |d upper_frac| vs %s:\n", YCOL))
print(stats[order(pair)])
f <- lm(yv ~ abs_d_upper + pair, d); n0 <- lm(yv ~ pair, d)
cf <- summary(f)$coefficients["abs_d_upper", ]
cat(sprintf("\npooled: beta=%.3f, p=%s, partial R2=%.1f%%\n", cf["Estimate"], fmt_p(cf["Pr(>|t|)"]),
            100 * (1 - sum(residuals(f)^2) / sum(residuals(n0)^2))))
sw <- spear(d$abs_d_upper, d$wasserstein)
cat(sprintf("simple index vs Wasserstein: Spearman rho=%.2f, p=%s\n", sw$estimate, fmt_p(sw$p.value)))

# P as "2 x 10^-9" in plotmath, "0.48" when >= 0.001
p_math <- function(p) vapply(p, function(q) {
  if (q < 1e-3) {
    e <- floor(log10(q)); m <- round(q / 10^e)
    if (m == 10) { m <- 1; e <- e + 1 }   # 9.7e-10 -> 1 x 10^-9, not 10 x 10^-10
    sprintf("%.0f %%*%% 10^%d", m, e)
  } else sprintf("%.2g", q)
}, "")
lab <- stats[, .(pair, lab_n = sprintf("italic(n) == %d", n),
                 lab_r = sprintf("rho == '%.2f' * ',' ~ italic(P) == %s", rho, p_math(p)))]
lab[, `:=`(x = min(d$abs_d_upper), y = Inf)]   # top of the panel, in the headroom above the data

p_a <- ggplot(d, aes(abs_d_upper, yv)) +
  geom_point(aes(colour = pair), size = 0.9, alpha = 0.7, stroke = 0) +
  geom_smooth(aes(colour = pair), method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.45, show.legend = FALSE) +
  geom_text(data = lab, aes(x, y, label = lab_n), parse = TRUE, hjust = 0, vjust = 1.3,
           size = (fontsize - 1) / ggplot2::.pt, colour = "grey30") +
  geom_text(data = lab, aes(x, y, label = lab_r), parse = TRUE, hjust = 0, vjust = 2.6,
           size = (fontsize - 1) / ggplot2::.pt, colour = "grey30") +
  scale_colour_manual(values = PAIR_COL, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.4, 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.32))) +   # headroom for the n / rho text
  facet_wrap(~ pair, nrow = 2, labeller = as_labeller(STRIP_LAB)) +
  labs(x = "Difference in sampled depth between the two samples\n(|Δ upper-layer share of excitatory neurons|, L2/3 IT / (L2/3 IT + L6 IT))",
       y = Y_LAB) +
  theme_cns +
  theme(legend.position = "none", panel.spacing.x = unit(4, "pt"))

# ---- (b) forest plot: per-pair slope + pooled fixed-effects (meta-analytic) ----
per_pair <- d[, {
  m <- lm(yv ~ abs_d_upper)
  ci <- confint(m)["abs_d_upper", ]
  .(estimate = coef(m)["abs_d_upper"], lo = ci[1], hi = ci[2],
    p = summary(m)$coefficients["abs_d_upper", "Pr(>|t|)"], n = .N, r2 = summary(m)$r.squared)
}, by = pair]
POOLED <- "Pooled (fixed effects)"
ci_p <- confint(f)["abs_d_upper", ]
pooled_row <- data.table(pair = POOLED, estimate = unname(cf["Estimate"]), lo = ci_p[1], hi = ci_p[2],
                         p = unname(cf["Pr(>|t|)"]), n = nrow(d),
                         r2 = 1 - sum(residuals(f)^2) / sum(residuals(n0)^2))
forest <- rbind(per_pair, pooled_row)
forest[, is_pooled := pair == POOLED]
forest[, row_lab := factor(fifelse(is_pooled, POOLED, PAIR_LAB[as.character(pair)]),
                           levels = rev(c(unname(PAIR_LAB), POOLED)))]
forest[, lab := sprintf("%.2f [%.2f, %.2f]", estimate, lo, hi)]

p_b <- ggplot(forest, aes(estimate, row_lab)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = pair, linewidth = is_pooled),
               width = 0, orientation = "y", show.legend = FALSE) +
  geom_point(aes(colour = pair, shape = is_pooled, size = is_pooled)) +
  scale_colour_manual(values = c(PAIR_COL, setNames("black", POOLED)), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 2.8), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 0.8), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.9))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
  labs(x = sprintf("Effect of sampled-depth difference on inhibitory composition %s (β, 95%% CI)",
                   c(tv = "difference", disc = "discordance", corr = "agreement", icc = "agreement")[[Y]]),
       y = NULL) +
  theme_cns +
  theme(axis.text.y = element_text(hjust = 0),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ---- statistics table to the right of the forest plot (same row order) -------
tab <- forest[, .(row_lab, beta_ci = lab, P = fmt_p(p), R2 = sprintf("%.0f%%", 100 * r2))]
cols <- c(beta_ci = "β [95% CI]", P = "P", R2 = "R²")
# column centres chosen so the whitespace between columns is equal: the beta column is
# ~3x wider than P and R2, so equal centre spacing would leave unequal gaps
xs <- c(1, 2.0, 2.7)
tab_long <- rbindlist(lapply(seq_along(cols), function(j)
  data.table(row_lab = tab$row_lab, col = names(cols)[j], x = xs[j], value = tab[[names(cols)[j]]])))
p_tab <- ggplot(tab_long, aes(x, row_lab, label = value)) +
  geom_text(size = (fontsize - 1.5) / ggplot2::.pt, colour = "grey20", hjust = 0.5) +
  annotate("text", x = xs, y = length(levels(forest$row_lab)) + 0.9,
           label = unname(cols), size = (fontsize - 1) / ggplot2::.pt, fontface = "bold", hjust = 0.5) +
  scale_x_continuous(limits = c(0.25, 3.25), expand = c(0, 0),
                     breaks = xs, labels = rep(" ", length(cols))) +
  scale_y_discrete(limits = levels(forest$row_lab), expand = expansion(add = c(0.6, 1.0))) +
  labs(tag = " ", x = " ") +
  theme_void(base_size = fontsize) +
  theme(plot.margin = margin(2, 4, 2, 4),
        axis.title.x = element_text(size = fontsize),
        axis.text.x = element_text(size = fontsize))

fwd <- p_b + labs(tag = "b") + theme(plot.margin = margin(2, 0, 2, 4))
p_a <- p_a + labs(tag = "a")
p_b <- (fwd | p_tab) + plot_layout(widths = c(2.2, 1.9))   # aligned panels keep table rows on forest rows

fig <- (free(p_a, type = "panel") / p_b) +
  plot_layout(heights = c(1.15, 0.95)) +
  plot_annotation() &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, OUT_STEM), height_in = 7.2, width_mm = 183)
