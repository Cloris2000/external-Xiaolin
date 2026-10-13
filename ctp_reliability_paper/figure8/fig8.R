# Figure 8. Does sequencing depth explain per-donor discordance?
# Question (user, 2026-10-05): for the same donor measured in two snRNA-seq
# datasets, is cross-dataset disagreement in cell-type composition larger
# when the two datasets sequenced that donor to very different depths? This
# is SEQUENCING depth (median UMIs per nucleus), not the cortical SAMPLING
# depth already covered by Figure 3 (a laminar-position proxy; a different
# QC axis that happens to share the word "depth" -- see fig8_data.py header).
#
# Per-donor discordance (ALL cell types, user decision 2026-10-05, not just
# inhibitory as Figure 3's agree_inh): mean absolute residual CLR across every
# shared supertype, residual = that donor's CLR minus the cell type's own
# across-donor mean in its own dataset (so a systematic between-dataset
# offset in one cell type's mean proportion is not counted as discordance --
# only each donor's deviation from their own dataset's norm is compared).
#
# Panels, one per the five independent pairs used throughout this paper
# (common/reliability_stats.py): per-donor mean |residual CLR| (y) against
# that donor's depth gap, log2(median UMI in B / median UMI in A) (x). A
# Spearman rho + p is printed per pair (not assumed linear: discordance vs.
# a ratio is not expected to be a straight line, and N per pair is modest).
#
# Inputs (from fig8_data.py): data/fig8_pair_donor.tsv
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig8.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure8"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

PAIRS <- c("Green x Mathys", "Green x Multiome", "Mathys x Multiome",
          "Multiome x PsychAD", "BrainSCOPE CMC x PsychAD")
PAIR_COL <- c("Green x Mathys" = "#3C5488", "Green x Multiome" = "#8491B4",
             "Mathys x Multiome" = "#00A087", "Multiome x PsychAD" = "#B5651D",
             "BrainSCOPE CMC x PsychAD" = "#E64B35")

d <- fread(file.path(HERE, "data/fig8_pair_donor.tsv"))
d[, pair := factor(pair, levels = PAIRS)]
stopifnot(!anyNA(d$pair))

# per-pair Spearman rho of depth gap vs discordance, printed and annotated
stats <- d[, {
  s <- spear(log2_umi_ratio, mean_abs_resid_clr)
  .(n = .N, rho = unname(s$estimate), p = s$p.value)
}, by = pair]
cat("\ndepth gap (log2 UMI ratio) vs per-donor mean |residual CLR|, by pair:\n")
print(stats[order(pair)])

lab <- stats[, .(pair, lab = sprintf("n=%d\nrho=%+.2f, p=%s", n, rho, fmt_p(p)))]
lab <- merge(lab, d[, .(x = min(log2_umi_ratio) + 0.95 * (max(log2_umi_ratio) - min(log2_umi_ratio)),
                       y = max(mean_abs_resid_clr)), by = pair],
            by = "pair")

# ---- does |depth gap| predict discordance, more formally? -----------------
# (1) per-pair Wilcoxon rank-sum: large-gap tertile vs small-gap tertile
#     (unpaired groups, non-parametric -- discordance is not Normal and the
#     two tertiles are different donors, not a before/after pair).
d[, abs_gap := abs(log2_umi_ratio)]
d[, tertile := cut(abs_gap, quantile(abs_gap, c(0, 1/3, 2/3, 1)),
                   labels = c("small gap", "mid gap", "large gap"),
                   include.lowest = TRUE), by = pair]
wilcox_stats <- d[tertile != "mid gap", {
  w <- wilcox.test(mean_abs_resid_clr ~ droplevels(tertile), exact = FALSE)
  .(n_small = sum(tertile == "small gap"), n_large = sum(tertile == "large gap"),
    median_small = median(mean_abs_resid_clr[tertile == "small gap"]),
    median_large = median(mean_abs_resid_clr[tertile == "large gap"]),
    W = unname(w$statistic), p = w$p.value)
}, by = pair]
cat("\nWilcoxon rank-sum, large-gap vs small-gap tertile, by pair:\n")
print(wilcox_stats[order(pair)])

# (2) pooled test across all five pairs: does |depth gap| predict discordance
#     once each pair's own baseline level is accounted for? lme4 is broken in
#     this environment (Matrix ABI mismatch segfaults lmer -- checked
#     2026-10-07, not reinstalled per the project's no-unapproved-installs
#     rule), and five pairs is too few levels to estimate a random intercept
#     well regardless, so this uses pair as a FIXED effect (equivalent to a
#     one-way ANCOVA / de-meaning discordance and log-gap within each pair)
#     rather than a mixed model. log(abs_gap) is used, not abs_gap itself, so
#     a pair with one extreme donor cannot dominate the slope; +1 pseudo-fold
#     (0 log2-fold) keeps a perfectly-matched donor (abs_gap = 0) finite.
d[, log_abs_gap := log2(abs_gap + 1)]
m_pooled <- lm(mean_abs_resid_clr ~ log_abs_gap + pair, data = d)
m_null <- lm(mean_abs_resid_clr ~ pair, data = d)
a_test <- anova(m_null, m_pooled)
cf <- summary(m_pooled)$coefficients["log_abs_gap", ]
cat(sprintf("\npooled fixed-effects model (pair as covariate), all %d donors across 5 pairs:\n",
           nrow(d)))
cat(sprintf("  |depth gap| (log2) slope on discordance: beta=%+.4f, SE=%.4f, p=%s\n",
           cf["Estimate"], cf["Std. Error"], fmt_p(cf["Pr(>|t|)"])))
cat(sprintf("  F-test for dropping |depth gap| entirely: F=%.2f, p=%s\n",
           a_test$F[2], fmt_p(a_test$`Pr(>F)`[2])))
fwrite(wilcox_stats[order(pair)], file.path(HERE, "data/fig8_wilcoxon_tertile.tsv"), sep = "\t")

p_a <- ggplot(d, aes(log2_umi_ratio, mean_abs_resid_clr)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = pair), size = 1.1, alpha = 0.75, stroke = 0) +
  geom_smooth(aes(colour = pair), method = "loess", se = FALSE, linewidth = 0.5,
              span = 1.5, show.legend = FALSE) +
  geom_text(data = lab, aes(x, y, label = lab), hjust = 1, vjust = 1,
           size = (fontsize - 1) / ggplot2::.pt, colour = "grey30", lineheight = 0.9) +
  scale_colour_manual(values = PAIR_COL, name = NULL) +
  facet_wrap(~ pair, nrow = 1) +
  labs(x = expression(log[2]*"(median UMI, dataset B / dataset A)"),
       y = "Donor-level cross-dataset disagreement\n(mean |CLR difference|, all shared cell types)",
       subtitle = "Each point: one donor. Low = the donor's cell-type pattern looks the same in both datasets; high = the two datasets disagree") +
  theme_cns +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = fontsize - 1, colour = "grey30"))

# ---- panel b: forest plot, per-pair slope + pooled fixed-effects estimate --
# Slope of discordance on log2(|depth gap| + 1) (same x-transform as the
# pooled model above), so each pair's own effect and the pooled summary sit
# on one comparable scale. This is the direct visual counterpart of the
# Wilcoxon/pooled-model tests above: five individually noisy per-pair slopes
# (most with a CI crossing 0) versus one summary estimate, in the classic
# meta-analysis forest-plot layout.
per_pair_slope <- d[, {
  m <- lm(mean_abs_resid_clr ~ log_abs_gap)
  ci <- confint(m)["log_abs_gap", ]
  .(estimate = coef(m)["log_abs_gap"], lo = ci[1], hi = ci[2],
    p = summary(m)$coefficients["log_abs_gap", "Pr(>|t|)"], n = .N)
}, by = pair]

pooled_ci <- confint(m_pooled)["log_abs_gap", ]
pooled_row <- data.table(
  pair = "Pooled (fixed-effects)", estimate = unname(cf["Estimate"]),
  lo = pooled_ci[1], hi = pooled_ci[2], p = unname(cf["Pr(>|t|)"]), n = nrow(d))

forest <- rbind(per_pair_slope, pooled_row)
forest[, row_lab := factor(pair, levels = rev(c(PAIRS, "Pooled (fixed-effects)")))]
forest[, is_pooled := pair == "Pooled (fixed-effects)"]
forest[, lab := sprintf("n=%d   beta=%+.3f [%+.3f, %+.3f]   p=%s",
                       n, estimate, lo, hi, fmt_p(p))]

p_b <- ggplot(forest, aes(estimate, row_lab)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = pair, linewidth = is_pooled),
                height = 0, show.legend = FALSE) +
  geom_point(aes(colour = pair, shape = is_pooled, size = is_pooled)) +
  geom_text(aes(label = lab), hjust = 0, vjust = -1.1,
           size = (fontsize - 1.5) / ggplot2::.pt, colour = "grey30") +
  scale_colour_manual(values = c(PAIR_COL, "Pooled (fixed-effects)" = "black"), guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.8, "TRUE" = 3), guide = "none") +
  scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 0.9), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.38))) +
  labs(x = expression("Slope of discordance on "*log[2]*"(|depth gap| + 1)"), y = NULL) +
  theme_cns +
  theme(axis.text.y = element_text(hjust = 0),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

fig <- (p_a / p_b) +
  plot_layout(heights = c(1.5, 1.1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, "figure8_depth_vs_discordance"), height_in = 6.5, width_mm = 220)
