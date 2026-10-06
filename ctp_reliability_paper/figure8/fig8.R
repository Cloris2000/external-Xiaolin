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
lab <- merge(lab, d[, .(x = max(log2_umi_ratio), y = max(mean_abs_resid_clr)), by = pair],
            by = "pair")

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
       y = "Donor mean |residual CLR|\n(all shared supertypes)") +
  theme_cns +
  theme(legend.position = "none")

# per-pair summary: median |residual CLR| in donors with the biggest vs.
# smallest depth gap (top/bottom tertile by |log2 UMI ratio|), so the panel
# above's per-donor scatter has a simple paired comparison alongside it
d[, abs_gap := abs(log2_umi_ratio)]
d[, tertile := cut(abs_gap, quantile(abs_gap, c(0, 1/3, 2/3, 1)),
                   labels = c("small gap", "mid gap", "large gap"),
                   include.lowest = TRUE), by = pair]
summ <- d[, .(median_resid = median(mean_abs_resid_clr),
             q25 = quantile(mean_abs_resid_clr, 0.25),
             q75 = quantile(mean_abs_resid_clr, 0.75)), by = .(pair, tertile)]

p_b <- ggplot(summ, aes(tertile, median_resid, colour = pair, group = pair)) +
  geom_line(linewidth = 0.4, alpha = 0.7) +
  geom_pointrange(aes(ymin = q25, ymax = q75), size = 0.35, fatten = 2.5) +
  scale_colour_manual(values = PAIR_COL, name = NULL) +
  labs(x = "Donor's |sequencing-depth gap| between the two datasets (tertile)",
       y = "Median donor mean |residual CLR|\n(IQR)") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_cns +
  theme(legend.position = "bottom", legend.justification = "center")

fig <- (p_a / p_b) +
  plot_layout(heights = c(1.5, 1.3)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, "figure8_depth_vs_discordance"), height_in = 6.5, width_mm = 220)
