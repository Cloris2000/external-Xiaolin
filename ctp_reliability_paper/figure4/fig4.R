# Figure 4 (DRAFT). Deconvolution works where snRNA-seq is reliable, and a
# reliability-only threshold (no bulk data in the selection) picks out a
# deconvolution-ready set of supertypes.
#   a  per supertype: out-of-fold CelMod accuracy (bulk vs matched snRNA, median
#      over 6 pairs) against snRNA-seq reliability (symmetric ICC, median over
#      5 matched-donor pairs -- the paper's decided metric, common/data/
#      reliability_summary.tsv; same metric as Fig 1)
#   b  out-of-fold CelMod r by group (Core Supertypes = lower bound of the donor-
#      bootstrap 95% CI of the median snRNA-seq ICC (common/reliability_ci.py) > 0; vs non-core supertypes,
#      vs subclasses; user decision 2026-10-27, was median ICC > 0), per bulk x snRNA pair -- direct test of the panel a
#      relationship, selection uses no bulk data; cohort strip labels are the
#      display labels of common/cohort_names.tsv
# Supplement (not plotted here): per-pair Spearman rho(ICC, out-of-fold r),
#   data/fig4_supp_rho_by_pair.tsv, shows the panel a relationship holds in
#   each of the 6 pairs (rho 0.36-0.63, all p < 1e-5).
# Redesign 2026-10-14 (dropped old panel b, ROSMAP-triad triangulation
# reliability; switched panel a off the old 6-pair icc_person_a metric; the
# per-pair rho forest plot was panel c briefly and moved to the supplement
# because panel b already carries the per-cohort message).
# Taxonomy: SEA-AD DFC 2026. Inputs from fig4_data.py.
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig4.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure4"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures"); dir.create(FIG, showWarnings = FALSE)

# ---- a ---------------------------------------------------------------------
library(ggrepel)
# Allen-style hues, matching Figure 1's class colours (user decision 2026-10-05
# for Fig 1): excitatory = red/orange, inhibitory = blue/purple, non-neuronal = green.
CLASS3_COL <- c("Excitatory" = "#D1491C", "Inhibitory" = "#2F4FA8", "Non-neuronal" = "#2E8B3C")
a <- fread(file.path(DATA, "fig4a_icc_vs_deconv.tsv"))[!is.na(bulk_r)]
a[, key := factor(fifelse(class_id == "Glia", "Non-neuronal", class_id), levels = names(CLASS3_COL))]
# filled = Core Supertype (lower bound of the bootstrap 95% CI of the median ICC > 0), open = non-core
ICC_POS_SHP <- c("Core Supertype" = 16, "Non-core Supertype" = 1)
a[, icc_pos := factor(fifelse(ci_lo > 0, names(ICC_POS_SHP)[1], names(ICC_POS_SHP)[2]),
                      levels = names(ICC_POS_SHP))]
cat(sprintf("a: Core Supertype (CI lower bound > 0) in %d/%d supertypes\n", sum(a$ci_lo > 0), nrow(a)))
THRESH <- 0.3
ta <- spear(a$median_icc_sym, a$bulk_r)
cat(sprintf("a: n=%d rho=%.2f p=%s\n", nrow(a), ta$estimate, fmt_p(ta$p.value)))
p_a <- ggplot(a[order(as.integer(key))], aes(median_icc_sym, bulk_r)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", fill = "grey85",
              linewidth = 0.4, linetype = "dashed") +
  geom_point(aes(colour = key, shape = icc_pos), size = 1.3, stroke = 0.45) +
  geom_text_repel(data = a[median_icc_sym >= THRESH & bulk_r >= THRESH],
                  aes(label = cell_type), size = fontsize / .pt * 0.75,
                  min.segment.length = 0, segment.size = 0.25,
                  max.overlaps = Inf, seed = 1, ylim = c(NA, 0.73)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = fontsize / .pt,
           label = sprintf("ρ = %.2f, p = %s", ta$estimate, fmt_p(ta$p.value))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  scale_colour_manual(values = CLASS3_COL, name = NULL) +
  scale_shape_manual(values = ICC_POS_SHP, name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)),
         shape = guide_legend(order = 2, override.aes = list(colour = "grey30"))) +
  labs(x = "snRNA-seq reliability of the supertype\n(median ICC across matched-donor pairs)",
       y = "Deconvolution accuracy of the supertype\n(median out-of-fold r, bulk vs snRNA-seq)") +
  theme_cns + theme(legend.position = "inside", legend.position.inside = c(1, 0),
                    legend.justification = c(1, 0), legend.spacing.y = unit(0, "pt"),
                    legend.margin = margin(0, 0, 0, 0))

# ---- b ---------------------------------------------------------------------
b <- fread(file.path(DATA, "fig4b_oof_by_concordance.tsv"))
# Subclasses stay in data/fig4b_oof_by_concordance.tsv and fig4b_pvals.tsv (resolution
# question, not this panel's claim; user decision 2026-10-27) but are not drawn.
SET_LEV <- c("Non-core Supertype", "Core Supertype")   # as written by fig4_data.py
SET_LAB <- c("Non-core Supertype", "Core Supertype")   # same wording as panel a
b <- b[set %in% SET_LEV]
b[, set_shape := factor(set, levels = SET_LEV, labels = SET_LAB)]
b[, set := factor(set, levels = SET_LEV)]
stopifnot(!anyNA(b$set))

# Displayed dataset names are the display_label of common/cohort_names.tsv
# (same convention as fig8.R), looked up by the wp3 source_key of each side;
# the "bulk x snRNA" pair string stays a data key only.
cohorts <- fread(file.path(HERE, "..", "common", "cohort_names.tsv"))
canon <- function(key) {
  nm <- cohorts$display_label[cohorts$source_key == key]
  stopifnot(length(nm) == 1)
  nm
}
SN_KEY <- c("Green_DFC" = "Green_DFC_person", "Mathys_DFC" = "Mathys_DFC_person",
            "Multiome_DFC" = "Multiome2025_DLPFC_DFC_person")
wrap_name <- function(nm, width = 15) {   # greedy wrap at underscores for ~17 mm strips
  parts <- strsplit(nm, "_")[[1]]
  done <- character(0); cur <- parts[1]
  for (p in parts[-1]) {
    if (nchar(cur) + 1 + nchar(p) <= width) cur <- paste0(cur, "_", p)
    else { done <- c(done, paste0(cur, "_")); cur <- p }
  }
  paste(c(done, cur), collapse = "\n")
}
pairs <- unique(b[, .(pair, bulk, snRNA)])
pairs[, lab := vapply(seq_len(.N), function(i) paste0(
  wrap_name(canon(bulk[i]), 13), "\n×\n",
  wrap_name(canon(if (snRNA[i] %in% names(SN_KEY)) SN_KEY[[snRNA[i]]] else snRNA[i]), 13)), "")]
PAIR_LAB <- setNames(pairs$lab, pairs$pair)
b[, pair := factor(PAIR_LAB[pair], levels = PAIR_LAB)]

# Mann-Whitney one-sided (Core Supertypes > comparison), computed in
# fig4_data.py; drawn as brackets above each facet's points. x = 2 is Core Supertype.
pv <- fread(file.path(DATA, "fig4b_pvals.tsv"))
pv[, pair := factor(PAIR_LAB[pair], levels = PAIR_LAB)]
pv <- pv[comparison == "Non-core Supertype"]
stopifnot(!anyNA(pv$pair), nrow(pv) == uniqueN(b$pair))
pv[, x0 := match(comparison, SET_LEV)]
pv[, sig := fcase(p < 1e-3, "***", p < 1e-2, "**", p < 0.05, "*", default = "n.s.")]
ymax <- b[, .(ymax = max(r, na.rm = TRUE)), by = pair]
pv <- merge(pv, ymax, by = "pair")
setorder(pv, pair, -x0)
pv[, bracket_y := ymax + (seq_len(.N) - 1) * 0.11 + 0.06, by = pair]

p_b <- ggplot(b, aes(set, r)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_jitter(aes(shape = set_shape), colour = "grey30", width = 0.2, size = 0.9, alpha = 0.8, stroke = 0.35) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black", linewidth = 0.3, width = 0.55) +
  geom_segment(data = pv, aes(x = x0, xend = 2, y = bracket_y, yend = bracket_y),
              colour = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(data = pv, aes(x = (x0 + 2) / 2, y = bracket_y, label = sig),
           colour = "black", size = fontsize / .pt, vjust = -0.15, inherit.aes = FALSE) +
  facet_wrap(~ pair, nrow = 1) +
  scale_shape_manual(values = setNames(unname(ICC_POS_SHP)[2:1], SET_LAB), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Out-of-fold r, deconvolved vs snRNA-seq") +
  guides(shape = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
                    legend.text = element_text(size = fontsize * 0.85),
                    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                    strip.text = element_text(size = fontsize * 0.6, lineheight = 0.9),
                    strip.clip = "off",
                    plot.caption = element_text(size = fontsize * 0.7, colour = "grey30", hjust = 0),
                    panel.spacing.x = unit(3, "pt"))

fig <- wrap_plots(A = free(p_a, type = "panel"), B = free(p_b, type = "panel"), nrow = 1) +
  plot_layout(widths = c(1, 1.4)) + plot_annotation(tag_levels = "a")
save_fig(fig, file.path(FIG, "figure4_resolution_draft"), height_in = 4.2)
