# Diagnostic (not a manuscript figure): are the subclass values in Fig 1a right?
#   top row: the per-pair values behind each subclass square and its
#            top supertype, labelled by pair, with the plotted median (black
#            diamond) and IQR (black bar) overlaid
#   bottom:  subclass CLR scatters for the three pairs with donor-level data;
#            each panel prints the statistic recomputed from the points and the
#            value in 114_pair_percelltype_subclass.csv
# Metric: "icc" (default; symmetric ICC = 2 cov / (var_a + var_b), the paper's metric
#   since 2026-10-05), "rho" (Spearman), "r" (Pearson), or "icc_a" (icc_person_a =
#   reference-dependent statistic of the 114 sweep). The other versions are
#   for comparison only; the top supertypes are fixed to the ones Fig 1 shows,
#   so the ranking is not re-derived.
# Inputs: data/fig1a_pairs.tsv, data/fig1bc_clr.tsv (from fig1_data.py)
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript check_fig1a_subclass.R [icc|rho|r|icc_a]

suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork) })
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure1"
source(file.path(dirname(HERE), "common", "theme_cns.R"))

METRIC <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "icc"
stopifnot(METRIC %in% c("icc", "rho", "r", "icc_a"))
stat_fun <- switch(METRIC,
  icc = function(ac, bc) 2 * cov(ac, bc) / (var(ac) + var(bc)),
  rho = function(ac, bc) cor(ac, bc, method = "spearman"),
  r   = function(ac, bc) cor(ac, bc),
  icc_a = function(ac, bc) cov(ac, bc) / var(ac))
YLAB <- switch(METRIC,
  icc = "Symmetric ICC per pair (2 cov / (var_a + var_b))",
  rho = "Spearman \u03c1 per pair (from donor CLR)",
  r   = "Pearson r per pair (r_cross_study, 114)",
  icc_a = "ICC per pair (icc_person_a = cov / var(a), 114)")
SUFFIX <- if (METRIC == "icc") "" else paste0("_", METRIC)

SUBCLASSES <- c("Immune", "Sst", "L2/3 IT")
TOP <- c(Immune = "Micro-PVM_2_3-SEAAD", Sst = "Sst_25", `L2/3 IT` = "L2/3 IT_13")
# study-pair colours as in Fig 3 (cnsplots "Cell")
PAIR_COL <- c("Green x Mathys" = "#2F7E8F", "Green x Multiome" = "#E1A22E",
              "Mathys x Multiome" = "#8B6FA8", "Multiome x PsychAD" = "#C84C3A",
              "BrainSCOPE CMC x PsychAD" = "#5F9862",
              "PsychAD x Ruzicka (MSSM2 x MSSM1)" = "#4E5A8A")

a <- fread(file.path(HERE, "data/fig1a_pairs.tsv"))
a <- a[cell_type %in% c(SUBCLASSES, TOP)]
a[, subclass := factor(subclass, levels = SUBCLASSES)]
a[, cell_type := factor(cell_type, levels = c(rbind(SUBCLASSES, TOP)))]
a[, pair := factor(pair, levels = names(PAIR_COL))]
a[, val := switch(METRIC, icc = icc_sym, rho = rho, r = r, icc_a = icc_person_a)]
sm <- a[, .(n = .N, median = median(val), q25 = quantile(val, .25), q75 = quantile(val, .75)),
        by = .(subclass, cell_type)]
cat(sprintf("metric = %s; values drawn as the squares/circles:\n", METRIC)); print(sm, digits = 3)

p1 <- ggplot(a, aes(cell_type, val)) +
  geom_linerange(data = sm, aes(cell_type, ymin = q25, ymax = q75), inherit.aes = FALSE,
                 colour = "black", linewidth = 1.6, alpha = 0.25) +
  geom_point(aes(colour = pair), size = 2, position = position_nudge(x = 0.22)) +
  ggrepel::geom_text_repel(aes(label = sprintf("%.2f", val), colour = pair), size = 2.3,
            nudge_x = 0.5, direction = "y", hjust = 0, segment.size = 0.2,
            min.segment.length = 0, show.legend = FALSE, seed = 1) +
  geom_point(data = sm, aes(cell_type, median), inherit.aes = FALSE,
             shape = 18, size = 3.2, colour = "black") +
  geom_text(data = sm, aes(cell_type, median, label = sprintf("median %.2f\n(n = %d)", median, n)),
            inherit.aes = FALSE, size = 2.3, hjust = 1, nudge_x = -0.1, colour = "black") +
  scale_colour_manual(values = PAIR_COL, name = NULL) +
  scale_x_discrete(expand = expansion(add = c(1.4, 1.2))) +
  facet_wrap(~ subclass, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = YLAB) +
  guides(colour = guide_legend(nrow = 2)) +
  theme_cns + theme(legend.position = "bottom",
                    axis.text.x = element_text(angle = 20, hjust = 1))

clr <- fread(file.path(HERE, "data/fig1bc_clr.tsv"))[level == "subclass" & cell_type %in% SUBCLASSES]
clr[, cell_type := factor(cell_type, levels = SUBCLASSES)]
clr[, pair := factor(pair, levels = names(PAIR_COL))]
st <- clr[, { ac <- clr_a - mean(clr_a); bc <- clr_b - mean(clr_b)
              list(v_pts = stat_fun(ac, bc), n = .N) }, by = .(pair, cell_type)]
st <- merge(st, a[, .(pair, cell_type, v_114 = val)], by = c("pair", "cell_type"))
st[, lab := sprintf("from points %.3f\n114 table  %.3f\nn = %d", v_pts, v_114, n)]
cat(sprintf("\nsubclass %s recomputed from the plotted donors vs 114:\n", METRIC))
print(st[, .(pair, cell_type, v_pts, v_114, n)], digits = 4)
stopifnot(max(abs(st$v_pts - st$v_114)) < 1e-9)

p2 <- ggplot(clr, aes(clr_a, clr_b)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_point(aes(colour = pair), size = 0.8, alpha = 0.7, stroke = 0, show.legend = FALSE) +
  geom_text(data = st, aes(label = lab), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
            size = 2.3, inherit.aes = FALSE) +
  scale_colour_manual(values = PAIR_COL) +
  facet_grid(pair ~ cell_type, scales = "free", labeller = labeller(pair = label_wrap_gen(16))) +
  labs(x = "CLR, dataset a (ln)", y = "CLR, dataset b (ln)") +
  theme_cns

fig <- (p1 / p2) + plot_layout(heights = c(1.1, 2)) + plot_annotation(tag_levels = "a")
save_fig(fig, file.path(HERE, paste0("figures/check_fig1a_subclass", SUFFIX)),
         height_in = 10, width_mm = 200)
