# Figure 9, version 5 (user 2026-11-10): as v4, but panel c shows all 12 test sets
# (bulk cohort x snRNA-seq ground truth) as columns, grouped by cohort; fig9_v4.R untouched.
# and their outputs are left untouched.
#   a  draft schematic (fig9_schematic.R)
#   b  S1 / S2 / S3 per test cohort: median Spearman rho over Core Supertypes on common
#      donors, strategies selected on Spearman (strategy_scores.py METRIC=spearman ->
#      fig9_by_bulk_spearman.tsv)
#   c  benchmark dot matrix: training model x test cohort, dot = median Spearman rho over
#      Core Supertypes (fig9_heatmap_spearman.tsv). Per cohort: mean over the cohort's
#      INDEPENDENT test sets (test dataset not in the training label); a cell with only
#      same-dataset tests shows the 5-fold CV value as an open dot. Bar = mean over the
#      model's independent test cohorts; rows sorted by it.
#   d  per Core Supertype, the recommended model (ROSMAP all 4): mean rho over independent
#      aging/AD test cohorts (MSBB, Diverse; ROSMAP is CV for this model, so excluded) vs
#      over adult/psychiatric test cohorts (CMC, HBCC, NABEC)
#      (fig9_heatmap_spearman_celltype.tsv; mean over ground truths within a cohort).
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig9_v5.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure9"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures")
PT <- 6.5
STRAT_COL <- c(S1 = "#7B8C9E", S2 = "#2F7E8F", S3 = "#E1A22E", prod = "black")
STRAT_LAB <- c(S1 = "S1 paired", S2 = "S2 global", S3 = "S3 class-matched")
KIND_SHP <- c(cross = 16, oof = 1)
BULK_ORD <- c("ROSMAP", "MSBB_FP", "Diverse_Rush2_DLPFC", "CommonMind_MSSM", "NIMH_HBCC", "NABEC")
SHORT <- c(ROSMAP = "ROSMAP", MSBB_FP = "MSBB (FP)", Diverse_Rush2_DLPFC = "Diverse Rush",
           CommonMind_MSSM = "CMC", NIMH_HBCC = "HBCC", NABEC = "NABEC")
COH1 <- c(ROSMAP = "ROSMAP", MSBB_FP = "MSBB", Diverse_Rush2_DLPFC = "Diverse",
          CommonMind_MSSM = "CMC", NIMH_HBCC = "HBCC", NABEC = "NABEC")
CONS_LAB <- c("ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2" = "ROSMAP: all 4",
              "ROSMAP__Green+Mathys+Multiome2025_DLPFC" = "ROSMAP: Green+Mathys+Diverse_mult.",
              "NIMH_HBCC__PsychAD_HBCC+Catching_2026_HBCC" = "HBCC: PsychAD+Catching")
S2_MODEL <- "ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2"
CLS_COL <- c(Excitatory = NEU, Inhibitory = "#8491B4", "Non-neuronal" = GLIA)
small <- theme(axis.text = element_text(size = PT), axis.title = element_text(size = 8),
               strip.text = element_text(size = PT), legend.text = element_text(size = PT),
               legend.title = element_text(size = PT), legend.key.size = unit(0.3, "cm"))

# ---- a: schematic --------------------------------------------------------------
source(file.path(HERE, "fig9_schematic.R"))
pa <- make_schematic(STRAT_COL, pt = 5)

# ---- b: strategies per test cohort (Spearman) -----------------------------------
by <- fread(file.path(DATA, "fig9_by_bulk_spearman.tsv"))
sc <- melt(by, id.vars = c("bulk", "cls", "kind"), measure.vars = paste0(names(STRAT_LAB), "_med_core"),
           variable.name = "strategy", value.name = "r")[!is.na(r) & bulk %in% BULK_ORD]
sc[, strategy := factor(sub("_med_core", "", strategy), names(STRAT_LAB))]
sc[, coh := factor(SHORT[bulk], rev(SHORT[BULK_ORD]))]
sc[, y := as.numeric(coh) + c(S1 = 0.2, S2 = 0, S3 = -0.2)[as.character(strategy)]]
rng <- sc[, .(lo = min(r), hi = max(r), y = as.numeric(coh[1])), by = coh]
pb <- ggplot(sc) +
  geom_segment(data = rng, aes(x = lo, xend = hi, y = y, yend = y), colour = GREYS[1], linewidth = 0.4) +
  geom_point(aes(r, y, colour = strategy, shape = kind), size = 1.8, stroke = 0.6) +
  scale_y_continuous(breaks = seq_along(levels(sc$coh)), labels = levels(sc$coh),
                     expand = expansion(add = 0.4)) +
  scale_colour_manual(values = STRAT_COL[names(STRAT_LAB)], labels = STRAT_LAB, name = NULL) +
  scale_shape_manual(values = KIND_SHP, name = NULL,
                     labels = c(cross = "other snRNA-seq", oof = "same snRNA-seq (CV)")) +
  labs(x = "Median Spearman ρ (Core Supertypes)", y = "Test cohort") +
  guides(colour = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1)) +
  theme_cns + small + theme(legend.position = "bottom", legend.box = "vertical",
                            legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))

# ---- c: benchmark dot matrix, 12 test sets ---------------------------------------
# Column = one test set (bulk cohort x snRNA-seq ground truth), grouped by cohort as in
# figure9_supp_matrix. Dot = median Spearman rho over Core Supertypes for that cell
# (fig9_heatmap_spearman.tsv). Open circle = test dataset is part of the training label
# (5-fold CV): filled with its value, marked by a black centre dot. Black ring = best
# model for that test set among independent cells.
# Bar = mean over test cohorts of the cohort's mean over independent test sets (each
# cohort counts once, as in the strategy selection); rows sorted by it. Teal band = the
# recommended model (S2 pick on Spearman: ROSMAP all 4).
SN_ORD <- c("Green", "Mathys", "Multiome2025_DLPFC", "CUIMC2", "PsychAD_MSBB", "PsychAD_ROSMAP",
            "BrainSCOPE_CMC", "PsychAD_HBCC", "Catching_2026_HBCC", "Catching_2026_NABEC")
sn_lab <- function(s) sub("Multiome2025_DLPFC", "Diverse_multiome", s)
tt <- fread(file.path(DATA, "fig9_heatmap_spearman.tsv"))
tt[, `:=`(tb = sub("__.*", "", trainer), ts = sub(".*?__", "", trainer),
          ub = sub("__.*", "", target), us = sub(".*?__", "", target))]
tt[, cons := grepl("+", ts, fixed = TRUE)]
if (!"cv" %in% names(tt)) tt[, cv := FALSE]
tt[is.na(cv), cv := FALSE]
tt <- tt[ub %in% BULK_ORD & tb %in% BULK_ORD & (!cons | trainer %in% names(CONS_LAB))]
tt[, indep := !(diag | cv)]
ov <- tt[indep == TRUE, .(m = mean(med_core)), by = .(trainer, tb, ts, cons, ub)][
  , .(overall = mean(m), n_coh = .N), by = .(trainer, tb, ts, cons)]
ov[, lab := ifelse(cons, CONS_LAB[trainer], paste0(COH1[tb], ": ", sn_lab(ts)))]
ov[, lab_f := factor(lab, ov[order(overall)]$lab)]
tt <- ov[, .(trainer, lab_f)][tt, on = "trainer"]
tt[, coh := factor(COH1[ub], COH1[BULK_ORD])]
tt[, test := factor(sn_lab(us), unique(sn_lab(SN_ORD)))]
tt[indep == TRUE, best := med_core == max(med_core), by = target]
tt[is.na(best), best := FALSE]
tt[indep == TRUE, rk := frank(-med_core, ties.method = "min"), by = target]
ov <- tt[indep == TRUE, .(top3 = sum(rk <= 3), n_test = .N), by = trainer][ov, on = "trainer"]
fill_sc <- scale_fill_gradient2(low = "#2F4FA8", mid = "white", high = "#B2182B", midpoint = 0,
                                limits = c(-0.1, 0.3), oob = scales::squish, name = "Spearman ρ")
rec_row <- ov[trainer == S2_MODEL]$lab_f
pc_dot <- ggplot(tt, aes(test, lab_f)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(rec_row) - 0.5,
           ymax = as.numeric(rec_row) + 0.5, fill = STRAT_COL[["S2"]], alpha = 0.15) +
  geom_point(aes(size = pmax(med_core, 0), fill = med_core), shape = 21, colour = "grey45",
             stroke = 0.25) +
  geom_point(data = tt[indep == FALSE], aes(shape = "5-fold CV (test set in training data)"),
             size = 0.55, colour = "black") +
  geom_point(data = tt[best == TRUE], aes(size = pmax(med_core, 0),
             shape = "Best model for this test set"), colour = "black", stroke = 0.75) +
  scale_shape_manual(values = c("5-fold CV (test set in training data)" = 16,
                                "Best model for this test set" = 1), name = NULL) +
  scale_size_area(max_size = 3.6, limits = c(0, 0.3), oob = scales::squish, guide = "none") +
  fill_sc +
  facet_grid(. ~ coh, scales = "free_x", space = "free_x") +
  labs(x = "Test set: snRNA-seq ground truth, by bulk cohort", y = "Training model") +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(2.4, "cm"), barheight = unit(0.22, "cm"),
                                title.position = "left", title.vjust = 0.9),
         shape = guide_legend(order = 2, override.aes = list(size = c(0.9, 2.4), stroke = c(0.5, 0.75)))) +
  theme_cns + small + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust = 1, size = PT - 0.5),
                            strip.text = element_text(size = PT), strip.background = element_blank(),
                            panel.spacing = unit(2, "pt"),
                            panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
                            legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
pc_bar <- ggplot(ov, aes(overall, lab_f, fill = cons)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(rec_row) - 0.5,
           ymax = as.numeric(rec_row) + 0.5, fill = STRAT_COL[["S2"]], alpha = 0.15) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", overall)), hjust = -0.15, size = 5.5 / .pt) +
  scale_fill_manual(values = c(`TRUE` = STRAT_COL[["S2"]], `FALSE` = GREYS[1]), guide = "none") +
  scale_x_continuous(limits = c(0, 0.31), breaks = c(0, 0.1, 0.2), expand = expansion(mult = 0)) +
  labs(x = "Mean ρ over\ntest cohorts", y = NULL) +
  theme_cns + small + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                            axis.line.y = element_blank())
pc_top <- ggplot(ov, aes(0, lab_f)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(rec_row) - 0.5,
           ymax = as.numeric(rec_row) + 0.5, fill = STRAT_COL[["S2"]], alpha = 0.15) +
  geom_text(aes(label = sprintf("%d / %d", top3, n_test), colour = cons), size = PT / .pt) +
  scale_colour_manual(values = c(`TRUE` = STRAT_COL[["S2"]], `FALSE` = "grey25"), guide = "none") +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = NULL) +
  labs(x = "Top 3 in\ntest sets", y = NULL) +
  theme_cns + small + theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
                            axis.line = element_blank())
pc <- pc_dot + pc_bar + pc_top + plot_layout(widths = c(6, 1.3, 0.55))

# ---- d: per supertype, aging/AD vs adult/psychiatric test cohorts -------------------
tax <- fread("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv")
tax[, cls := fcase(grepl("Glutamatergic", class_label), "Excitatory",
                   grepl("GABAergic", class_label), "Inhibitory", default = "Non-neuronal")]
pt <- fread(file.path(DATA, "fig9_heatmap_spearman_celltype.tsv"))[trainer == S2_MODEL]
pt[, ub := sub("__.*", "", target)]
pt <- pt[ub %in% BULK_ORD, .(rho = mean(rho)), by = .(ub, cell_type)]
AGE <- c("MSBB_FP", "Diverse_Rush2_DLPFC"); PSY <- c("CommonMind_MSSM", "NIMH_HBCC", "NABEC")
sx <- pt[, .(x = mean(rho[ub %in% AGE]), y = mean(rho[ub %in% PSY]),
             nx = sum(ub %in% AGE), ny = sum(ub %in% PSY)), by = cell_type]
sx <- sx[nx == length(AGE) & ny == length(PSY)]
sx <- tax[, .(cell_type = supertype_id, label = supertype_label, cls)][sx, on = "cell_type", nomatch = 0]
sx[, cls := factor(cls, names(CLS_COL))]
st <- cor.test(sx$x, sx$y, method = "spearman", exact = FALSE)
below <- mean(sx$y < sx$x)
lim <- range(c(sx$x, sx$y, 0)) + c(-0.02, 0.02)
lab_pts <- rbind(sx[order(-(x + y))][1:3], sx[order(y - x)][1:2])   # best overall + largest drops
pd <- ggplot(sx, aes(x, y)) +
  geom_hline(yintercept = 0, colour = GREYS[1], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GREYS[1], linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_point(aes(fill = cls), shape = 21, colour = "white", stroke = 0.2, size = 1.9, alpha = 0.9) +
  ggrepel::geom_text_repel(data = lab_pts, aes(label = label), size = 5.5 / .pt, min.segment.length = 0,
                           segment.size = 0.2, segment.colour = "grey50", box.padding = 0.35,
                           max.overlaps = Inf, seed = 1) +
  annotate("text", x = lim[1] + 0.01, y = lim[2] - 0.005, hjust = 0, vjust = 1, size = PT / .pt,
           label = sprintf("Spearman ρ = %.2f\n%.0f%% of supertypes below the line", st$estimate, 100 * below)) +
  scale_fill_manual(values = CLS_COL, name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "ρ in aging/AD test cohorts\n(MSBB, Diverse)",
       y = "ρ in adult/psychiatric test cohorts\n(CMC, HBCC, NABEC)") +
  guides(fill = guide_legend(override.aes = list(size = 2.2))) +
  theme_cns + small + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- assemble --------------------------------------------------------------------
TAG <- theme(plot.tag = element_text(size = 10, face = "bold", colour = "black"))
top_row <- (pa + labs(tag = "a") + wrap_elements(full = pb) + labs(tag = "b") +
              plot_layout(widths = c(1.5, 1))) & TAG
# v5 main = a | b over c (full width: 12 test sets need the page). Panel d is saved as a
# standalone file until its role is decided (user question 2026-11-10).
fig <- (wrap_elements(full = top_row) / (wrap_elements(full = pc) + labs(tag = "c")) +
          plot_layout(heights = c(2.5, 3.1))) & TAG
ggsave(file.path(FIG, "figure9_supertype_scatter_standalone.png"), pd, width = 95 / 25.4,
       height = 95 / 25.4, dpi = 300, bg = "white")
ggsave(file.path(FIG, "figure9_supertype_scatter_standalone.pdf"), pd, width = 95 / 25.4,
       height = 95 / 25.4, bg = "white", device = cairo_pdf)
save_fig(fig, file.path(FIG, "figure9_celmod_strategies_v5"), height_in = 5.8)
ggsave(file.path(FIG, "figure9_celmod_strategies_v5.svg"), fig, width = 183 / 25.4, height = 5.8,
       bg = "white", device = svglite::svglite)
print(ov[order(-overall), .(lab, overall = round(overall, 3), n_coh, top3, n_test)])
print(tt[best == TRUE, .(target, lab_f, med_core = round(med_core, 3))][order(target)])
message(sprintf("d: n = %d supertypes; Spearman %.2f (p = %.1e); below diagonal %.0f%%",
                nrow(sx), st$estimate, st$p.value, 100 * below))
print(sx[, .(below = round(mean(y < x), 2), med_drop = round(median(y - x), 3)), by = cls])
