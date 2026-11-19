# Figure 9, version 4 (user 2026-11-10). Everything on Spearman rho; fig9.R / _v2 / _v3
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
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig9_v4.R
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

# ---- c: benchmark dot matrix ---------------------------------------------------------
tt <- fread(file.path(DATA, "fig9_heatmap_spearman.tsv"))
tt[, `:=`(tb = sub("__.*", "", trainer), ts = sub(".*?__", "", trainer), ub = sub("__.*", "", target))]
tt[, cons := grepl("+", ts, fixed = TRUE)]
if (!"cv" %in% names(tt)) tt[, cv := FALSE]
tt[is.na(cv), cv := FALSE]
tt <- tt[ub %in% BULK_ORD & tb %in% BULK_ORD & (!cons | trainer %in% names(CONS_LAB))]
tt[, indep := !(diag | cv)]
dm <- tt[, if (any(indep)) .(rho = mean(med_core[indep]), open = FALSE)
           else .(rho = mean(med_core), open = TRUE), by = .(trainer, tb, ts, cons, ub)]
ov <- dm[open == FALSE, .(overall = mean(rho), n_coh = .N), by = .(trainer, tb, ts, cons)]
ov[, lab := ifelse(cons, CONS_LAB[trainer],
                   paste0(COH1[tb], ": ", sub("Multiome2025_DLPFC", "Diverse_multiome", ts)))]
ov[, lab_f := factor(lab, ov[order(overall)]$lab)]
dm <- ov[, .(trainer, lab_f)][dm, on = "trainer"]
dm[, coh := factor(COH1[ub], COH1[BULK_ORD])]
fill_sc <- scale_fill_gradient2(low = "#2F4FA8", mid = "white", high = "#B2182B", midpoint = 0,
                                limits = c(-0.1, 0.3), oob = scales::squish, name = "Spearman ρ")
pc_dot <- ggplot(dm, aes(coh, lab_f)) +
  geom_point(data = dm[open == FALSE], aes(size = pmax(rho, 0), fill = rho), shape = 21,
             colour = "grey30", stroke = 0.25) +
  geom_point(data = dm[open == TRUE], aes(size = pmax(rho, 0)), shape = 21, fill = NA,
             colour = "black", stroke = 0.45) +
  scale_size_area(max_size = 4.2, limits = c(0, 0.3), oob = scales::squish, guide = "none") +
  fill_sc +
  scale_x_discrete(position = "top") +
  labs(x = "Test cohort", y = "Training model") +
  guides(fill = guide_colourbar(barwidth = unit(2.4, "cm"), barheight = unit(0.22, "cm"),
                                title.position = "top")) +
  theme_cns + small + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                            axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = PT),
                            panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
                            legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
pc_bar <- ggplot(ov, aes(overall, lab_f, fill = cons)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", overall)), hjust = -0.15, size = 5.5 / .pt) +
  scale_fill_manual(values = c(`TRUE` = STRAT_COL[["S2"]], `FALSE` = GREYS[1]), guide = "none") +
  scale_x_continuous(limits = c(0, 0.31), breaks = c(0, 0.1, 0.2), expand = expansion(mult = 0),
                     position = "top") +
  labs(x = "Mean ρ", y = NULL) +
  theme_cns + small + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                            axis.line.y = element_blank())
pc <- pc_dot + pc_bar + plot_layout(widths = c(6, 2.8))

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
bot_row <- (wrap_elements(full = pc) + labs(tag = "c") + wrap_elements(full = pd) + labs(tag = "d") +
              plot_layout(widths = c(1.35, 1))) & TAG
fig <- (wrap_elements(full = top_row) / bot_row + plot_layout(heights = c(2.5, 3.3))) & TAG
save_fig(fig, file.path(FIG, "figure9_celmod_strategies_v4"), height_in = 6.0)
ggsave(file.path(FIG, "figure9_celmod_strategies_v4.svg"), fig, width = 183 / 25.4, height = 6.0,
       bg = "white", device = svglite::svglite)
print(ov[order(-overall), .(lab, overall = round(overall, 3), n_coh)])
message(sprintf("d: n = %d supertypes; Spearman %.2f (p = %.1e); below diagonal %.0f%%",
                nrow(sx), st$estimate, st$p.value, 100 * below))
print(sx[, .(below = round(mean(y < x), 2), med_drop = round(median(y - x), 3)), by = cls])
