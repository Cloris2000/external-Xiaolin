# Figure 9 (DRAFT). Which CelMod training strategy generalises across bulk cohorts?
# Strategies (celmod_strategies/strategy_scores.py, decisions 2026-11-04):
#   S1 paired   train on the target bulk cohort's own (other) snRNA-seq
#   S2 global   one recipe for every cohort, chosen leave-one-cohort-family-out
#   S3 matched  as S2 within diagnosis class (aging/AD vs adult/psychiatric)
#   prod        Shreejoy's production trio, ROSMAP x Green+Mathys+Multiome (reference)
# Ground truth = a DIFFERENT snRNA-seq dataset of the target's own people where the
# cohort has one ("cross"); otherwise the same dataset, out of fold ("oof", open symbols,
# favours S1). Every comparison within a unit is on one common donor set (all
# strategies' training people removed). Metric: median Pearson r over the 86 Fig 4 Core
# Supertypes (ci_lo > 0).
#   a  empty slot for the BioRender schematic
#   b  trainer x target median Core r (scores.csv, each trainer's own hold-out donors)
#   c  per bulk cohort: S1 / S2 / S3 / prod on common donors
#   d  per Core supertype: r of the best strategy vs median symmetric ICC (link to Fig 4)
# Inputs: figure9/data/fig9_*.tsv from celmod_strategies/strategy_scores.py.
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig9.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure9"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures"); dir.create(FIG, showWarnings = FALSE)

STRAT_LEV <- c("S1", "S2", "S3", "prod")
STRAT_LAB <- c(S1 = "S1 paired", S2 = "S2 global", S3 = "S3 class-matched", prod = "Production trio")
STRAT_COL <- c(S1 = "#7B8C9E", S2 = "#2F7E8F", S3 = "#E1A22E", prod = "black")  # fig5 palette
CLS_LAB <- c(aging_AD = "Aging / AD", adult_psych = "Adult / psychiatric")
KIND_SHP <- c(cross = 16, oof = 1)                                              # filled vs open

pair_lab <- function(p) {          # "ROSMAP__Green+Mathys" -> "ROSMAP: Green+Mathys"
  b <- sub("__.*", "", p); s <- sub(".*?__", "", p)
  s <- vapply(strsplit(s, "+", fixed = TRUE), function(m) paste(cohort_label(m), collapse = "+"), "")
  paste0(cohort_label(b), ": ", s)
}

# ---- b: trainer x target -------------------------------------------------------
# Square matrix: one bulk region per cohort (strategy_scores.REP_BULK), rows and columns in
# the SAME order -- data bank (AMP-AD, PsychENCODE, NABEC), then cohort, then snRNA-seq
# dataset -- so the circled diagonal is each pair's own out-of-fold score. Consensus models
# (several snRNA-seq labels on one bulk) are an extra block of rows. Marginal strips:
# top = how easy each test set is (median over single-dataset models, diagonal excluded);
# right = how good each model is (median over test sets, diagonal excluded).
BULK_ORD <- c("ROSMAP", "MSBB_FP", "Diverse_Rush2_DLPFC", "CommonMind_MSSM", "NIMH_HBCC", "NABEC")
BANK <- c(ROSMAP = "AMP-AD", MSBB_FP = "AMP-AD", Diverse_Rush2_DLPFC = "AMP-AD",
          CommonMind_MSSM = "PsychENCODE", NIMH_HBCC = "PsychENCODE", NABEC = "NABEC")
SN_ORD <- c("Green", "Mathys", "Multiome2025_DLPFC", "CUIMC2", "PsychAD_MSBB", "PsychAD_ROSMAP",
            "BrainSCOPE_CMC", "PsychAD_HBCC", "Catching_2026_HBCC", "Catching_2026_NABEC")
SHORT <- c(ROSMAP = "ROSMAP", MSBB_FP = "MSBB (FP)", Diverse_Rush2_DLPFC = "Diverse Rush",
           CommonMind_MSSM = "CMC", NIMH_HBCC = "HBCC", NABEC = "NABEC")
SHORT["MSBB_FP"] <- "MSBB\n(FP)"
STRIP <- setNames(ifelse(BANK[BULK_ORD] == SHORT[BULK_ORD], SHORT[BULK_ORD],
                         sprintf("%s\n%s", BANK[BULK_ORD], SHORT[BULK_ORD])), BULK_ORD)  # bank + cohort; one line when equal (NABEC)
CSTRIP <- setNames(sub(" (FP)", "", unname(SHORT[BULK_ORD]), fixed = TRUE), BULK_ORD)                      # column strips
CONS_LAB <- c("ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2" = "ROSMAP: all 4 (S2)",
              "ROSMAP__Green+Mathys+Multiome2025_DLPFC" = "ROSMAP: Green+Mathys+Diverse_multiome",
              "NIMH_HBCC__PsychAD_HBCC+Catching_2026_HBCC" = "HBCC: PsychAD+Catching")
tt <- fread(file.path(DATA, "fig9_heatmap_spearman.tsv"))   # Spearman (user 2026-11-04); heatmap_spearman.py
tt[, `:=`(tb = sub("__.*", "", trainer), ts = sub(".*?__", "", trainer),
          ub = sub("__.*", "", target), us = sub(".*?__", "", target))]
tt <- tt[ub %in% BULK_ORD & tb %in% BULK_ORD & (!grepl("+", ts, fixed = TRUE) | trainer %in% names(CONS_LAB))]
pairs <- unique(tt[!grepl("+", ts, fixed = TRUE), .(b = tb, s = ts)])
pairs <- pairs[order(match(b, BULK_ORD), match(s, SN_ORD))]
sn_lab <- function(s) cohort_label(s)
tt[, cons := grepl("+", ts, fixed = TRUE)]
tt[, row_f := factor(ifelse(cons, "Consensus\nmodels", STRIP[tb]), c(STRIP, "Consensus\nmodels"))]
tt[, col_f := factor(STRIP[ub], STRIP)]   # bank + cohort, as the row strips
lev_tick <- unique(c(sn_lab(pairs$s), CONS_LAB))
tt[, tr_tick := factor(ifelse(cons, CONS_LAB[trainer], sn_lab(ts)), rev(lev_tick))]
tt[, tg_tick := factor(sn_lab(us), lev_tick)]
fill_sc <- scale_fill_gradient2(low = "#2F4FA8", mid = "white", high = "#B2182B", midpoint = 0,
                                limits = c(-0.1, 0.35), oob = scales::squish,
                                name = "Median Spearman \u03c1\n(Core Supertypes)")
if (!"cv" %in% names(tt)) tt[, cv := FALSE]
tt[is.na(cv), cv := FALSE]
# circled = scored against a dataset that is part of the training label: the 5-fold
# out-of-fold diagonal, and the consensus models' out-of-fold cells (cons_oof.py).
# Cells with < 50 test donors after the person hold-out are left empty (user, 2026-11-04).
tt[, circled := diag | cv]
CIRC_LAB <- "5-fold cross-validation\n(test dataset in training label)"
pb_main <- ggplot(tt, aes(tg_tick, tr_tick, fill = med_core)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = tt[circled == TRUE], aes(shape = CIRC_LAB), size = 2.6, stroke = 0.5,
             colour = "black") +
  scale_shape_manual(values = 1, name = NULL) +
  facet_grid(row_f ~ col_f, scales = "free", space = "free", switch = "both") + fill_sc +
  labs(x = "Test set: bulk RNA-seq cohort and snRNA-seq ground truth",
       y = "Training set: bulk RNA-seq cohort and snRNA-seq reference") +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3.5, "cm"), barheight = unit(0.3, "cm"), title.position = "top"),
         shape = guide_legend(order = 2)) +
  theme_cns + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                    axis.text.y = element_text(size = 7), panel.spacing = unit(2, "pt"),
                    strip.text = element_text(size = 7), strip.text.y.left = element_text(angle = 0),
                    strip.placement = "outside", legend.position = "bottom",
                    legend.box = "horizontal", legend.justification = "center", legend.spacing.x = unit(0.8, "cm"),
                    strip.text.x.bottom = element_text(size = 6),
                    axis.title = element_text(size = 8),
                    legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
                    axis.line = element_blank(), axis.ticks = element_blank())
single <- tt[cons == FALSE & diag == FALSE]
top <- single[, .(m = median(med_core)), by = .(col_f, tg_tick)]
right <- tt[diag == FALSE, .(m = median(med_core)), by = .(row_f, tr_tick)]
pb_top <- ggplot(top, aes(tg_tick, m)) +
  geom_segment(aes(xend = tg_tick, y = 0, yend = m), colour = GREYS[2], linewidth = 0.4) +
  geom_point(size = 1, colour = GREYS[3]) +
  facet_grid(. ~ col_f, scales = "free_x", space = "free_x") +
  scale_y_continuous(breaks = c(0, 0.2), limits = c(0, NA)) +
  labs(x = NULL, y = "test-set\nmedian r") +
  theme_cns + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                    strip.text = element_blank(), panel.spacing = unit(1.5, "pt"),
                    axis.title.y = element_text(size = 5), axis.text.y = element_text(size = 5))
pb_right <- ggplot(right, aes(m, tr_tick)) +
  geom_segment(aes(yend = tr_tick, x = 0, xend = m), colour = GREYS[2], linewidth = 0.4) +
  geom_point(size = 1, colour = GREYS[3]) +
  facet_grid(row_f ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(breaks = c(0, 0.2), limits = c(0, NA)) +
  labs(y = NULL, x = "model\nmedian r") +
  theme_cns + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    strip.text = element_blank(), panel.spacing = unit(1.5, "pt"),
                    axis.title.x = element_text(size = 5), axis.text.x = element_text(size = 5))
pb <- pb_main   # marginal lollipops (pb_top, pb_right) dropped on request, 2026-11-04

# ---- c: strategies per bulk cohort ----------------------------------------------
by <- fread(file.path(DATA, "fig9_by_bulk.tsv"))
cl <- melt(by, id.vars = c("bulk", "cls", "kind"), measure.vars = paste0(STRAT_LEV, "_med_core"),
           variable.name = "strategy", value.name = "r")[!is.na(r)]
cl[, strategy := factor(sub("_med_core", "", strategy), STRAT_LEV)]
cl[, bulk_lab := factor(cohort_label(bulk), cohort_label(by[order(cls, bulk)]$bulk))]
cl[, cls := factor(CLS_LAB[cls], CLS_LAB)]
pc <- ggplot(cl, aes(strategy, r, group = bulk_lab)) +
  geom_line(colour = GREYS[1], linewidth = 0.3) +
  geom_point(aes(colour = strategy, shape = kind), size = 1.6) +
  facet_wrap(~cls, nrow = 1) +
  scale_colour_manual(values = STRAT_COL, labels = STRAT_LAB, name = NULL) +
  scale_shape_manual(values = KIND_SHP, name = "Ground truth",
                     labels = c(cross = "other snRNA-seq", oof = "same snRNA-seq, out of fold")) +
  scale_x_discrete(labels = c(S1 = "S1", S2 = "S2", S3 = "S3", prod = "Prod")) +
  labs(x = NULL, y = "median r, Core Supertypes\n(per bulk cohort)") +
  theme_cns + theme(legend.position = "bottom", legend.box = "vertical")

# ---- d: per supertype, best strategy vs ICC --------------------------------------
summ <- fread(file.path(DATA, "fig9_summary.tsv"))
best <- cl[, .(m = median(r)), by = strategy][order(-m)]$strategy[1]
ct <- fread(file.path(DATA, "fig9_celltype_r.tsv"))
ci <- fread(file.path(HERE, "..", "common", "data", "reliability_ci_summary.tsv"))[level == "supertype"]
ci[, ct := gsub("[/ ]", "_", cell_type)]
d <- ct[strategy == best, .(r = median(r)), by = cell_type][ci, on = .(cell_type = ct), nomatch = 0]
d[, set := factor(ifelse(ci_lo > 0, "Core Supertype", "Non-core Supertype"),
                  c("Core Supertype", "Non-core Supertype"))]
sp <- spear(d$median_icc_sym, d$r)
pd <- ggplot(d, aes(median_icc_sym, r)) +
  geom_hline(yintercept = 0, colour = GREYS[1], linewidth = 0.3) +
  geom_point(aes(shape = set), colour = "grey30", size = 1.1) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 2.6,
           label = sprintf("rho = %.2f\np = %s, n = %d", sp$estimate, fmt_p(sp$p.value), nrow(d))) +
  labs(x = "median symmetric ICC",
       y = sprintf("CelMod r, %s\n(median over cohorts)", STRAT_LAB[[as.character(best)]])) +
  theme_cns + theme(legend.position = "bottom")

# a: draft schematic (fig9_schematic.R, user 2026-11-09); layout guide for BioRender
source(file.path(HERE, "fig9_schematic.R"))
pa <- make_schematic(STRAT_COL)
pc <- pc + guides(colour = guide_legend(nrow = 2, order = 1), shape = guide_legend(nrow = 2, order = 2)) +
  theme(legend.text = element_text(size = 6), legend.title = element_text(size = 6))
pd <- pd + guides(shape = guide_legend(nrow = 2)) +
  theme(legend.text = element_text(size = 6), legend.position = "bottom")
# Supplementary heatmap only (user 2026-11-04): no CV legend entry; the colour bar sits
# in the empty lower-left corner (under the row labels) instead of its own band.
# one text size everywhere on the supplementary heatmap (user 2026-11-04)
SUPP_PT <- 6.5
pbs <- pb + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = SUPP_PT),
                  axis.text.y = element_text(size = SUPP_PT),
                  strip.text = element_text(size = SUPP_PT),
                  strip.text.x.bottom = element_text(size = SUPP_PT),
                  strip.text.y.left = element_text(angle = 0, size = SUPP_PT),
                  axis.title = element_text(size = 8),   # axis titles larger (user 2026-11-04)
                  legend.title = element_text(size = SUPP_PT),
                  legend.text = element_text(size = SUPP_PT))
# supplementary only: drop the "ROSMAP: " prefix on the trio row (user 2026-11-04)
pbs <- pbs + scale_y_discrete(labels = function(x)
  sub("^ROSMAP: Green\\+Mathys\\+Diverse_multiome$", "Green+Mathys+Diverse_multiome", x))
supp_leg <- cowplot::get_legend(pbs + guides(shape = "none"))
supp <- pbs + guides(shape = "none", fill = "none") +
  inset_element(supp_leg, left = 0.02, bottom = 0.05, right = 0.33, top = 0.22,
                align_to = "full", on_top = TRUE)
save_fig(supp, file.path(FIG, "figure9_supp_matrix"), height_in = 5.9)
ggsave(file.path(FIG, "figure9_supp_matrix.svg"), supp, width = 183 / 25.4, height = 5.9,
       bg = "white", device = svglite::svglite)   # SVG copy (user 2026-11-04)

# ---- c (user 2026-11-09): supertype x test cohort map for the recommended model -------
# Spearman rho of the ROSMAP all-4 consensus model (Green+Mathys+Diverse_multiome+CUIMC2)
# for each Core Supertype in each test cohort (one bulk region per cohort, as panel b).
# Cohorts with two snRNA-seq ground truths: mean rho over the two. ROSMAP is 5-fold
# cross-validated (the model trained on these donors), as the circled cells of panel b.
S2_MODEL <- "ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2"
tax <- fread("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv")
tax[, cls := fcase(grepl("Glutamatergic", class_label), "Excitatory",
                   grepl("GABAergic", class_label), "Inhibitory", default = "Non-neuronal")]
tax[, ord := .I]
pt <- fread(file.path(DATA, "fig9_heatmap_spearman_celltype.tsv"))[trainer == S2_MODEL]
pt[, ub := sub("__.*", "", target)]
pt <- pt[ub %in% BULK_ORD, .(rho = mean(rho), n_sets = .N, cv = any(cv)), by = .(ub, cell_type)]
pt <- tax[, .(cell_type = supertype_id, label = supertype_label, cls, ord)][pt, on = "cell_type", nomatch = 0]
pt[, cls := factor(cls, c("Excitatory", "Inhibitory", "Non-neuronal"))]
pt[, st := factor(label, unique(pt[order(cls, ord)]$label))]
ROW_LAB <- setNames(ifelse(BULK_ORD == "ROSMAP", "ROSMAP (CV)", unname(SHORT[BULK_ORD])), BULK_ORD)
ROW_LAB <- sub("\n", " ", ROW_LAB)
pt[, coh := factor(ROW_LAB[ub], rev(ROW_LAB))]
pcm <- ggplot(pt, aes(st, coh, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(. ~ cls, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2F4FA8", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-0.2, 0.6), oob = scales::squish, name = "Spearman \u03c1") +
  labs(x = "Core Supertype", y = "Test cohort") +
  theme_cns + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
                    axis.text.y = element_text(size = SUPP_PT), strip.text = element_text(size = SUPP_PT),
                    axis.title = element_text(size = 8), panel.spacing = unit(2, "pt"),
                    legend.title = element_text(size = SUPP_PT), legend.text = element_text(size = SUPP_PT),
                    legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.5, "cm"),
                    axis.line = element_blank(), axis.ticks = element_blank())

# ---- compact strategy comparison (user 2026-11-09): main figure ----------------------
# One row per representative test cohort (bank order, as the heatmap); S1 / S2 / S3 as
# coloured dots on the common donor set of fig9_by_bulk.tsv (strategy_scores.py, Pearson).
# Open symbols = same-dataset out-of-fold ground truth (cohort has one snRNA-seq dataset).
# All 11 bulk regions + production trio + ICC link go to figure9_supp_strategies.
STRAT3 <- c("S1", "S2", "S3")
sc <- cl[bulk %in% BULK_ORD & strategy %in% STRAT3]
sc[, strategy := factor(as.character(strategy), STRAT3)]
sc[, coh := factor(sub("\n", " ", SHORT[bulk]), rev(sub("\n", " ", SHORT[BULK_ORD])))]
sc[, y := as.numeric(coh) + c(S1 = 0.2, S2 = 0, S3 = -0.2)[as.character(strategy)]]
rng <- sc[, .(lo = min(r), hi = max(r), y = as.numeric(coh[1])), by = coh]
pcs <- ggplot(sc) +
  geom_segment(data = rng, aes(x = lo, xend = hi, y = y, yend = y), colour = GREYS[1], linewidth = 0.4) +
  geom_point(aes(r, y, colour = strategy, shape = kind), size = 1.8, stroke = 0.6) +
  scale_y_continuous(breaks = seq_along(levels(sc$coh)), labels = levels(sc$coh),
                     expand = expansion(add = 0.4)) +
  scale_colour_manual(values = STRAT_COL[STRAT3], labels = STRAT_LAB[STRAT3], name = NULL) +
  scale_shape_manual(values = KIND_SHP, name = NULL,
                     labels = c(cross = "other snRNA-seq", oof = "same snRNA-seq (CV)")) +
  labs(x = "Median Pearson r (Core Supertypes)", y = "Test cohort") +
  guides(colour = guide_legend(order = 1, ncol = 1), shape = guide_legend(order = 2, ncol = 1)) +
  theme_cns + theme(axis.text = element_text(size = SUPP_PT), axis.title = element_text(size = 8),
                    legend.text = element_text(size = SUPP_PT), legend.position = "right",
                    legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(2, "pt"))

# Main figure (user 2026-11-09): a = schematic | b = strategy comparison, c = heatmap
# exactly as figure9_supp_matrix, d = supertype x cohort map for the recommended model.
pa <- pa + labs(tag = "a")
pbs_t <- pcs + labs(tag = "b")
pct <- wrap_elements(full = supp) + labs(tag = "c")
pdt <- wrap_elements(full = pcm) + labs(tag = "d")
top_row <- (pa + pbs_t + plot_layout(widths = c(1.4, 1))) &   # wider schematic (user 2026-11-09)
  theme(plot.tag = element_text(size = 10, face = "bold", colour = "black"))
fig <- (wrap_elements(full = top_row) / pct / pdt + plot_layout(heights = c(2.5, 5.9, 2.2))) &
  theme(plot.tag = element_text(size = 10, face = "bold", colour = "black"))   # one tag style
save_fig(fig, file.path(FIG, "figure9_celmod_strategies_draft"), height_in = 10.6)

# ---- supplementary strategy detail: all bulk regions, + production trio, ICC link ---------
supp_str <- (pc + labs(tag = "a")) + (pd + labs(tag = "b")) + plot_layout(widths = c(1.6, 1))
save_fig(supp_str, file.path(FIG, "figure9_supp_strategies"), height_in = 3.6)

message("best strategy for d: ", best)
print(summ)
