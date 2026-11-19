# Figure 9, version 2 (user 2026-11-10). Same data and message as fig9.R, compact visuals;
# fig9.R and its outputs (figure9_celmod_strategies_draft.*) are left untouched.
#   a  draft schematic (fig9_schematic.R)
#   b  S1 / S2 / S3 per test cohort (as fig9.R; Pearson, common donors, fig9_by_bulk.tsv)
#   c  which models generalise: one row per training model, one dot per test set
#      (median Spearman rho over Core Supertypes; fig9_heatmap_spearman.tsv, the data of
#      figure9_supp_matrix). Open = test dataset in the training label (5-fold CV);
#      bar = median over the filled (independent) test sets only.
#   d  which supertypes transfer: Spearman rho of the ROSMAP all-4 consensus for each of
#      the 86 Core Supertypes, per test cohort (fig9_heatmap_spearman_celltype.tsv; mean
#      over the two ground truths where a cohort has two; ROSMAP = 5-fold CV). Lines:
#      the 4 supertypes with rho > 0.2 in all six cohorts.
# Full grids stay in the supplement: figure9_supp_matrix (c), figure9_supp_supertype (d).
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig9_v2.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure9"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures")
PT <- 6.5
STRAT_COL <- c(S1 = "#7B8C9E", S2 = "#2F7E8F", S3 = "#E1A22E", prod = "black")
STRAT_LAB <- c(S1 = "S1 paired", S2 = "S2 global", S3 = "S3 class-matched")
KIND_SHP <- c(cross = 16, oof = 1)
BULK_ORD <- c("ROSMAP", "MSBB_FP", "Diverse_Rush2_DLPFC", "CommonMind_MSSM", "NIMH_HBCC", "NABEC")
BANK <- c(ROSMAP = "AMP-AD", MSBB_FP = "AMP-AD", Diverse_Rush2_DLPFC = "AMP-AD",
          CommonMind_MSSM = "PsychENCODE", NIMH_HBCC = "PsychENCODE", NABEC = "NABEC")
SHORT <- c(ROSMAP = "ROSMAP", MSBB_FP = "MSBB (FP)", Diverse_Rush2_DLPFC = "Diverse Rush",
           CommonMind_MSSM = "CMC", NIMH_HBCC = "HBCC", NABEC = "NABEC")
COH_COL <- setNames(c("#B2182B", "#E08214", "#C51B7D", "#2166AC", "#5AAED6", "#1B9E77"), BULK_ORD)
SN_ORD <- c("Green", "Mathys", "Multiome2025_DLPFC", "CUIMC2", "PsychAD_MSBB", "PsychAD_ROSMAP",
            "BrainSCOPE_CMC", "PsychAD_HBCC", "Catching_2026_HBCC", "Catching_2026_NABEC")
CONS_LAB <- c("ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2" = "ROSMAP: all 4 (S2)",
              "ROSMAP__Green+Mathys+Multiome2025_DLPFC" = "ROSMAP: Green+Mathys+Diverse_multiome",
              "NIMH_HBCC__PsychAD_HBCC+Catching_2026_HBCC" = "HBCC: PsychAD+Catching")
small <- theme(axis.text = element_text(size = PT), axis.title = element_text(size = 8),
               strip.text = element_text(size = PT), legend.text = element_text(size = PT),
               legend.title = element_text(size = PT), legend.key.size = unit(0.3, "cm"))

# ---- a: schematic --------------------------------------------------------------
source(file.path(HERE, "fig9_schematic.R"))
pa <- make_schematic(STRAT_COL, pt = 5)

# ---- b: strategies per test cohort (identical to fig9.R) -----------------------
by <- fread(file.path(DATA, "fig9_by_bulk.tsv"))
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
  labs(x = "Median Pearson r (Core Supertypes)", y = "Test cohort") +
  guides(colour = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1)) +
  theme_cns + small + theme(legend.position = "bottom", legend.box = "vertical",
                            legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))

# ---- c: which training models generalise -----------------------------------------
tt <- fread(file.path(DATA, "fig9_heatmap_spearman.tsv"))
tt[, `:=`(tb = sub("__.*", "", trainer), ts = sub(".*?__", "", trainer),
          ub = sub("__.*", "", target), us = sub(".*?__", "", target))]
tt[, cons := grepl("+", ts, fixed = TRUE)]
tt <- tt[ub %in% BULK_ORD & tb %in% BULK_ORD & (!cons | trainer %in% names(CONS_LAB))]
if (!"cv" %in% names(tt)) tt[, cv := FALSE]
tt[is.na(cv), cv := FALSE]
tt[, circled := diag | cv]
grp_lab <- function(b) sub("Diverse Rush", "Diverse", sub(" \\(FP\\)", "", SHORT[b]))
GRP_LEV <- c(grp_lab(BULK_ORD), "Consensus")
tt[, grp := factor(ifelse(cons, "Consensus", grp_lab(tb)), GRP_LEV)]
CONS_SHORT <- c("ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2" = "ROSMAP: all 4 (S2)",
                "ROSMAP__Green+Mathys+Multiome2025_DLPFC" = "ROSMAP: 3 (no CUIMC2)",
                "NIMH_HBCC__PsychAD_HBCC+Catching_2026_HBCC" = "HBCC: PsychAD+Catching")
tt[, model := ifelse(cons, CONS_SHORT[trainer], sub("Multiome2025_DLPFC", "Diverse_multiome", ts))]
mod_ord <- unique(tt[order(grp, cons, match(ts, SN_ORD), match(trainer, names(CONS_LAB)))]$model)
tt[, model_f := factor(model, rev(unique(mod_ord)))]
tt[, ub := factor(ub, BULK_ORD)]
med <- tt[circled == FALSE, .(m = median(med_core)), by = .(grp, model_f)]
pc <- ggplot(tt, aes(med_core, model_f)) +
  geom_vline(xintercept = 0, colour = GREYS[1], linewidth = 0.3) +
  geom_point(data = tt[circled == FALSE], aes(colour = ub), size = 1.5, alpha = 0.85,
             position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  geom_point(data = tt[circled == TRUE], aes(colour = ub), shape = 1, size = 1.5, stroke = 0.5) +
  geom_errorbarh(data = med, aes(xmin = m, xmax = m, y = model_f), inherit.aes = FALSE,
                 height = 0.75, linewidth = 0.7, colour = "black") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = COH_COL, labels = SHORT, name = "Test cohort") +
  labs(x = "Median Spearman ρ\n(Core Supertypes)", y = "Training model") +
  guides(colour = guide_legend(nrow = 3, title.position = "top", override.aes = list(size = 1.8, alpha = 1))) +
  theme_cns + small + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = PT),
                            strip.background = element_blank(), panel.spacing = unit(3, "pt"),
                            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
                            legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- d: which supertypes transfer (recommended model) ----------------------------
S2_MODEL <- "ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2"
tax <- fread("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv")
tax[, cls := fcase(grepl("Glutamatergic", class_label), "Excitatory",
                   grepl("GABAergic", class_label), "Inhibitory", default = "Non-neuronal")]
pt <- fread(file.path(DATA, "fig9_heatmap_spearman_celltype.tsv"))[trainer == S2_MODEL]
pt[, ub := sub("__.*", "", target)]
pt <- pt[ub %in% BULK_ORD, .(rho = mean(rho)), by = .(ub, cell_type)]
pt <- tax[, .(cell_type = supertype_id, label = supertype_label, cls)][pt, on = "cell_type", nomatch = 0]
COH_X <- setNames(c("ROSMAP\n(CV)", "MSBB", "Diverse", "CMC", "HBCC", "NABEC"), BULK_ORD)
pt[, coh := factor(COH_X[ub], COH_X)]
pt[, cls := factor(cls, c("Excitatory", "Inhibitory", "Non-neuronal"))]
robust <- pt[, .(ok = all(rho > 0.2) && .N == length(BULK_ORD)), by = label][ok == TRUE]$label
CLS_COL <- c(Excitatory = NEU, Inhibitory = "#8491B4", "Non-neuronal" = GLIA)
set.seed(1); pt[, xj := as.numeric(coh) + runif(.N, -0.17, 0.17)]   # one jitter for dots + rings
rl <- pt[label %in% robust]
pdd <- ggplot(pt, aes(coh, rho)) +
  geom_hline(yintercept = 0, colour = GREYS[1], linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, colour = "grey35", fill = NA) +
  geom_point(aes(xj, rho, colour = cls), size = 0.9, alpha = 0.75) +
  geom_point(data = rl, aes(xj, rho, shape = "ρ > 0.2 in all six cohorts"), size = 1.7, stroke = 0.5,
             colour = "black") +
  scale_shape_manual(values = 1, name = NULL) +
  scale_colour_manual(values = CLS_COL, name = NULL) +
  labs(x = "Test cohort", y = "Spearman ρ per Core Supertype\n(S2 model: ROSMAP all 4)") +
  guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.8, alpha = 1)),
         shape = guide_legend(order = 2)) +
  theme_cns + small + theme(legend.position = "bottom", legend.box = "vertical",
                            legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))

# ---- assemble --------------------------------------------------------------------
TAG <- theme(plot.tag = element_text(size = 10, face = "bold", colour = "black"))
top_row <- (pa + labs(tag = "a") + wrap_elements(full = pb) + labs(tag = "b") +
              plot_layout(widths = c(1.5, 1))) & TAG
bot_row <- (wrap_elements(full = pc) + labs(tag = "c") + wrap_elements(full = pdd) + labs(tag = "d") +
              plot_layout(widths = c(1.25, 1))) & TAG
fig <- (wrap_elements(full = top_row) / bot_row + plot_layout(heights = c(2.5, 3.2))) & TAG
save_fig(fig, file.path(FIG, "figure9_celmod_strategies_v2"), height_in = 5.9)
ggsave(file.path(FIG, "figure9_celmod_strategies_v2.svg"), fig, width = 183 / 25.4, height = 5.9,
       bg = "white", device = svglite::svglite)

# ---- supplement: the full supertype x cohort grid (former panel d) -----------------
tax[, ord := .I]
ps <- tax[, .(label = supertype_label, ord)][pt, on = "label"]
ps[, st := factor(label, unique(ps[order(cls, ord)]$label))]
ps[, cohr := factor(sub("\n", " ", as.character(coh)), rev(sub("\n", " ", COH_X)))]
psup <- ggplot(ps, aes(st, cohr, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(. ~ cls, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2F4FA8", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-0.2, 0.6), oob = scales::squish, name = "Spearman ρ") +
  labs(x = "Core Supertype", y = "Test cohort") +
  theme_cns + small + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
                            axis.line = element_blank(), axis.ticks = element_blank(),
                            panel.spacing = unit(2, "pt"))
save_fig(psup, file.path(FIG, "figure9_supp_supertype"), height_in = 2.6)
message("robust supertypes (rho > 0.2 in all six): ", paste(robust, collapse = ", "))
