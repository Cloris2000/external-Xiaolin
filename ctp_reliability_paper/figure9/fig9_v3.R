# Figure 9, version 3 (user 2026-11-10): one summary number per mark, so each panel makes
# a single point. fig9.R (draft) and fig9_v2.R are left untouched.
#   a  draft schematic (fig9_schematic.R)
#   b  S1 / S2 / S3 per test cohort (as fig9.R; Pearson, common donors, fig9_by_bulk.tsv)
#   c  which models generalise: each training model's median Spearman rho over its
#      INDEPENDENT test sets (fig9_heatmap_spearman.tsv; cells whose test dataset is part of
#      the training label -- diagonal and consensus CV cells -- excluded), ranked
#   d  which cell classes transfer: ROSMAP all-4 consensus, median Spearman rho over the
#      Core Supertypes of each class, per test cohort (fig9_heatmap_spearman_celltype.tsv;
#      mean over the two ground truths where a cohort has two; ROSMAP = 5-fold CV)
# Full grids: figure9_supp_matrix (c), figure9_supp_supertype (d; written by fig9_v2.R).
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig9_v3.R
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

# ---- c: which training models generalise (one number per model) -------------------
tt <- fread(file.path(DATA, "fig9_heatmap_spearman.tsv"))
tt[, `:=`(tb = sub("__.*", "", trainer), ts = sub(".*?__", "", trainer), ub = sub("__.*", "", target))]
tt[, cons := grepl("+", ts, fixed = TRUE)]
if (!"cv" %in% names(tt)) tt[, cv := FALSE]
tt[is.na(cv), cv := FALSE]
tt <- tt[ub %in% BULK_ORD & tb %in% BULK_ORD & (!cons | trainer %in% names(CONS_LAB)) & !(diag | cv)]
SN_SHORT <- function(s) sub("Multiome2025_DLPFC", "Diverse_multiome", s)
mc <- tt[, .(m = median(med_core), n = .N), by = .(trainer, tb, ts, cons)]
mc[, lab := ifelse(cons, CONS_LAB[trainer], paste0(sub(" \\(FP\\)", "", SHORT[tb]), ": ", SN_SHORT(ts)))]
mc[, lab := sub("ROSMAP: Green\\+Mathys\\+Diverse_multiome", "ROSMAP: Green+Mathys+Diverse_mult.", lab)]
mc[, lab_f := factor(lab, mc[order(m)]$lab)]
mc[, kind := factor(ifelse(cons, "Consensus (multi-dataset)", "Single dataset"),
                    c("Consensus (multi-dataset)", "Single dataset"))]
pc <- ggplot(mc, aes(m, lab_f, colour = kind)) +
  geom_segment(aes(x = 0, xend = m, yend = lab_f), linewidth = 0.5) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%.2f", m)), hjust = -0.35, size = PT / .pt, show.legend = FALSE) +
  scale_colour_manual(values = c(STRAT_COL[["S2"]], GREYS[1]), name = NULL) +
  scale_x_continuous(limits = c(0, 0.235), breaks = c(0, 0.1, 0.2), expand = expansion(mult = c(0, 0))) +
  labs(x = "Median Spearman ρ\n(independent test sets)", y = "Training model") +
  theme_cns + small + theme(legend.position = "none")   # teal = multi-dataset consensus (stated in the legend text)

# ---- d: which cell classes transfer (recommended model) ---------------------------
S2_MODEL <- "ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2"
tax <- fread("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv")
tax[, cls := fcase(grepl("Glutamatergic", class_label), "Excitatory",
                   grepl("GABAergic", class_label), "Inhibitory", default = "Non-neuronal")]
pt <- fread(file.path(DATA, "fig9_heatmap_spearman_celltype.tsv"))[trainer == S2_MODEL]
pt[, ub := sub("__.*", "", target)]
pt <- pt[ub %in% BULK_ORD, .(rho = mean(rho)), by = .(ub, cell_type)]
pt <- tax[, .(cell_type = supertype_id, cls)][pt, on = "cell_type", nomatch = 0]
cm <- pt[, .(m = median(rho), n = .N), by = .(ub, cls)]
COH_X <- setNames(c("ROSMAP\n(CV)", "MSBB", "Diverse", "CMC", "HBCC", "NABEC"), BULK_ORD)
cm[, x := match(ub, BULK_ORD)]
CLS_COL <- c(Excitatory = NEU, Inhibitory = "#8491B4", "Non-neuronal" = GLIA)
ends <- cm[x == length(BULK_ORD)]
ends[, ylab := m]
ends[cls == "Excitatory", ylab := m - 0.012]; ends[cls == "Inhibitory", ylab := m + 0.012]
pd <- ggplot(cm, aes(x, m, colour = cls)) +
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "#EEE8F4", alpha = 0.6) +
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "#E4F1EC", alpha = 0.6) +
  annotate("text", x = 2, y = 0.335, label = "Aging / AD", size = PT / .pt, colour = "grey30") +
  annotate("text", x = 5, y = 0.335, label = "Adult / psychiatric", size = PT / .pt, colour = "grey30") +
  geom_hline(yintercept = 0, colour = GREYS[1], linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  geom_text(data = ends, aes(x = 6.6, y = ylab, label = cls), hjust = 0, size = PT / .pt) +
  scale_colour_manual(values = CLS_COL, guide = "none") +
  scale_x_continuous(breaks = seq_along(BULK_ORD), labels = COH_X, expand = expansion(add = 0)) +
  scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.3, 0.1)) +
  coord_cartesian(xlim = c(0.5, 6.5), clip = "off") +
  labs(x = "Test cohort", y = "Median Spearman ρ\n(Core Supertypes in class)") +
  theme_cns + small + theme(plot.margin = margin(5.5, 48, 5.5, 5.5))

# ---- assemble --------------------------------------------------------------------
TAG <- theme(plot.tag = element_text(size = 10, face = "bold", colour = "black"))
top_row <- (pa + labs(tag = "a") + wrap_elements(full = pb) + labs(tag = "b") +
              plot_layout(widths = c(1.5, 1))) & TAG
bot_row <- (wrap_elements(full = pc) + labs(tag = "c") + wrap_elements(full = pd) + labs(tag = "d") +
              plot_layout(widths = c(1.1, 1))) & TAG
fig <- (wrap_elements(full = top_row) / bot_row + plot_layout(heights = c(2.5, 2.9))) & TAG
save_fig(fig, file.path(FIG, "figure9_celmod_strategies_v3"), height_in = 5.6)
ggsave(file.path(FIG, "figure9_celmod_strategies_v3.svg"), fig, width = 183 / 25.4, height = 5.6,
       bg = "white", device = svglite::svglite)
print(mc[order(-m), .(lab, m = round(m, 3), n)])
print(dcast(cm, ub ~ cls, value.var = "m")[order(match(ub, BULK_ORD))])
