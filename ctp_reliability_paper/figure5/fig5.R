# Figure 5 (DRAFT). Transferring CelMod across cohorts.
# Message: a model trained on ROSMAP keeps most or all of its accuracy for Core
# Supertypes when applied to another cohort's bulk (MSBB), and where it does not
# (NIMH HBCC) every training cohort loses about as much, so the loss is the target's.
# Grouping is the Fig 4 Core Supertype set (reliability-only: lower bound of the
# donor-bootstrap 95% CI of the median symmetric ICC > 0, no bulk data used), not the
# PI's CORE portability flag, which is partly selected on cross-cohort transfer and
# would make this figure circular (user decision 2026-10-28).
# Transfer is scored only where the target has its own matched snRNA-seq: ROSMAP,
# MSBB (IFG, STG) and NIMH HBCC (the PI's 192 design). In-sample pairs (a model
# applied to its own training donors) are excluded, as in 192_retention_clean.csv.
#   a  empty slot for the BioRender design schematic (drawn by hand)
#   b  ROSMAP-trained model vs the target's own same-cohort model, per supertype r
#      against the target's snRNA labels, in MSBB IFG, MSBB STG, NIMH HBCC
#   c  summary of b: median r, same-cohort vs ROSMAP-trained, Core vs Non-core, with
#      Core retention (median r transferred / median r same-cohort)
#   d  Core retention for every train -> target pair, bootstrap 95% CI over supertypes
# Supplementary figures:
#   figure5_supp_allpairs   all-pairs median r heatmap; per-supertype cost of transfer
#   figure5_supp_controls   (PI's 34 CORE supertypes, scores precomputed on that set)
#                           reference-profile NNLS vs CelMod; shared-trainer inflation
# Inputs from fig5_data.py. Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig5.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure5"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures"); dir.create(FIG, showWarnings = FALSE)

# Cohort names are the paper's display labels (common/cohort_names.tsv via cohort_label());
# the pooled "MSBB + HBCC" model is not a cohort and passes through unchanged.
# Styling follows Fig 4b: Core vs Non-core Supertype is shape only, filled vs open circle
# in grey30 (fig4.R ICC_POS_SHP), so circles are reserved for that. Training cohorts get
# other shapes and colours. Colours come from the palettes already used in earlier
# figures, not theme_cns ACCENT (used by no earlier figure): training cohorts take the
# muted categorical set of Figs 3/8 (PAIR_COL) with the pooled model in black (as Fig 8's
# pooled estimate); the target-trained reference model takes Fig 3's neutral slate
# (MODEL_COL). ROSMAP-trained = ROSMAP's colour, diamond, in both c and d. The cell-class
# colours (Figs 1/2/4) and Nature blue (NEU, = neuronal) are not used here.
TRAIN_COL <- c(ROSMAP = "#2F7E8F", MSBB = "#E1A22E", HBCC = "#8B6FA8", `MSBB + HBCC` = "black")
TARGET_COL <- "#7B8C9E"
SET_LEV <- c("Core Supertype", "Non-core Supertype")
SET_SHP <- setNames(c(16, 1), SET_LEV)                     # as fig4.R ICC_POS_SHP
SET_FILL <- setNames(c(GREYS[3], "white"), SET_LEV)       # box fill, supp gap panel
SRC_LEV <- c("ROSMAP", "MSBB", "HBCC", "MSBB + HBCC"); TGT_LEV <- c("ROSMAP", "MSBB_IFG", "MSBB_STG", "HBCC")
COHORT_OF <- c(ROSMAP = "ROSMAP", MSBB_IFG = "MSBB", MSBB_STG = "MSBB", HBCC = "HBCC")
relabel <- function(d) d[, `:=`(source_lab = cohort_label(source_lab), target_lab = cohort_label(target_lab))][]

m <- relabel(fread(file.path(DATA, "fig5b_retention_by_type.tsv"))[!is.na(core_supertype)])
m[, set := factor(ifelse(core_supertype, SET_LEV[1], SET_LEV[2]), levels = SET_LEV)]
pairs <- unique(fread(file.path(DATA, "fig5_retention_by_set.tsv"))[, .(source, target)])  # out-of-sample pairs only
mo <- merge(m, pairs, by = c("source", "target"))[!is.na(r_transfer) & !is.na(r_same)]

# ---- b: ROSMAP-trained vs same-cohort, per supertype ----------------------
XT_LEV <- c("MSBB_IFG", "MSBB_STG", "HBCC")
b <- mo[source == "srcROSMAP" & target_lab %in% XT_LEV]
b[, target_lab := factor(target_lab, levels = XT_LEV)]
print(b[, .(med_same = median(r_same), med_transfer = median(r_transfer), n = .N), by = .(target_lab, set)])
p_b <- ggplot(b, aes(r_same, r_transfer)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(aes(shape = set), colour = "grey30", size = 1.1, alpha = 0.8, stroke = 0.35) +
  facet_wrap(~ target_lab, nrow = 1) +
  scale_shape_manual(values = SET_SHP, breaks = SET_LEV, name = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6)) + scale_y_continuous(breaks = c(-0.3, 0, 0.3, 0.6)) +
  labs(x = "Agreement with snRNA-seq (r), target-trained model",
       y = "Agreement with snRNA-seq (r),\nROSMAP-trained model") +
  coord_equal() +
  theme_cns + theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(0, 1),
                    legend.background = element_blank(), legend.key.size = unit(8, "pt"))

# ---- c: summary of b -------------------------------------------------------
cs <- relabel(fread(file.path(DATA, "fig5_retention_by_set.tsv")))[source == "srcROSMAP" & target_lab %in% XT_LEV]
cs[, `:=`(target_lab = factor(target_lab, levels = XT_LEV), scope = factor(scope, levels = rev(SET_LEV)))]
csl <- melt(cs, id.vars = c("target_lab", "scope", "retention"), measure.vars = c("med_same", "med_transfer"),
            variable.name = "model", value.name = "r")
csl[, model := factor(model, levels = c("med_same", "med_transfer"),
                      labels = c("Target-trained model", "ROSMAP-trained model"))]
p_cs <- ggplot(cs, aes(y = scope)) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_segment(aes(x = med_same, xend = med_transfer), colour = "grey60", linewidth = 0.5) +
  geom_point(data = csl, aes(x = r, shape = model, fill = model), colour = "black", size = 2.2, stroke = 0.4) +
  geom_text(data = cs[scope == SET_LEV[1]], aes(x = pmax(med_same, med_transfer), label = sprintf("%.0f%%", 100 * retention)),
            hjust = -0.35, size = fontsize / .pt) +
  facet_wrap(~ target_lab, ncol = 1) +
  scale_shape_manual(values = c(22, 23), name = NULL) +           # target-trained = square, ROSMAP-trained = diamond (as in d)
  scale_fill_manual(values = c(TARGET_COL, TRAIN_COL[["ROSMAP"]]), name = NULL) +
  scale_y_discrete(labels = function(x) sub(" Supertype", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Median agreement with snRNA-seq (r)", y = NULL) +
  theme_cns + theme(legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 0))

# ---- d: Core retention, every pair -----------------------------------------
dr <- relabel(fread(file.path(DATA, "fig5_core_retention_ci.tsv")))
dr[, `:=`(target_lab = factor(target_lab, levels = rev(TGT_LEV)),
          source_lab = factor(source_lab, levels = SRC_LEV))]
print(dr[, .(source_lab, target_lab, retention, ci_lo, ci_hi)], digits = 2)
SRC_COL <- TRAIN_COL[SRC_LEV]   # CI lines and points
SRC_SHP <- setNames(c(23, 22, 24, 25), SRC_LEV)  # no circles: those mean Core / Non-core
pd_ <- position_dodge(width = 0.6)
p_dr <- ggplot(dr, aes(retention, target_lab, group = factor(source_lab, levels = rev(SRC_LEV)))) +  # ROSMAP on top within each row
  geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, colour = source_lab), height = 0, linewidth = 0.4, position = pd_) +
  geom_point(aes(fill = source_lab, shape = source_lab), colour = "black", size = 1.8, stroke = 0.3, position = pd_) +
  scale_colour_manual(values = SRC_COL, name = "Trained on") +
  scale_fill_manual(values = SRC_COL, name = "Trained on") +
  scale_shape_manual(values = SRC_SHP, name = "Trained on") +
  labs(x = "Core Supertype retention", y = "Applied to") +
  theme_cns + theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0))

# ---- supplementary: all pairs, median r ------------------------------------------------
off <- mo[, .(r = median(r_transfer), n = .N), by = .(source_lab, target_lab, set)][, diag := FALSE]
same <- unique(m[, .(target_lab, cell_type, set, r_same)])[!is.na(r_same)][
  , .(r = median(r_same), n = .N), by = .(target_lab, set)][, `:=`(source_lab = COHORT_OF[target_lab], diag = TRUE)]
cc <- rbind(off, same, use.names = TRUE)
cc[, `:=`(source_lab = factor(source_lab, levels = rev(SRC_LEV)), target_lab = factor(target_lab, levels = TGT_LEV))]
print(dcast(cc, set + source_lab ~ target_lab, value.var = "r"), digits = 2)
grid <- CJ(source_lab = factor(rev(SRC_LEV), levels = rev(SRC_LEV)), target_lab = factor(TGT_LEV, levels = TGT_LEV),
           set = factor(SET_LEV, levels = SET_LEV))
p_c <- ggplot(cc, aes(target_lab, source_lab)) +
  geom_tile(data = grid, fill = "grey95", colour = "white", linewidth = 0.6) +
  geom_tile(aes(fill = r), colour = "white", linewidth = 0.6) +
  geom_tile(data = cc[diag == TRUE], fill = NA, colour = "black", linewidth = 0.5, width = 0.94, height = 0.94) +
  geom_text(aes(label = sprintf("%.2f", r), fontface = ifelse(diag, "italic", "plain")),
            size = fontsize / .pt, colour = "black") +
  facet_wrap(~ set, nrow = 1) +
  scale_fill_gradient(low = "#F7F7F7", high = GREYS[2], limits = c(0, 0.3), oob = scales::squish,
                      name = "Median r", guide = guide_colourbar(title.position = "left", title.vjust = 1)) +
  scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
  scale_y_discrete(labels = function(x) sub(" \\+ ", " +\n", x)) +
  labs(x = "Applied to (bulk cohort)", y = "Trained on") +
  coord_fixed() +
  theme_cns + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    legend.position = "bottom", legend.key.height = unit(6, "pt"),
                    legend.key.width = unit(30, "pt"), legend.margin = margin(0, 0, 0, 0),
                    plot.margin = margin(2, 4, 2, 4))

# ---- supplementary: cost of transfer ---------------------------------------------------
PAIR_LEV <- c("ROSMAP>MSBB_IFG", "ROSMAP>MSBB_STG", "ROSMAP>HBCC", "MSBB>ROSMAP", "MSBB>HBCC",
              "HBCC>ROSMAP", "HBCC>MSBB_IFG", "HBCC>MSBB_STG", "MSBB + HBCC>ROSMAP")
mo[, gap := r_transfer - r_same]
mo[, pair := factor(paste(source_lab, target_lab, sep = ">"), levels = PAIR_LEV)]
stopifnot(!anyNA(mo$pair))
pv <- mo[, .(p = wilcox.test(gap[set == SET_LEV[1]], gap[set == SET_LEV[2]])$p.value,
             med_core = median(gap[set == SET_LEV[1]]), med_non = median(gap[set == SET_LEV[2]])), by = pair]
print(pv[order(pair)], digits = 2)
pv[, lab := ifelse(p < 0.001, "p < 0.001", sprintf("p = %.2g", p))]
p_d <- ggplot(mo, aes(pair, gap, fill = set)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_boxplot(width = 0.7, outlier.size = 0.4, outlier.stroke = 0.2, linewidth = 0.3,
               position = position_dodge(width = 0.8)) +
  geom_text(data = pv, aes(pair, max(mo$gap) + 0.1, label = lab), inherit.aes = FALSE, size = fontsize / .pt) +
  scale_fill_manual(values = SET_FILL, name = NULL) +
  scale_x_discrete(labels = function(x) sub(">", "\n→\n", x)) +
  labs(x = "Trained on → applied to", y = "Cost of transfer per supertype\n(r transferred − r same-cohort)") +
  theme_cns + theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

# ---- supplementary: reference profiles vs CelMod ----------------------------
s1 <- fread(file.path(DATA, "fig5c_profile_vs_celmod.tsv"))
s1[, arm_lab := fcase(arm == "CelMod (Multiome-trained)", "CelMod\n(trained on matched pairs)",
                      arm == "profile_hier", "Reference profiles\n(hierarchical NNLS)",
                      arm == "profile_flat", "Reference profiles\n(flat NNLS)")]
cl <- melt(s1, id.vars = "arm_lab", measure.vars = c("disatt_green_core", "h1_rel_core"),
           variable.name = "metric", value.name = "value")
cl[, metric := factor(metric, levels = c("disatt_green_core", "h1_rel_core"),
                      labels = c("Disattenuated r\nvs Green", "Triangulation\nreliability"))]
cl[, arm_lab := factor(arm_lab, levels = c("Reference profiles\n(flat NNLS)", "Reference profiles\n(hierarchical NNLS)", "CelMod\n(trained on matched pairs)"))]
p_s1 <- ggplot(cl, aes(value, arm_lab, fill = grepl("CelMod", arm_lab))) +
  geom_col(width = 0.6, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.15, size = fontsize / .pt) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = TRAIN_COL[["ROSMAP"]], `FALSE` = "grey85"), guide = "none") +  # CelMod here is ROSMAP-trained
  scale_x_continuous(expand = expansion(mult = c(0, 0.45)), n.breaks = 3) +
  labs(x = "Score on the 34 CORE supertypes", y = NULL) +
  theme_cns

# ---- supplementary: shared-trainer inflation --------------------------------
s2 <- fread(file.path(DATA, "fig5d_trainer_inflation.tsv"))
cat(sprintf("supp: median inflation all %+.3f, CORE %+.3f\n", median(s2$inflation), median(s2[core == TRUE]$inflation)))
p_s2 <- ggplot(s2, aes(inflation, fill = core)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_histogram(binwidth = 0.05, colour = "white", linewidth = 0.2, boundary = 0) +
  geom_vline(xintercept = median(s2$inflation), colour = "black", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = fontsize / .pt,
           label = sprintf("median %+.3f\nn = %d supertypes", median(s2$inflation), nrow(s2))) +
  scale_fill_manual(values = c(`TRUE` = GREYS[3], `FALSE` = GREYS[1]),  # PI CORE set
                    labels = c(`TRUE` = "CORE", `FALSE` = "other"), name = NULL) +
  labs(x = "Agreement of two bulk legs: shared training cohort − different cohorts (Δr)",
       y = "Supertypes") +
  theme_cns + theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(0, 1))

# ---- assemble --------------------------------------------------------------
# a: empty slot (reserved for the BioRender schematic); a blank ggplot, unlike plot_spacer(), still gets a panel letter
blank <- ggplot() + theme_void() + theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))
design <- "
AAA
BBB
CDD
"
fig <- wrap_plots(A = blank, B = free(p_b, type = "panel"), C = free(p_cs, type = "panel"),
                  D = free(p_dr, type = "panel"), design = design) +
  plot_layout(heights = c(0.9, 1, 1.1)) + plot_annotation(tag_levels = "a")
save_fig(fig, file.path(FIG, "figure5_transfer"), height_in = 8)

design_a <- "
A
B
"
sup_a <- wrap_plots(A = p_c, B = p_d, design = design_a) +
  plot_layout(heights = c(1.1, 0.9)) + plot_annotation(tag_levels = "a")
save_fig(sup_a, file.path(FIG, "figure5_supp_allpairs"), height_in = 6.5)

design_s <- "
AB
CC
"
sup <- wrap_plots(A = free(p_s1, type = "panel"), B = plot_spacer(),  # untagged: an empty slot, not a panel
                   C = free(p_s2, type = "panel"), design = design_s) +
  plot_layout(heights = c(1, 0.8)) + plot_annotation(tag_levels = "a")
save_fig(sup, file.path(FIG, "figure5_supp_controls"), height_in = 5.5)
