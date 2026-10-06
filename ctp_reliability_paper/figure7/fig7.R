# Figure 7. Cross-region reliability of cell-type proportions, same donors,
# same dataset, two cortical regions instead of two datasets. Companion to
# Figure 1 (same-region, cross-dataset reliability): same three lineages
# (Immune = the DFC 2026 subclass holding microglia, Sst, L2/3 IT), same
# metric, same CLR/pseudocount convention.
#   (a) one panel per lineage: symmetric ICC of the subclass and of each of
#       its supertypes, both cohorts overlaid at the same x position (colour)
#       so the same cell type's reliability in Multiome2025 (DLPFC vs STG)
#       and in SEA-AD Multiregion (DFC vs MTG) can be compared directly.
#       Supertypes are ordered by their SEA-AD Multiregion ICC (native calls,
#       same donors measured on one pipeline -- the cleaner of the two
#       cohorts); the one supertype SEA-AD does not have (L2/3 IT_14-SEAAD,
#       Multiome2025-only) is placed last in its facet, and Monocyte
#       (SEA-AD-only) has no Multiome2025 point. Two cohorts:
#         Multiome2025   AMP-AD Diverse multiome 2025, DLPFC vs STG, same 398
#                        donors (BRISC-labelled, bridged to registry TLP ids)
#         SEA-AD Multiregion   DFC vs MTG, same 80 donors, native SEA-AD
#                        supertype calls (no label transfer needed)
#       One pair per cohort (not five independent pairs as in Fig 1), so
#       there is no across-pair IQR to show (fig7_data.py still writes a
#       donor-bootstrap q25/q75 per point; not plotted here since a
#       single-pair CI is not the same quantity as Fig 1's IQR).
#   (b) all ~150 DFC-2026 supertypes (not just the three target lineages),
#       ranked by SEA-AD Multiregion icc_sym, both cohorts overlaid by
#       colour; only the most and least reliable few are labelled (as Fig 1
#       panel a). Point estimate only, no bootstrap (fig7_data.py skips it
#       for this pass to keep runtime down across ~150 types x 2 cohorts).
# Metric: symmetric ICC = 2 cov / (var_a + var_b) of matched donors' CLR
# values between the two regions, identical formula to Figure 1 and
# common/reliability_stats.py. Taxonomy: SEA-AD DFC 2026 for the target
# lineages; a handful of supertypes found only in other SEA-AD regions
# (hippocampal, entorhinal, a few rare cortical states) are not in that
# taxonomy and are pooled into "Other" for the subclass-level CLR denominator
# only -- none belong to Immune/Sst/L2/3 IT (see fig7_data.py).
# Inputs (from fig7_data.py):
#   data/fig7a_pairs.tsv   per cohort x level x cell type: icc_sym, bootstrap
#                          median/q25/q75
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig7.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(ggrepel)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure7"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)
TAX <- "/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv"

SUBCLASSES <- c("Immune", "Sst", "L2/3 IT")
COHORTS <- c("Multiome2025", "SEA-AD Multiregion")
COHORT_LAB <- c("Multiome2025" = "Multiome RADC (DLPFC vs STG)",
                "SEA-AD Multiregion" = "SEA-AD Multiregion (DFC vs MTG)")
COHORT_COL <- c("Multiome2025" = "#B5651D", "SEA-AD Multiregion" = "#3C5488")

LEVEL_SHP <- c("Subclass" = 22, "Supertype" = 21)
LEVEL_SIZE <- c("Subclass" = 2.3, "Supertype" = 1.4)

# same Allen-style class hues as Figure 1 (user decision 2026-10-05); used
# here only on the taxonomy side (class_of), not for point colour, since
# colour now encodes cohort
CLS_COL <- c("Excitatory" = "#D1491C", "Inhibitory" = "#2F4FA8", "Non-neuronal" = "#2E8B3C")

tax <- fread(TAX)
sc_class <- unique(tax[, .(subclass = subclass_label, class_id)])
stopifnot(!anyDuplicated(sc_class$subclass))
tax <- tax[, .(cell_type = supertype_label, class_id, reactive_state)]
class_of <- function(class_id) {
  factor(fifelse(class_id == "Glia", "Non-neuronal", class_id), levels = names(CLS_COL))
}

# ---- (a) per-cohort cross-region ICC, subclass vs supertypes ---------------
a <- fread(file.path(HERE, "data/fig7a_pairs.tsv"))
a <- merge(a, tax, by = "cell_type", all.x = TRUE)
sc_lookup <- setNames(sc_class$class_id, sc_class$subclass)
a[level == "subclass", `:=`(class_id = unname(sc_lookup[as.character(subclass)]),
                            reactive_state = FALSE)]
stopifnot(!anyNA(a$class_id))
a[, cls := class_of(class_id)]
stopifnot(!anyNA(a$cls))
a[, lvl := factor(fifelse(level == "subclass", "Subclass", "Supertype"),
                  levels = names(LEVEL_SHP))]
a[, subclass := factor(subclass, levels = SUBCLASSES)]
a[, cohort := factor(cohort, levels = COHORTS)]

sub_est <- a[level == "subclass", .(cohort, subclass, sub_icc = icc_sym)]
bold_subclass <- function(x)
  as.expression(lapply(x, function(s) if (s %in% SUBCLASSES) bquote(bold(.(s))) else s))

# x order within each subclass facet: by SEA-AD Multiregion icc_sym (the
# native-label, same-pipeline cohort); supertypes SEA-AD does not have
# (L2/3 IT_14-SEAAD, Multiome2025-only) are placed last in their facet
ord <- a[cohort == "SEA-AD Multiregion", .(cell_type, subclass, level, ord_icc = icc_sym)]
ord <- merge(unique(a[, .(cell_type, subclass, level)]), ord,
            by = c("cell_type", "subclass", "level"), all.x = TRUE)
ord[is.na(ord_icc), ord_icc := -Inf]
setorder(ord, subclass, level, -ord_icc)
a[, cell_type := factor(cell_type, levels = unique(ord$cell_type))]

dot_panel <- function(legend = FALSE) {
  ggplot(a, aes(cell_type, icc_sym)) +
    geom_hline(data = sub_est, aes(yintercept = sub_icc, colour = cohort),
               linewidth = 0.3, linetype = "dashed", show.legend = FALSE) +
    geom_point(aes(colour = cohort, shape = lvl, size = lvl),
              position = position_dodge(width = 0.5), stroke = 0.6,
              fill = "white") +
    scale_colour_manual(values = COHORT_COL, labels = COHORT_LAB, name = NULL) +
    scale_shape_manual(values = LEVEL_SHP, name = NULL) +
    scale_size_manual(values = LEVEL_SIZE, name = NULL) +
    scale_x_discrete(labels = bold_subclass) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    facet_wrap(~ subclass, nrow = 1, scales = "free_x") +
    labs(x = NULL, y = "Cross-region intraclass\ncorrelation (ICC)") +
    guides(size = "none",
           shape = guide_legend(nrow = 1, override.aes = list(
             size = unname(LEVEL_SIZE), fill = "white", colour = "grey30")),
           colour = guide_legend(nrow = 2, override.aes = list(
             shape = 16, size = 2.2))) +
    theme_cns +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = if (legend) "bottom" else "none",
          legend.justification = "center",
          legend.margin = margin(0, 0, 0, 0),
          legend.box = "vertical")
}
p_a <- dot_panel(legend = TRUE)

cat("\nsubclass icc_sym by cohort:\n")
print(sub_est[order(cohort, subclass)])

# ---- (b) all supertypes, ranked by SEA-AD Multiregion ICC ------------------
N_LABEL <- 8   # label the top/bottom N_LABEL by SEA-AD Multiregion ICC (the
               # ranking axis), one label per cell type on its SEA-AD point
b <- fread(file.path(HERE, "data/fig7_all_supertypes.tsv"))
b[, cohort := factor(cohort, levels = COHORTS)]

ord_b <- b[cohort == "SEA-AD Multiregion", .(cell_type, ord_icc = icc_sym)]
ord_b <- merge(unique(b[, .(cell_type)]), ord_b, by = "cell_type", all.x = TRUE)
ord_b[is.na(ord_icc), ord_icc := -Inf]
setorder(ord_b, -ord_icc)
ord_b[, rank := .I]
b <- merge(b, ord_b[, .(cell_type, rank)], by = "cell_type")

# label the most/least reliable supertypes by the ranking axis (SEA-AD
# Multiregion); one label per type, placed at its SEA-AD point so it is not
# duplicated at the (different) Multiome2025 value of the same cell type
n_b <- nrow(ord_b)
lab_top <- b[cohort == "SEA-AD Multiregion"][order(-icc_sym)][1:N_LABEL]
lab_bot <- b[cohort == "SEA-AD Multiregion"][order(icc_sym)][1:N_LABEL]

p_b <- ggplot(b, aes(rank, icc_sym)) +
  geom_point(aes(colour = cohort), size = 0.9, alpha = 0.75, stroke = 0) +
  geom_text_repel(data = lab_top, aes(label = cell_type),
                  colour = "black", size = fontsize / ggplot2::.pt, seed = 1,
                  min.segment.length = 0, segment.size = 0.22, segment.colour = "grey60",
                  box.padding = 0.15, point.padding = 0.15, max.overlaps = Inf,
                  hjust = 0, direction = "y", nudge_x = n_b * 0.18 - lab_top$rank,
                  nudge_y = 0.85 - lab_top$icc_sym) +
  geom_text_repel(data = lab_bot, aes(label = cell_type),
                  colour = "black", size = fontsize / ggplot2::.pt, seed = 1,
                  min.segment.length = 0, segment.size = 0.22, segment.colour = "grey60",
                  box.padding = 0.15, point.padding = 0.15, max.overlaps = Inf,
                  hjust = 1, direction = "y", nudge_x = n_b * 0.7 - lab_bot$rank,
                  nudge_y = -0.1 - lab_bot$icc_sym) +
  scale_colour_manual(values = COHORT_COL, labels = COHORT_LAB, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08)),
                     limits = c(-0.15, 0.9)) +
  labs(x = "Supertypes, ranked by SEA-AD Multiregion ICC",
       y = "Cross-region intraclass\ncorrelation (ICC)") +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  theme_cns +
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig <- (p_a / p_b) +
  plot_layout(heights = c(1.6, 1.3)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, "figure7_cross_region_reliability"), height_in = 6.5)
