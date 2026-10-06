# Figure 1. Supertype proportions are reproducible across snRNA-seq datasets for
# a limited set of supertypes, and those are more reproducible than their subclass.
#   (a) all 150 DFC 2026 supertypes ranked by median symmetric ICC across the five
#       independent matched-donor pairs (bar = IQR); labelled: the top N_LABEL and
#       the three supertypes highlighted in (b-d). Formerly Fig 2a (merged
#       2026-10-05; the nine Green x Mathys scatters of old Fig 2b go to the
#       supplement, script figure2/fig2.R).
#   (b) forest plots for three lineages (Immune = the DFC 2026 subclass that
#       holds microglia, Sst, L2/3 IT): symmetric ICC of the subclass and of each
#       of its supertypes across the five pairs: point = median, bar = IQR; the
#       dashed line marks the subclass median.
#   (c) donor-level CLR, Green vs Mathys (227 shared ROSMAP/RADC donors): subclass
#       next to its exemplar supertype = the highest median ICC among supertypes
#       that beat the subclass in all but at most one pair (consistency rule,
#       2026-10-05; Sst_25 loses only Mathys x Multiome, L2/3 IT_13 wins all five;
#       Multiome x PsychAD was used before and gives the Sst subclass an
#       unrepresentatively low value, 0.26 vs its median 0.35).
#   (d) same for BrainSCOPE CMC vs PsychAD (54 shared MSSM donors). PsychAD x
#       Ruzicka was tried first and dropped: Ruzicka has ~500 non-neuronal
#       nuclei per donor, so rare supertypes sit on the pseudocount floor.
#
# Metric (decision 2026-10-05): symmetric ICC = 2 cov / (var_a + var_b) of matched
# donors' CLR values (= Lin's concordance on centred data), recomputed from the
# donor caches by common/reliability_stats.py, which also checks r / icc_sym /
# icc_person_a against Shreejoy's 114 tables. Not icc_person_a (reference-dependent).
# Taxonomy: SEA-AD DFC 2026, reactive (-SEAAD) states included.
# Inputs (from fig1_data.py):
#   data/fig1a_pairs.tsv       per-pair icc_sym / rho / r / icc_person_a at both levels
#   data/fig1bc_clr.tsv        donor-level CLR for the pairs, both levels
#   common/data/reliability_summary.tsv  median ICC per type: source of (a), and
#                              (b) must agree with it
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig1.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(ggrepel)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure1"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)
WL <- file.path(dirname(HERE), "common", "data", "reliability_summary.tsv")
TAX <- "/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv"

SUBCLASSES <- c("Immune", "Sst", "L2/3 IT")
MIN_PAIRS_SHOW <- 4          # the highlighted supertype is measured in >= 4 of the 5 pairs
N_LABEL <- 5                 # (a): label the top N_LABEL plus the three (b-d) supertypes
PAIR_B <- "Green x Mathys"                # keys into fig1bc_clr.tsv$pair, not labels
PAIR_C <- "BrainSCOPE CMC x PsychAD"
LAB_B <- c(a = "Green RADC", b = "Mathys RADC")            # axis labels
LAB_C <- c(a = "BrainSCOPE CMC_MSSM", b = "MSSM_PsychAD")

# subclass = large filled square, supertype = small open circle (reactive
# -SEAAD states are included but not marked)
LEVEL_SHP <- c("Subclass" = 22, "Supertype" = 21)
LEVEL_SIZE <- c("Subclass" = 2.3, "Supertype" = 1.4)
PT_SIZE <- 1.4

# Cell class (taxonomy class_id), Allen-style hues (user decision 2026-10-05):
# excitatory = red/orange, inhibitory = blue/purple, non-neuronal = green.
# One hue per class so the 150-point ranking stays legible. NB: this differs
# from the Nature blue/red (neuronal/non-neuronal) used in Figs 3-6.
CLS_COL <- c("Excitatory" = "#D1491C", "Inhibitory" = "#2F4FA8", "Non-neuronal" = "#2E8B3C")

tax <- fread(TAX)
sc_class <- unique(tax[, .(subclass = subclass_label, class_id)])
stopifnot(!anyDuplicated(sc_class$subclass))
tax <- tax[, .(cell_type = supertype_label, class_id, reactive_state)]
class_of <- function(class_id) {
  factor(fifelse(class_id == "Glia", "Non-neuronal", class_id), levels = names(CLS_COL))
}

# ---- (a) per-pair ICC, subclass vs supertypes -------------------------------
a <- fread(file.path(HERE, "data/fig1a_pairs.tsv"))
stopifnot(uniqueN(a$pair) == 5)
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

sm <- a[, .(n_pairs = .N, median_icc = median(icc_sym),
            q25 = quantile(icc_sym, 0.25), q75 = quantile(icc_sym, 0.75)),
        by = .(subclass, level, cell_type, cls, lvl)]

# the medians must be the ones in the shared table
wl <- fread(WL)[, .(level, cell_type, wl_median = median_icc_sym, wl_n = n_pairs)]
chk <- merge(sm, wl, by = c("level", "cell_type"))
stopifnot(nrow(chk) == nrow(sm),
          max(abs(chk$median_icc - chk$wl_median)) < 1e-9,
          all(chk$n_pairs == chk$wl_n))

# order: subclass first, then supertypes by median ICC
setorder(sm, subclass, level, -median_icc)   # "subclass" < "supertype" alphabetically
sm[, cell_type := factor(cell_type, levels = unique(cell_type))]
a[, cell_type := factor(cell_type, levels = levels(sm$cell_type))]
sub_med <- sm[level == "subclass", .(subclass, sub_median = median_icc)]

# exemplar per lineage: highest median ICC among supertypes measured in
# >= MIN_PAIRS_SHOW pairs that beat their subclass in all but at most one of the
# pairs they share (the same one-pair allowance as MIN_PAIRS_SHOW)
sub_pair <- a[level == "subclass", .(subclass, pair, icc_sub = icc_sym)]
wins <- merge(a[level == "supertype"], sub_pair, by = c("subclass", "pair"))[
  , .(n_shared = .N, n_win = sum(icc_sym > icc_sub)), by = .(subclass, cell_type)]
cand <- merge(sm[level == "supertype" & n_pairs >= MIN_PAIRS_SHOW], wins,
              by = c("subclass", "cell_type"))
TOP <- cand[n_win >= n_shared - 1][order(-median_icc)][, .SD[1], by = subclass]
if (nrow(TOP) < length(SUBCLASSES)) {
  miss <- setdiff(SUBCLASSES, as.character(TOP$subclass))
  warning("no supertype beats its subclass in all but one pair for: ", paste(miss, collapse = ", "),
          "; falling back to highest median")
  TOP <- rbind(TOP, cand[subclass %in% miss][order(-median_icc)][, .SD[1], by = subclass])
}
cat("exemplar rule: wins in all but at most one shared pair, then highest median\n")
print(cand[order(subclass, -median_icc)][, .(subclass, cell_type, median_icc = round(median_icc, 3),
                                             n_win, n_shared)][, head(.SD, 4), by = subclass])
cat("subclass vs top supertype (median ICC over 5 pairs):\n")
print(merge(sub_med, TOP[, .(subclass, top = as.character(cell_type), top_median = median_icc)]))

# subclass: filled with its class colour; supertypes: white fill
sm[, fill_col := fifelse(level == "subclass", unname(CLS_COL[as.character(cls)]), "white")]
bold_subclass <- function(x)   # x-axis labels: subclass in bold, supertypes plain
  as.expression(lapply(x, function(s) if (s %in% SUBCLASSES) bquote(bold(.(s))) else s))

p_b <- ggplot(sm, aes(cell_type, median_icc)) +
  geom_hline(data = sub_med, aes(yintercept = sub_median),
             colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_linerange(aes(ymin = q25, ymax = q75, colour = cls), linewidth = 0.45) +
  geom_point(aes(colour = cls, fill = fill_col, shape = lvl, size = lvl), stroke = 0.5) +
  scale_colour_manual(values = CLS_COL, name = NULL) +
  scale_fill_identity() +
  scale_shape_manual(values = LEVEL_SHP, name = NULL) +
  scale_size_manual(values = LEVEL_SIZE, name = NULL) +
  scale_x_discrete(labels = bold_subclass) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  facet_wrap(~ subclass, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = "Intraclass correlation (ICC)\nacross five snRNA-seq pairs") +
  # only the subclass/supertype key lives here (below panel b); the class
  # colour key is drawn once under panel d
  guides(colour = "none", size = "none",
         shape = guide_legend(nrow = 1, override.aes = list(
           size = unname(LEVEL_SIZE), fill = c("grey40", "white"), colour = "grey40"))) +
  theme_cns +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom", legend.justification = "center",
        legend.margin = margin(0, 0, 0, 0))

# ---- (a) all supertypes ranked by median ICC (formerly Fig 2a) ---------------
rk <- fread(WL)[level == "supertype"]
rk <- merge(rk, tax, by = "cell_type", all.x = TRUE)
stopifnot(!anyNA(rk$class_id))
rk[, cls := class_of(class_id)]
setorder(rk, -median_icc_sym)
rk[, rank := .I]
LABEL <- union(rk[n_pairs >= MIN_PAIRS_SHOW][1:N_LABEL, cell_type], as.character(TOP$cell_type))
cat(sprintf("\n(a) %d supertypes (%d reactive); median ICC %.3f; %d with ICC >= 0.30; labelled: %s\n",
            nrow(rk), rk[reactive_state == TRUE, .N], median(rk$median_icc_sym),
            rk[median_icc_sym >= 0.30, .N], paste(LABEL, collapse = ", ")))
lab_a <- rk[cell_type %in% LABEL][order(rank)]

p_a <- ggplot(rk[order(as.integer(cls))], aes(rank, median_icc_sym)) +
  geom_linerange(aes(ymin = q25, ymax = q75), colour = GREYS[1], linewidth = 0.28) +
  geom_point(aes(colour = cls, fill = cls), shape = 21, size = 1.15, stroke = 0.28) +
  geom_text_repel(data = lab_a, aes(label = cell_type, colour = cls),
                  size = fontsize / ggplot2::.pt, show.legend = FALSE, seed = 1,
                  min.segment.length = 0, segment.size = 0.22, segment.colour = "grey60",
                  box.padding = 0.22, point.padding = 0.18, max.overlaps = Inf,
                  hjust = 0, direction = "y", nudge_x = 24 - lab_a$rank, ylim = c(0.32, NA)) +
  scale_colour_manual(values = CLS_COL, guide = "none") +
  scale_fill_manual(values = CLS_COL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.45))) +   # headroom for the labels
  labs(x = "Supertypes", y = "Intraclass correlation (ICC)\nacross five snRNA-seq pairs") +
  theme_cns + theme(legend.position = "none")

# ---- (c, d) donor-level CLR scatters -----------------------------------------
clr <- fread(file.path(HERE, "data/fig1bc_clr.tsv"))
show <- rbind(data.table(subclass = SUBCLASSES, cell_type = SUBCLASSES, level = "subclass"),
              TOP[, .(subclass = as.character(subclass), cell_type = as.character(cell_type),
                      level = "supertype")])
show <- merge(show, unique(sm[, .(cell_type = as.character(cell_type), cls)]), by = "cell_type")
show[, subclass := factor(subclass, levels = SUBCLASSES)]
setorder(show, subclass, level)          # subclass, then its top supertype, per lineage
clr <- merge(clr, show, by = c("cell_type", "level"))
clr[, cell_type := factor(cell_type, levels = show$cell_type)]

scatter_panel <- function(pair_label, lab, legend = FALSE) {
  d <- clr[pair == pair_label]
  stopifnot(uniqueN(d$cell_type) == nrow(show))
  st <- d[, { ac <- clr_a - mean(clr_a); bc <- clr_b - mean(clr_b)
              list(icc = 2 * cov(ac, bc) / (var(ac) + var(bc)), n = .N) }, by = cell_type]
  st[, lab := sprintf("ICC = %.2f", icc)]
  cat(sprintf("\n%s (n = %s donors):\n", pair_label, paste(unique(st$n), collapse = "/")))
  print(st[, .(cell_type, icc = round(icc, 3))])
  lims <- d[, {
    r <- range(c(clr_a, clr_b), na.rm = TRUE)
    pad <- max(diff(r) * 0.06, 0.15)
    list(clr_a = c(r[1] - pad, r[2] + pad), clr_b = c(r[1] - pad, r[2] + pad))
  }, by = cell_type]
  ggplot(d, aes(clr_a, clr_b)) +
    geom_blank(data = lims, aes(clr_a, clr_b), inherit.aes = FALSE) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55",
                linewidth = 0.3, linetype = "dashed") +
    geom_point(aes(colour = cls), size = 0.9, alpha = 0.7, stroke = 0, show.legend = legend) +
    geom_text(data = st, aes(label = lab), x = -Inf, y = Inf, hjust = -0.08, vjust = 1.25,
              size = fontsize / ggplot2::.pt, colour = "black", inherit.aes = FALSE) +
    scale_colour_manual(values = CLS_COL, name = NULL, drop = FALSE,
                        guide = if (legend) guide_legend(
                          nrow = 1, override.aes = list(size = PT_SIZE, alpha = 1)) else "none") +
    facet_wrap(~ cell_type, nrow = 1, scales = "free") +
    labs(x = sprintf("Centred log-ratio of cell-type proportion (%s)", lab["a"]),
         y = sprintf("Centred log-ratio of\ncell-type proportion \n (%s)", lab["b"])) +
    theme_cns +
    theme(legend.position = if (legend) "bottom" else "none",
          legend.justification = "center",
          panel.spacing.x = unit(7, "pt"), strip.clip = "off")
}
p_c <- scatter_panel(PAIR_B, LAB_B)
p_d <- scatter_panel(PAIR_C, LAB_C, legend = TRUE)

fig <- (p_a / p_b / p_c / p_d) +
  plot_layout(heights = c(1.25, 1.35, 1, 1), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, "figure1_snrna_reliability"), height_in = 10)
