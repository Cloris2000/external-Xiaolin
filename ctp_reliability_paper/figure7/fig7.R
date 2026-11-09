# Figure 7. Dataset identity, not cortical region, drives cell-type-proportion
# discordance across snRNA-seq cohorts. For the three lineages used throughout
# (Immune = the DFC 2026 subclass holding microglia, Sst, L2/3 IT), compares
# two symmetric-ICC distributions for the SAME cell type:
#   - snRNA-seq pairs: the eight matched-donor pairs (2026-10-29; were five),
#     different datasets, same region (Figure 1's quantity; median + IQR
#     across the pairs). "Five" below refers to the earlier pair set. NB the
#     dataset side now comes from the registry composition counts, while the
#     region side (fig7_data.py) still reads the wp3 caches / atlas cells.
#   - region pairs (pooled): the Multiome RADC pair (AMP-AD Diverse multiome
#     2025, DLPFC vs STG, 398 donors) and the SEA-AD Multiregion pair (DFC vs
#     MTG, 80 donors), combined into one estimate per subclass by Fisher-z
#     inverse-variance pooling (2026-10-13 decision, pool_icc()) -- the
#     standard meta-analytic treatment of two independent correlation-like
#     estimates (cf. metafor's measure = "ZCOR"), weighting the larger
#     Multiome cohort by its greater precision rather than averaging the two
#     cohorts unweighted.
#   (a) subclass-level forest: one row-pair per lineage (snRNA-seq pairs vs
#       pooled region pairs) -- the headline comparison, clean and
#       uncluttered. A two-sided Fisher-z z-test comparing the pooled region
#       ICC against the pooled snRNA-seq-pair ICC is computed and printed to
#       the console (2026-10-27; see the test block below) but not drawn on
#       the panel.
#   (b) supertype-level scatter (2026-10-13 redesign, replaces the dense
#       per-supertype forest rows which did not fit ~37 cell types legibly;
#       2026-10-21: y switched from two cohort-coloured points per supertype
#       to one Fisher-z pooled point, since cohort-specific values are not of
#       interest here): x = ICC across different datasets (one region held
#       fixed; snRNA-seq pairs' median), y = ICC across different regions
#       (one dataset held fixed; Multiome + SEA-AD pool_icc()), one point per
#       supertype, faceted by lineage, diagonal reference line (y = x).
#       Points above the diagonal are supertypes that are MORE reliable
#       across regions than across datasets -- i.e. dataset identity, not
#       region, is the bigger source of discordance for that cell type.
# If the cross-region ICC (y) consistently sits above the cross-dataset ICC
# (x) -- panel a's subclass rows, and most points above the diagonal in panel
# b -- that is evidence that which dataset a sample comes from drives more of
# the cross-cohort discordance in cell-type proportions than which cortical
# region it was sampled from.
#
# Metric: symmetric ICC = 2 cov / (var_a + var_b) of matched donors' CLR
# values, identical formula and pseudocount across all three distributions
# (common/reliability_stats.py for the snRNA-seq pairs, fig7_data.py for the
# same-dataset region pairs). Taxonomy: SEA-AD DFC 2026.
#
# Former all-supertype forest (one triplet of rows per supertype, ~37 rows
# across the three facets) and the even-earlier all-~150-supertype ranked
# scatter are both dropped (2026-10-13): neither scaled to the full
# supertype set legibly. Panel (b)'s scatter carries the same supertype-level
# evidence at full resolution.
#
# Inputs:
#   data/fig7a_pairs.tsv                     same-dataset region-pair icc_sym
#                                             + donor-bootstrap q25/q75, per
#                                             cohort x level x cell type
#   common/data/reliability_summary.tsv      same-region snRNA-seq-pair
#                                             median_icc_sym + q25/q75 across
#                                             the five pairs, per level x cell
#                                             type (Figure 1's summary table)
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig7.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(ggrepel)
  library(systemfonts)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure7"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

# Register real Arial (found bundled with an unrelated ANSYS/ROOT module on
# Trillium, not a system font) so ggplot/svglite render actual Arial glyphs
# instead of fontconfig's "sans" fallback (Liberation Sans). This only points
# at the font file for rendering here -- do not embed/outline this font if
# the SVG is later shared, since it is Monotype-licensed font software
# bundled for that software's own use, not freely redistributable.
ARIAL_DIR <- "/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/root/6.14.04/fonts"
register_font(name = "Arial",
             plain = file.path(ARIAL_DIR, "arial.ttf"),
             bold = file.path(ARIAL_DIR, "arialbd.ttf"),
             italic = file.path(ARIAL_DIR, "ariali.ttf"),
             bolditalic = file.path(ARIAL_DIR, "arialbi.ttf"))
theme_cns <- theme_cns + theme(text = element_text(family = "Arial"))
WL <- file.path(dirname(HERE), "common", "data", "reliability_summary.tsv")

SUBCLASSES <- c("Immune", "Sst", "L2/3 IT")
N_LABEL <- 3   # per facet, label the most extreme |pooled region - dataset| gaps
LABEL_NUDGE <- 0.12   # ICC units; how far panel (c) labels are pushed off the diagonal

# pair counts (n=8 snRNA-seq pairs, n=2 region pairs) deliberately not in the
# legend: they are given in the panel (a) schematic (2026-10-27)
COMPARISONS <- c("Different datasets, same region",
                 "Different regions, same dataset")
COMP_COL <- setNames(c("#2E8B3C", "#44638C"), COMPARISONS)

# Fisher-z inverse-variance pooling of two independent icc_sym estimates
# (standard meta-analytic treatment of a correlation-like statistic, e.g.
# metafor's measure = "ZCOR"): z = atanh(icc), var(z) = 1 / (n - 3),
# weight = 1 / var. Weights the 398-donor Multiome cohort against the
# 80-donor SEA-AD cohort by precision rather than averaging them unweighted.
pool_icc <- function(icc, n) {
  z <- atanh(icc); w <- n - 3
  z_hat <- sum(w * z) / sum(w)
  se <- sqrt(1 / sum(w))
  data.table(icc = tanh(z_hat), q25 = tanh(z_hat - 1.96 * se), q75 = tanh(z_hat + 1.96 * se))
}

# ---- shared data: snRNA-seq pairs (Fig 1's table) + region pairs (fig7_data.py) --
wl <- fread(WL)[subclass %in% SUBCLASSES,
                .(level, subclass, cell_type, n_pairs,
                  icc = median_icc_sym, q25, q75)]
wl[, comparison := COMPARISONS[1]]

rp <- fread(file.path(HERE, "data/fig7a_pairs.tsv"))
rp <- rp[subclass %in% SUBCLASSES]
rp[, comparison := fifelse(cohort == "Multiome2025", "Multiome RADC\n(DLPFC vs STG)",
                           "SEA-AD Multiregion\n(DFC vs MTG)")]
rp <- rp[, .(level, subclass, cell_type, comparison,
            icc = icc_sym, q25, q75, n_donors)]

# ---- (a) subclass-level forest: snRNA-seq pairs vs pooled region pairs ------
region_pooled <- rp[level == "subclass", pool_icc(icc, n_donors), by = subclass]
region_pooled[, comparison := COMPARISONS[2]]
sub <- rbind(wl[level == "subclass", .(subclass, comparison, icc, q25, q75)],
            region_pooled[, .(subclass, comparison, icc, q25, q75)])
sub[, comparison := factor(comparison, levels = rev(COMPARISONS))]
sub[, subclass := factor(subclass, levels = rev(SUBCLASSES))]
stopifnot(nrow(sub) == length(SUBCLASSES) * length(COMPARISONS))

cat("\nMultiome + SEA-AD inputs to the pooled region estimate:\n")
print(rp[level == "subclass", .(subclass, comparison, icc = round(icc, 3), n_donors)][order(subclass)])

# ---- (a) statistical test: pooled region ICC vs pooled dataset ICC ----------
# (2026-10-27) Per subclass, two-sided z-test on the Fisher-z scale between
# the two inverse-variance pooled estimates -- the same pool_icc() weighting
# (w = n_donors - 3) applied to the five snRNA-seq pairs
# (common/data/pair_stats_percelltype.tsv) as to the two region pairs:
#   z = (z_region - z_dataset) / sqrt(se_region^2 + se_dataset^2)
# This is the standard test for a difference between two independent
# correlation-like estimates, extended to the pooled values. The panel still
# draws the snRNA-seq side as median + IQR across the five pairs (Figure 1's
# summary); the pooled dataset-side value is printed below for reference.
# Caveats (not corrected for): the five snRNA-seq pairs share datasets
# (Green, Multiome, PsychAD each appear in 2-3 pairs) so their pooled
# precision is somewhat overstated, and the Multiome DLPFC donors enter both
# sides of the comparison.
pp <- fread(file.path(dirname(HERE), "common", "data", "pair_stats_percelltype.tsv"))
pp <- pp[level == "subclass" & subclass %in% SUBCLASSES]
stopifnot(pp[, .N, by = subclass][, all(N == 8)])   # eight pairs since 2026-10-29
pool_z <- function(icc, n) {
  z <- atanh(icc); w <- n - 3
  list(z = sum(w * z) / sum(w), se = sqrt(1 / sum(w)))
}
tst <- merge(pp[, pool_z(icc_sym, n_donors), by = subclass],
             rp[level == "subclass", pool_z(icc, n_donors), by = subclass],
             by = "subclass", suffixes = c("_ds", "_rg"))
tst[, z_stat := (z_rg - z_ds) / sqrt(se_rg^2 + se_ds^2)]
tst[, p := 2 * pnorm(-abs(z_stat))]
# significant rows get asterisks (* P < 0.05, ** P < 0.01, *** P < 0.001);
# non-significant rows show the P value itself (2026-10-27 decision)
tst[, p_label := fcase(p < 0.001, "***",
                       p < 0.01,  "**",
                       p < 0.05,  "*",
                       default = sprintf("P = %.2f", p))]
tst[, subclass := factor(subclass, levels = rev(SUBCLASSES))]

cat("\npanel (a) test: pooled region vs pooled dataset ICC (Fisher-z z-test):\n")
print(tst[order(subclass), .(subclass, icc_dataset_pooled = round(tanh(z_ds), 3),
                             icc_region_pooled = round(tanh(z_rg), 3),
                             z = round(z_stat, 2), p = signif(p, 3))])

# The test result is printed above but NOT drawn on the panel (2026-10-27
# decision, after noting the blue bars are a 95% CI of the pooled estimate
# while the green bars are an IQR across pairs).
p_a <- ggplot(sub, aes(icc, subclass)) +
  geom_linerange(aes(xmin = q25, xmax = q75, colour = comparison),
                linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(aes(colour = comparison), size = 2.4, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3, linetype = "dashed") +
  scale_colour_manual(values = COMP_COL, name = NULL, breaks = COMPARISONS,
                      guide = guide_legend(ncol = 1)) +   # datasets first, then regions
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = "Symmetric intraclass correlation (ICC)", y = NULL) +
  theme_cns +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.margin = margin(0, 0, 0, 0))

cat("\nsubclass-level comparison (icc, [q25, q75]):\n")
print(sub[order(subclass, comparison), .(subclass, comparison, icc = round(icc, 3),
                                         q25 = round(q25, 3), q75 = round(q75, 3))])

# ---- (b) supertype-level scatter: same-region pairs vs pooled same-dataset pairs --
# y is the same Fisher-z pooling as panel (a)'s subclass rows (pool_icc()),
# applied per supertype instead -- not a cohort-coloured scatter any more
# (2026-10-21 decision): cohort-specific values are not of interest here,
# only each supertype's combined cross-region reliability.
region_pooled_st <- rp[level == "supertype", pool_icc(icc, n_donors), by = .(subclass, cell_type)]
st <- merge(wl[level == "supertype", .(subclass, cell_type, same_region_icc = icc)],
           region_pooled_st[, .(subclass, cell_type, same_dataset_icc = icc)],
           by = c("subclass", "cell_type"))
# Monocyte (SEA-AD-only) has no snRNA-seq-pairs row and drops out of this
# inner merge; L2/3 IT_14-SEAAD (Multiome-only) keeps its pool_icc() value,
# which for a single input is just that cohort's own icc_sym. Both documented
# in fig7_data.py.
st[, subclass := factor(subclass, levels = SUBCLASSES)]
st[, gap := same_dataset_icc - same_region_icc]

# label the top N_LABEL most extreme dataset-vs-region gaps per facet
lab <- st[order(-abs(gap)), head(.SD, N_LABEL), by = subclass]

lims <- range(c(st$same_region_icc, st$same_dataset_icc))
pad <- diff(lims) * 0.06
# The y = x reference line is a segment with explicit endpoints at the panel
# corners rather than geom_abline() (2026-10-27): geom_abline draws far past
# the panel and relies on the SVG clip-path to trim it, and the clip-paths
# are stripped below for Affinity, which left dashed diagonals running across
# the whole SVG. Endpoints = coord limits plus ggplot's default 5% expansion
# so the segment still spans the full panel.
diag_ends <- (lims + c(-pad, pad)) + c(-0.05, 0.05) * diff(lims + c(-pad, pad))

p_b <- ggplot(st, aes(same_region_icc, same_dataset_icc)) +
  annotate("segment", x = diag_ends[1], y = diag_ends[1],
           xend = diag_ends[2], yend = diag_ends[2],
           colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_point(colour = "black", size = 1.3, alpha = 0.8, stroke = 0) +
  # Labels are nudged perpendicular to the y = x diagonal, away from it (up-left
  # for points above, down-right for points below), with a larger repel
  # padding and darker leaders (2026-10-27): the labelled points are by
  # construction the furthest from the diagonal, and pushing the labels in
  # that same direction keeps each leader unambiguous within clusters such as
  # the leftmost Sst points.
  geom_text_repel(data = lab, aes(label = cell_type),
                  nudge_x = -LABEL_NUDGE * sign(lab$gap),
                  nudge_y =  LABEL_NUDGE * sign(lab$gap),
                  colour = "black", size = fontsize / ggplot2::.pt, seed = 1,
                  min.segment.length = 0, segment.size = 0.3, segment.colour = "grey40",
                  box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf,
                  # ggrepel defaults to max.time = 0.5 s of wall time, so on a
                  # loaded login node the layout can stop before it converges
                  # and differ between the PNG and SVG prints (2026-10-27: the
                  # Immune labels overlapped in the SVG only). A large budget
                  # makes the layout converge every time.
                  max.time = 10, max.iter = 1e5) +
  coord_equal(xlim = lims + c(-pad, pad), ylim = lims + c(-pad, pad)) +
  facet_wrap(~ subclass, nrow = 1) +
  labs(x = "Different datasets, same region\n(median ICC across 8 snRNA-seq pairs)",
       y = "Different regions, same dataset\n(Multiome + SEA-AD, pooled)") +
  theme_cns +
  theme(legend.position = "none")

n_above <- st[, .(n = .N, n_above = sum(gap > 0))]
cat("\nsupertypes with pooled same-dataset ICC > same-region ICC:\n")
print(n_above)

# row 1: blank space at left (45%, tagged "a") reserved for a schematic to be
# added manually afterward, panel (a)'s forest plot narrower at right (55%,
# tagged "b"); row 2: panel (b)'s scatter, full width, tagged "c" (2026-10-21:
# tags shifted down one letter so the reserved schematic space is itself "a").
# plot_spacer() has no panel/axes to anchor plot.tag against and the tag does
# not render on it; an empty ggplot with theme_void() (still a full plot
# object, just with blank geometry) does carry the tag correctly.
blank_a <- ggplot() + theme_void() + labs(tag = "a")
row1 <- blank_a + (p_a + labs(tag = "b")) + plot_layout(widths = c(0.45, 0.55))
fig <- (row1 / (p_b + labs(tag = "c"))) +
  plot_layout(heights = c(0.55, 1.3)) &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, "figure7_cross_region_reliability"), height_in = 5.2)

# SVG (vector, for manually adding the schematic in row 1's blank space):
# svglite (mamba-installed into r_env, 2026-10-21 -- install.packages() fails
# here with a missing-libz link error; use mamba for new R packages on
# Trillium) gives proper <text> strings instead of base R's svg() device,
# which positions each character individually and is harder to edit by hand.
# fixed = FALSE keeps absolute text positioning consistent with how theme_cns
# sizes text for the PDF/PNG (fixed = TRUE rescales for cross-platform font
# metrics, which is not needed here since this SVG is edited on this machine).
library(svglite)
svg_path <- file.path(FIG_DIR, "figure7_cross_region_reliability.svg")
svglite(svg_path, width = 183 / 25.4, height = 5.2, bg = "white", fix_text_size = FALSE)
print(fig)
dev.off()

# Affinity Designer's SVG importer does not apply svglite's shared <style>
# block (the ".svglite line, .svglite polyline, ..." CSS-class default of
# stroke: #000000; fill: none), so axis lines/ticks/panel borders -- which
# have no per-element stroke colour and rely entirely on that stylesheet --
# import with no visible stroke, while coloured geoms (points, which set
# stroke/fill inline) are unaffected. Inline the same defaults onto every
# bare line/polyline/path/rect/circle element so the SVG renders correctly
# without depending on CSS class cascading.
lines <- readLines(svg_path)
bare_tag <- "^\\s*<(line|polyline|path|rect|circle)\\b"
has_style <- grepl("style='", lines)
is_bare_el <- grepl(bare_tag, lines) & has_style
lines[is_bare_el] <- sub("style='", "style='stroke: #000000; fill: none; ", lines[is_bare_el])

# Every ggplot panel (faceted or not) is wrapped by svglite in a
# <g clip-path='url(#...)'> referencing a <clipPath><rect .../></clipPath>
# that defines the panel's data-area boundary. Affinity Designer's SVG
# importer renders that clip-path reference as a visible stroked box around
# each such group -- the "panel border" seen around every panel on import --
# regardless of the referenced rect's own style (confirmed: it affects every
# panel, not just facets, so it is the clip-path mechanism itself, not a
# facet-specific issue). Since this is a static published figure (not
# re-rendered at arbitrary sizes), dropping the clip-path attribute is safe:
# the only effect is that a geom drawn exactly at a panel's edge could in
# principle overflow it very slightly, which does not occur in this figure.
lines <- sub("<g clip-path='url\\(#[^)]*\\)'>", "<g>", lines)
writeLines(lines, svg_path)
message("wrote ", svg_path, " (inlined default stroke/fill, dropped clip-path refs for Affinity compatibility)")
