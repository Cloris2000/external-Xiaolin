# Figure 6 (DRAFT). The production classifier: CelMod trained on ROSMAP bulk with
# a consensus label over Green + Mathys + AMP-AD Diverse multiome (tag `prod`).
#   a  which snRNA-seq datasets contribute to each training donor's label
#   b  consensus label vs a single dataset's label, same donors, scored against
#      the held-out dataset (disattenuated r, CORE supertypes), three families
#   c  MSBB harness (never seen in training): production vs single-label references
#   d  cross-region demonstration: bulk-derived cross-region agreement of each
#      CORE supertype vs SEA-AD's snRNA-only prediction
#   e  the deliverable: deconvolved donors per bulk cohort
# Inputs from fig6_data.py. Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig6.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure6"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures"); dir.create(FIG, showWarnings = FALSE)

# ---- a ---------------------------------------------------------------------
a <- fread(file.path(DATA, "fig6a_label_legs.tsv"))
a[, combo := factor(combo, levels = a[order(n_legs, -n_donors)]$combo)]
a[, n_legs := factor(n_legs)]
p_a <- ggplot(a, aes(n_donors, combo, fill = n_legs)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = n_donors), hjust = -0.2, size = fontsize / .pt) +
  scale_fill_manual(values = c(`1` = "grey85", `2` = "#8491B4", `3` = NEU), name = "snRNA-seq\ndatasets per donor") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = sprintf("ROSMAP donors with a consensus label (n = %d)", sum(a$n_donors)), y = NULL) +
  guides(fill = guide_legend(nrow = 1, title.position = "left")) +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- b ---------------------------------------------------------------------
b <- fread(file.path(DATA, "fig6b_consensus_vs_single.tsv"))[scope == "CORE"]
b[, vs := fifelse(is.na(vs) | vs == "", NA_character_, vs)]  # empty strings in the TSV
bars <- b[is.na(vs)]
bars[, kind := ifelse(grepl("^consensus", arm_lab), "Consensus label", "Single-dataset label")]
bars[, arm_short := sub("^(consensus|single): ", "", arm_lab)]
tests <- b[!is.na(vs), .(family_lab, vs_short = sub(".*__solo_", "", vs), paired_diff, p)]
bars[, c("vs", "paired_diff", "p") := NULL]  # empty in the bar rows; the test rows carry them
bars <- merge(bars, tests, by.x = c("family_lab", "arm_short"), by.y = c("family_lab", "vs_short"), all.x = TRUE)
bars[, lab := ifelse(is.na(p), "", sprintf("p = %s", fmt_p(p)))]  # paired Wilcoxon vs the consensus arm
# singles at the top, consensus at the bottom of each facet (horizontal bars)
bars[, arm_short := factor(arm_short, levels = rev(unique(c(grep("\\+", arm_short, value = TRUE), sort(grep("\\+", arm_short, value = TRUE, invert = TRUE))))))]
print(bars[, .(family_lab, arm_short, median_r_disatt, paired_diff, p)])
p_b <- ggplot(bars, aes(median_r_disatt, arm_short, fill = kind)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = lab), hjust = -0.15, size = fontsize / .pt) +
  facet_wrap(~ family_lab, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Consensus label" = ACCENT, "Single-dataset label" = "grey85"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(y = "Training label", x = "Disattenuated r vs held-out dataset (CORE supertypes)") +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- c ---------------------------------------------------------------------
c <- fread(file.path(DATA, "fig6c_msbb_harness.tsv"))
c[, arm_lab := fcase(grepl("production", arm), "Production (consensus)", grepl("Green", arm), "Single: Green", grepl("Mathys", arm), "Single: Mathys")]
c[, arm_lab := factor(arm_lab, levels = c("Production (consensus)", "Single: Mathys", "Single: Green"))]
c[, scope := factor(scope, levels = c("CORE", "all"), labels = c("CORE (n = 34)", "All supertypes"))]
p_c <- ggplot(c, aes(median_reliability, arm_lab, fill = grepl("Production", arm_lab))) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", median_reliability)), hjust = -0.2, size = fontsize / .pt) +
  facet_wrap(~ scope, ncol = 1) +
  scale_fill_manual(values = c(`TRUE` = ACCENT, `FALSE` = "grey85"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4)), n.breaks = 4) +
  labs(y = "Reference applied to MSBB bulk", x = "Triangulation reliability (MSBB harness, median)") +
  theme_cns

# ---- d ---------------------------------------------------------------------
d <- fread(file.path(DATA, "fig6d_cross_region.tsv"))
d[, family := factor(family, levels = c("MSBB", "ROSMAP"), labels = c("MSBB: FP, IFG, STG", "ROSMAP: DLPFC, PCC"))]
st_d <- d[CORE == TRUE, {t <- spear(seaad_cross_region, bulk_cross_region)
                          .(lab = sprintf("CORE (n = %d)\nρ = %.2f, p = %.3f", .N, t$estimate, t$p.value))}, by = family]
print(st_d)
p_d <- ggplot(d, aes(seaad_cross_region, bulk_cross_region)) +
  geom_point(data = d[CORE == FALSE], colour = "grey80", size = 0.8, stroke = 0) +
  geom_smooth(data = d[CORE == TRUE], method = "lm", formula = y ~ x, colour = "grey30", fill = "grey85",
              linewidth = 0.4, linetype = "dashed") +
  geom_point(data = d[CORE == TRUE], fill = NEU, shape = 21, size = 1.5, stroke = 0.3, colour = "black") +
  geom_text(data = st_d, aes(label = lab), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
            size = fontsize / .pt, lineheight = 0.9, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  facet_wrap(~ family, nrow = 1) +
  labs(x = "Cross-region agreement predicted by SEA-AD snRNA-seq",
       y = "Cross-region agreement of\ndeconvolved proportions (r)") +
  theme_cns

# ---- e ---------------------------------------------------------------------
e <- fread(file.path(DATA, "fig6e_deliverable.tsv"))
GRP <- c(rnaseq = "Bulk RNA-seq, with matched snRNA-seq cohorts", wp1_bulk = "Bulk RNA-seq, added 2026-09",
         gtex_regions = "GTEx (polyA)", microarray = "Microarray")
e[, group := factor(GRP[group], levels = GRP)]
e[, cohort := factor(gsub("_", " ", cohort), levels = gsub("_", " ", e[order(group, n_donors)]$cohort))]
GRP_COL <- setNames(c(NEU, "#8491B4", "#B7B7B7", "#3D3D3D"), GRP)
p_e <- ggplot(e, aes(n_donors, cohort, fill = group)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n_donors), hjust = -0.2, size = (fontsize - 1) / .pt) +
  scale_fill_manual(values = GRP_COL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = sprintf("Donors with deconvolved composition (%d cohorts, %s donors)", nrow(e),
                   format(sum(e$n_donors), big.mark = ",")), y = NULL) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# 12-column grid: a and c are narrower than b and d
design <- "
AAAAABBBBBBB
CCCCCDDDDDDD
EEEEEEEEEEEE
"
fig <- wrap_plots(A = free(p_a, type = "panel"), B = free(p_b, type = "panel"),
                  C = free(p_c, type = "panel"), D = free(p_d, type = "panel"),
                  E = free(p_e, type = "panel"), design = design) +
  plot_layout(heights = c(1, 1, 1.5)) + plot_annotation(tag_levels = "a")
save_fig(fig, file.path(FIG, "figure6_production_draft"), height_in = 10)
