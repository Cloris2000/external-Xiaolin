# Figure 4 (DRAFT). Deconvolution works where snRNA-seq is reliable, and that
# reliability lives at supertype, not subclass, resolution.
#   a  per supertype: out-of-fold CelMod accuracy (bulk vs matched snRNA, median
#      over 6 pairs) against snRNA-seq ICC (Fig 2 metric, median over 6 pairs)
#   b  triangulation reliability (ROSMAP triad, bulk leg): each subclass against
#      its own child supertypes
#   c  out-of-fold CelMod r by resolution, per bulk x snRNA pair
# Taxonomy: SEA-AD DFC 2026. Inputs from fig4_data.py.
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig4.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure4"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures"); dir.create(FIG, showWarnings = FALSE)

# ---- a ---------------------------------------------------------------------
a <- fread(file.path(DATA, "fig4a_icc_vs_deconv.tsv"))[!is.na(bulk_r)]
a[, key := class_key(compartment, reactive)]
ta <- spear(a$median_icc, a$bulk_r)
cat(sprintf("a: n=%d rho=%.2f p=%s\n", nrow(a), ta$estimate, fmt_p(ta$p.value)))
p_a <- ggplot(a[order(as.integer(key))], aes(median_icc, bulk_r)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", fill = "grey85",
              linewidth = 0.4, linetype = "dashed") +
  geom_point(aes(fill = key, colour = key, shape = key), size = 1.3, stroke = 0.3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = fontsize / .pt,
           label = sprintf("ρ = %.2f, p = %s\nn = %d supertypes", ta$estimate,
                           fmt_p(ta$p.value), nrow(a))) +
  scale_fill_manual(values = CLASS_COL, name = NULL) +
  scale_colour_manual(values = CLASS_OUT, name = NULL) +
  scale_shape_manual(values = CLASS_SHP, name = NULL) +
  labs(x = "snRNA-seq reliability of the supertype\n(median ICC across matched-donor pairs)",
       y = "Deconvolution accuracy of the supertype\n(median out-of-fold r, bulk vs snRNA-seq)") +
  theme_cns + theme(legend.position = "inside", legend.position.inside = c(1, 0),
                    legend.justification = c(1, 0))

# ---- b ---------------------------------------------------------------------
b <- fread(file.path(DATA, "fig4b_subclass_vs_children.tsv"))
sub <- b[, .(rel_sub = first(reliability_subclass), med_child = median(reliability),
             n_child = .N, resolved_sub = first(resolved_subclass),
             n_child_res = sum(resolved)), by = subclass][order(med_child)]
sub[, subclass := factor(subclass, levels = subclass)]
b[, subclass := factor(subclass, levels = levels(sub$subclass))]
cat(sprintf("b: subclass below its median child in %d/%d; resolved: %d/%d subclasses, %d/%d supertypes\n",
            sub[rel_sub < med_child, .N], nrow(sub), sum(sub$resolved_sub), nrow(sub),
            sum(b$resolved), nrow(b)))
XL <- c(-0.6, 1.6)
p_b <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(data = b, aes(pmin(pmax(reliability, XL[1]), XL[2]), subclass,
                           fill = resolved), shape = 21, size = 1.2, stroke = 0.25,
             colour = "grey40", position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  geom_point(data = sub, aes(pmin(pmax(rel_sub, XL[1]), XL[2]), subclass),
             shape = 23, size = 2.4, fill = GLIA, colour = "black", stroke = 0.35) +
  scale_fill_manual(values = c(`TRUE` = NEU, `FALSE` = "white"),
                    labels = c(`TRUE` = "resolved supertype", `FALSE` = "unresolved supertype"),
                    name = NULL) +
  coord_cartesian(xlim = XL) +
  labs(x = "Triangulation reliability\n(ROSMAP triad, bulk leg; axis clipped)", y = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 1.8))) +
  theme_cns + theme(legend.position = "bottom", legend.direction = "vertical",
                    legend.margin = margin(0, 0, 0, 0))

# ---- c ---------------------------------------------------------------------
c <- fread(file.path(DATA, "fig4c_oof_by_resolution.tsv"))
c[, set := factor(set, levels = c("Subclasses", "Other supertypes", "CORE supertypes"))]
c[, pair := gsub(" x ", " ×\n", pair)]
SET_COL <- c("Subclasses" = "#B7B7B7", "Other supertypes" = "#8491B4", "CORE supertypes" = NEU)
p_c <- ggplot(c, aes(set, r, colour = set)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, stroke = 0) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black", linewidth = 0.3, width = 0.55) +
  facet_wrap(~ pair, nrow = 1) +
  scale_colour_manual(values = SET_COL, name = NULL) +
  labs(x = NULL, y = "Out-of-fold r, deconvolved vs snRNA-seq") +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
                    axis.text.x = element_blank(), axis.ticks.x = element_blank())

design <- "
AB
CC
"
fig <- wrap_plots(A = free(p_a, type = "panel"), B = free(p_b, type = "panel"),
                  C = free(p_c, type = "panel"), design = design) +
  plot_layout(widths = c(1, 1.15), heights = c(1.25, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig, file.path(FIG, "figure4_resolution_draft"), height_in = 7.4)
