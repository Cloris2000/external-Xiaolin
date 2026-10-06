# Figure 5 (DRAFT). Transferring CelMod across cohorts.
#   a  retention matrix: transferred model's median per-supertype r against the
#      target's own snRNA labels, divided by the same-cohort model's (same donors
#      and supertypes). In-sample diagonal excluded. CORE and all supertypes.
#   b  ROSMAP-trained model applied to MSBB: per-supertype r, transferred vs
#      same-cohort, CORE filled
#   c  reference-profile deconvolution (NNLS on SEA-AD profiles) vs CelMod, scored
#      in the ROSMAP triad
#   d  two bulk legs sharing a training cohort: inflation of their agreement
# Inputs from fig5_data.py. Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig5.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure5"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures"); dir.create(FIG, showWarnings = FALSE)

# ---- a ---------------------------------------------------------------------
a <- fread(file.path(DATA, "fig5a_retention.tsv"))
SRC_LEV <- c("ROSMAP", "MSBB", "NIMH HBCC", "MSBB + HBCC"); TGT_LEV <- c("ROSMAP", "MSBB IFG", "MSBB STG", "NIMH HBCC")
a[, source_lab := factor(source_lab, levels = rev(SRC_LEV))]
a[, target_lab := factor(target_lab, levels = TGT_LEV)]
a[, scope := factor(scope, levels = c("CORE", "all"), labels = c("CORE supertypes (n = 34)", "All supertypes (n = 149)"))]
grid <- CJ(source_lab = factor(rev(SRC_LEV), levels = rev(SRC_LEV)), target_lab = factor(TGT_LEV, levels = TGT_LEV),
           scope = levels(a$scope))
grid[, scope := factor(scope, levels = levels(a$scope))]
a <- merge(grid, a, by = c("source_lab", "target_lab", "scope"), all.x = TRUE)
p_a <- ggplot(a, aes(target_lab, source_lab)) +
  geom_tile(aes(fill = retention), colour = "white", linewidth = 0.6) +
  geom_tile(data = a[is.na(retention)], fill = "grey95", colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(retention), "", sprintf("%.2f", retention))),
            size = fontsize / .pt, colour = "black") +
  facet_wrap(~ scope, nrow = 1) +
  scale_fill_gradient2(low = "#F7F7F7", mid = "#B7C4D6", high = NEU, midpoint = 0.7,
                       limits = c(0.2, 1.2), name = "Retention",
                       guide = guide_colourbar(title.position = "top")) +
  scale_x_discrete(labels = function(x) sub(" ", "\n", x)) +
  scale_y_discrete(labels = function(x) sub(" \\+ ", " +\n", x)) +
  labs(x = "Applied to (bulk cohort)", y = "Trained on") +
  coord_fixed() +
  theme_cns + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    legend.position = "right", legend.key.height = unit(12, "pt"),
                    legend.key.width = unit(8, "pt"), legend.margin = margin(0, 0, 0, 0),
                    plot.margin = margin(2, 24, 2, 4))

# ---- b ---------------------------------------------------------------------
b <- fread(file.path(DATA, "fig5b_retention_by_type.tsv"))[source == "srcROSMAP" & target %in% c("MSBB_IFG", "MSBB_STG")]
b[, target_lab := factor(target_lab, levels = c("MSBB IFG", "MSBB STG"))]
st_b <- b[, .(med_x = median(r_same), med_t = median(r_transfer), n = .N), by = .(target_lab, core)]
print(st_b)
p_b <- ggplot(b, aes(r_same, r_transfer)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(aes(fill = core), shape = 21, size = 1.3, stroke = 0.3, colour = "grey30", alpha = 0.9) +
  facet_wrap(~ target_lab, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = NEU, `FALSE` = "white"),
                    labels = c(`TRUE` = "CORE supertype", `FALSE` = "other supertype"), name = NULL) +
  labs(x = "r, model trained on MSBB (same cohort)", y = "r, model trained on ROSMAP\n(transferred)") +
  coord_equal() +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- c ---------------------------------------------------------------------
c <- fread(file.path(DATA, "fig5c_profile_vs_celmod.tsv"))
c[, arm_lab := fcase(arm == "CelMod (Multiome-trained)", "CelMod\n(trained on matched pairs)",
                     arm == "profile_hier", "Reference profiles\n(hierarchical NNLS)",
                     arm == "profile_flat", "Reference profiles\n(flat NNLS)")]
cl <- melt(c, id.vars = "arm_lab", measure.vars = c("disatt_green_core", "h1_rel_core"),
           variable.name = "metric", value.name = "value")
cl[, metric := factor(metric, levels = c("disatt_green_core", "h1_rel_core"),
                      labels = c("Disattenuated r\nvs Green", "Triangulation\nreliability"))]
cl[, arm_lab := factor(arm_lab, levels = c("Reference profiles\n(flat NNLS)", "Reference profiles\n(hierarchical NNLS)", "CelMod\n(trained on matched pairs)"))]
p_c <- ggplot(cl, aes(value, arm_lab, fill = grepl("CelMod", arm_lab))) +
  geom_col(width = 0.6, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.15, size = fontsize / .pt) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = NEU, `FALSE` = "grey85"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45)), n.breaks = 3) +
  labs(x = "Score on the 34 CORE supertypes", y = NULL) +
  theme_cns

# ---- d ---------------------------------------------------------------------
d <- fread(file.path(DATA, "fig5d_trainer_inflation.tsv"))
cat(sprintf("d: median inflation all %+.3f, CORE %+.3f\n", median(d$inflation), median(d[core == TRUE]$inflation)))
p_d <- ggplot(d, aes(inflation, fill = core)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_histogram(binwidth = 0.05, colour = "white", linewidth = 0.2, boundary = 0) +
  geom_vline(xintercept = median(d$inflation), colour = "black", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, size = fontsize / .pt,
           label = sprintf("median %+.3f\nn = %d supertypes", median(d$inflation), nrow(d))) +
  scale_fill_manual(values = c(`TRUE` = NEU, `FALSE` = "grey80"),
                    labels = c(`TRUE` = "CORE", `FALSE` = "other"), name = NULL) +
  labs(x = "Agreement of two bulk legs: shared training cohort − different cohorts (Δr)",
       y = "Supertypes") +
  theme_cns + theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(0, 1))

design <- "
AA
BC
DD
"
fig <- wrap_plots(A = free(p_a, type = "panel"), B = free(p_b, type = "panel"),
                  C = free(p_c, type = "panel"), D = free(p_d, type = "panel"), design = design) +
  plot_layout(heights = c(1.1, 1.1, 0.7)) + plot_annotation(tag_levels = "a")
save_fig(fig, file.path(FIG, "figure5_transfer_draft"), height_in = 8.2)
