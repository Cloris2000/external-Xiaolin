# Figure 6 (DRAFT). The production classifier: CelMod trained on ROSMAP bulk with a
# consensus label over Green + Mathys + Diverse_multiome (tag `prod`).
# Message: averaging the snRNA-seq labels gives a better model than any one dataset,
# the production model holds up in cohorts it never saw, and it is applied to 22 cohorts.
# Conventions as Fig 5: Core Supertype = Fig 4 set (reliability-only; user decision
# 2026-10-28), drawn as filled vs open grey30 circles (fig4.R ICC_POS_SHP); the
# production / consensus model = #2F7E8F diamond (ROSMAP-trained colour in Fig 5); the
# single-dataset reference models = Fig 3's neutral slate #7B8C9E square; dataset and
# cohort names = common/cohort_names.tsv via cohort_label().
#   a  empty slot for the BioRender schematic of the production design (drawn by hand)
#   b  consensus label vs each single dataset's label, same ROSMAP donors and model,
#      scored against the snRNA-seq dataset held out of training (disattenuated r),
#      median over Core Supertypes; paired Wilcoxon per supertype
#   c  MSBB harness (never seen in training): production vs single-label references,
#      triangulation reliability, Core vs Non-core
#   d  cross-region agreement of deconvolved proportions vs SEA-AD's snRNA-only
#      prediction for the same supertype
#   e  the deliverable: deconvolved donors per bulk cohort
# Supplementary (figure6_supp_labels): snRNA-seq datasets behind each donor's label.
# Inputs from fig6_data.py. Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig6.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure6"
source(file.path(HERE, "..", "common", "theme_cns.R"))
DATA <- file.path(HERE, "data"); FIG <- file.path(HERE, "figures"); dir.create(FIG, showWarnings = FALSE)

PROD_COL <- "#2F7E8F"; REF_COL <- "#7B8C9E"                    # as Fig 5 (ROSMAP-trained / reference)
SET_LEV <- c("Core Supertype", "Non-core Supertype")
SET_SHP <- setNames(c(16, 1), SET_LEV)                         # as fig4.R ICC_POS_SHP
ds <- function(x) cohort_label(sub("_DFC$", "", x))            # snRNA dataset display name
fmt_p2 <- function(p) ifelse(p < 1e-3, "p < 0.001", sprintf("p = %.2g", p))

# ---- b: consensus vs single label ------------------------------------------
pb <- fread(file.path(DATA, "fig6b_per_type.tsv"))[core_supertype == TRUE & !is.na(r_disatt)]
tb <- fread(file.path(DATA, "fig6_core_supertype_tests.tsv"))[panel == "b"]
pb[, `:=`(held = ds(sub("exclude_", "", family)), kind = fifelse(grepl("__cons_", arm), "consensus", "single"))]
pb[, legs := sub("^(cons|solo)_", "", sub(".*__", "", arm))]
# the facet names the held-out dataset, so the consensus is of the other two
pb[, arm_lab := fifelse(kind == "consensus", "Consensus", ds(legs))]
sb <- pb[, .(r = median(r_disatt), n = .N), by = .(held, arm, arm_lab, kind)]
sb <- merge(sb, tb[, .(arm = other, p)], by = "arm", all.x = TRUE)
sb[, held := factor(paste("Held out:", held))]
sb[, arm_lab := factor(arm_lab, levels = c(sort(setdiff(unique(as.character(arm_lab)), "Consensus")), "Consensus"))]  # consensus on top
print(sb[order(held, kind)], digits = 3)
p_b <- ggplot(sb, aes(r, arm_lab)) +
  geom_segment(aes(x = 0, xend = r, yend = arm_lab), colour = "grey80", linewidth = 0.4) +
  geom_point(aes(shape = kind, fill = kind), colour = "black", size = 2.2, stroke = 0.4) +
  geom_text(data = sb[kind == "single"], aes(label = fmt_p2(p)), hjust = -0.25, size = fontsize / .pt) +
  facet_wrap(~ held, ncol = 1, scales = "free_y") +
  scale_shape_manual(values = c(consensus = 23, single = 22), labels = c("Consensus", "Single dataset"), name = NULL) +
  scale_fill_manual(values = c(consensus = PROD_COL, single = REF_COL), labels = c("Consensus", "Single dataset"), name = NULL) +
  scale_x_continuous(limits = c(0, NA), breaks = c(0, 0.2, 0.4), expand = expansion(mult = c(0, 0.55))) +
  labs(x = "Agreement with held-out snRNA-seq\n(disattenuated r, Core Supertypes)", y = "Training label") +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- c: MSBB harness --------------------------------------------------------
pc <- fread(file.path(DATA, "fig6c_per_type.tsv"))[!is.na(core_supertype) & !is.na(reliability)]
tc <- fread(file.path(DATA, "fig6_core_supertype_tests.tsv"))[panel == "c"]
pc[, set := factor(ifelse(core_supertype, SET_LEV[1], SET_LEV[2]), levels = SET_LEV)]
pc[, arm_lab := fcase(arm == "production consensus", "Production\n(consensus label)",
                      arm == "solo Green (reference)", "Green only",
                      arm == "solo Mathys (reference)", "Mathys only")]
sc <- pc[, .(rel = median(reliability), n = .N), by = .(arm, arm_lab, set)]
sc <- merge(sc, tc[, .(arm = other, p)], by = "arm", all.x = TRUE)
sc[set != SET_LEV[1], p := NA]                                 # tests are on Core Supertypes only
sc[, kind := fifelse(arm == "production consensus", "consensus", "single")]
sc[, arm_lab := factor(arm_lab, levels = c("Mathys only", "Green only", "Production\n(consensus label)"))]
print(sc[order(set, arm)], digits = 3)
p_c <- ggplot(sc, aes(rel, arm_lab)) +
  geom_segment(aes(x = 0, xend = rel, yend = arm_lab), colour = "grey80", linewidth = 0.4) +
  geom_point(aes(shape = kind, fill = kind), colour = "black", size = 2.2, stroke = 0.4) +
  geom_text(data = sc[!is.na(p)], aes(label = fmt_p2(p)), hjust = -0.25, size = fontsize / .pt) +
  facet_wrap(~ set, ncol = 1) +
  labs(y = "Training label") +
  scale_shape_manual(values = c(consensus = 23, single = 22), guide = "none") +
  scale_fill_manual(values = c(consensus = PROD_COL, single = REF_COL), guide = "none") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.45))) +
  labs(x = "Reliability in MSBB, unseen in training\n(triangulation, median over supertypes)") +
  theme_cns

# ---- d: cross-region agreement ---------------------------------------------
d <- fread(file.path(DATA, "fig6d_cross_region.tsv"))[!is.na(core_supertype)]
d[, set := factor(ifelse(core_supertype, SET_LEV[1], SET_LEV[2]), levels = SET_LEV)]
d[, fam_lab := fcase(family == "MSBB", "MSBB: FP, IFG, STG", family == "ROSMAP", "ROSMAP: DLPFC, PCC")]
st <- d[, {
  a <- spear(seaad_cross_region, bulk_cross_region); k <- spear(seaad_cross_region[core_supertype], bulk_cross_region[core_supertype])
  .(strip = sprintf("%s\nCore ρ = %.2f (%s)\nAll ρ = %.2f (%s)", fam_lab[1],
                    k$estimate, fmt_p2(k$p.value), a$estimate, fmt_p2(a$p.value)),
    n_core = sum(core_supertype), n_all = .N)
}, by = fam_lab]
print(st)
d <- merge(d, st, by = "fam_lab")
p_d <- ggplot(d, aes(seaad_cross_region, bulk_cross_region)) +
  geom_smooth(data = d[set == SET_LEV[1]], method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_point(aes(shape = set), colour = "grey30", size = 1.1, alpha = 0.8, stroke = 0.35) +
  facet_wrap(~ strip, nrow = 1) +
  scale_shape_manual(values = SET_SHP, breaks = SET_LEV, name = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Cross-region agreement predicted\nby SEA-AD snRNA-seq (r)",
       y = "Cross-region agreement of\ndeconvolved proportions (r)") +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
                    strip.clip = "off", panel.spacing.x = unit(14, "pt"))

# ---- e: the deliverable ----------------------------------------------------
e <- fread(file.path(DATA, "fig6e_deliverable.tsv"))
GRP_LAB <- c(rnaseq = "RNA-seq, with matched snRNA-seq cohorts", wp1_bulk = "RNA-seq, added 2026-09",
             gtex_regions = "GTEx (polyA)", microarray = "Microarray")
e[, `:=`(lab = cohort_label(cohort), grp = factor(GRP_LAB[group], levels = GRP_LAB))]
e[, lab := factor(lab, levels = e[order(grp, -n_donors)]$lab)]
p_e <- ggplot(e, aes(lab, n_donors, fill = grp)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n_donors), vjust = -0.4, size = fontsize / .pt * 0.85) +
  scale_fill_manual(values = setNames(c(GREYS[3], GREYS[2], GREYS[1], "white"), GRP_LAB), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, size = fontsize / .pt,
           label = sprintf("%d cohorts, %s donors", nrow(e), format(sum(e$n_donors), big.mark = ","))) +
  labs(x = NULL, y = "Donors with deconvolved\ncomposition") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_cns + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.86), legend.justification = c(1, 1),
                    axis.text.x = element_text(angle = 40, hjust = 1, size = fontsize * 0.85))

# ---- supplementary: label legs ---------------------------------------------
s <- fread(file.path(DATA, "fig6a_label_legs.tsv"))
s[, combo := sapply(strsplit(combo, " \\+ "), function(x) paste(ds(x), collapse = " + "))]
s[, combo := factor(combo, levels = s[order(n_legs, -n_donors)]$combo)]
s[, n_legs := factor(n_legs)]
p_s <- ggplot(s, aes(n_donors, combo, fill = n_legs)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = n_donors), hjust = -0.2, size = fontsize / .pt) +
  scale_fill_manual(values = c(`1` = GREYS[1], `2` = GREYS[2], `3` = GREYS[3]), name = "snRNA-seq datasets\nper donor") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = sprintf("ROSMAP donors with a consensus label (n = %d)", sum(s$n_donors)), y = NULL) +
  guides(fill = guide_legend(nrow = 1, title.position = "left")) +
  theme_cns + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

# ---- assemble --------------------------------------------------------------
# a: empty slot (BioRender); a blank ggplot, unlike plot_spacer(), still gets a panel letter
blank <- ggplot() + theme_void() + theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))
design <- "
AB
CD
EE
"
fig <- wrap_plots(A = blank, B = free(p_b, type = "panel"), C = free(p_c, type = "panel"),
                  D = free(p_d, type = "panel"), E = free(p_e, type = "panel"), design = design) +
  plot_layout(heights = c(1.25, 1, 0.75)) + plot_annotation(tag_levels = "a")
save_fig(fig, file.path(FIG, "figure6_production"), height_in = 9)
save_fig(p_s, file.path(FIG, "figure6_supp_labels"), height_in = 2.6)
