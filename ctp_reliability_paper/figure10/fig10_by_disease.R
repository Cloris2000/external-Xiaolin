# Figure 10 (DRAFT), one figure per disorder: disease effects on the 86 Fig 4 Core Supertypes,
# bulk RNA-seq (top) and snRNA-seq (bottom) as separate bar plots in one column, sharing the
# x axis (supertypes) and the y scale. Bars coloured by subclass, as shades of the class colour
# (Excitatory orange, Inhibitory blue, Non-neuronal green). Same numbers as fig10.R:
# random-effects meta beta +- SE, case vs control, outcome z-scored within cohort (SD units),
# y ~ case + age + sex + PMI; SCZ/BP adults 18-69, AD all ages. * = BH FDR < 0.05 within
# disease x modality. Donor counts are printed for the legend text, not drawn.
# Inputs: data/fig10_effects_meta.tsv, data/fig10_cohorts.tsv (fig10_data.py).
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig10_by_disease.R
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(HERE, "..", "common", "theme_cns.R"))

TAX <- "/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv"
DIS <- c(SCZ = "Schizophrenia", BD = "Bipolar disorder", AD = "Alzheimer's disease")
FILE <- c(SCZ = "SCZ", BD = "BP", AD = "AD")
CLASS_LAB <- c("Neuronal: Glutamatergic" = "Excitatory", "Neuronal: GABAergic" = "Inhibitory",
               "Non-neuronal and Non-neural" = "Non-neuronal")
MOD_LAB <- c(bulk = "Bulk RNA-seq (deconvolved)", snRNA = "snRNA-seq")
# light -> dark ramp per class, centred on the paper's class colours (#D1491C, #2F4FA8, #2E8B3C)
RAMP <- list(Excitatory = c("#F4BE9F", "#D1491C", "#6E230A"),
             Inhibitory = c("#AFC0EC", "#2F4FA8", "#16244F"),
             "Non-neuronal" = c("#ABDCAE", "#2E8B3C", "#124A19"))

m <- fread(file.path(HERE, "data", "fig10_effects_meta.tsv"))[model == "primary"]
co <- fread(file.path(HERE, "data", "fig10_cohorts.tsv"))[entered == TRUE]
tax <- fread(TAX)

# x order: class, then taxonomy order (supertype_id keeps subclasses contiguous)
ord <- tax[supertype_label %in% m$cell_type]
ord[, class := factor(CLASS_LAB[class_label], levels = CLASS_LAB)]
setorder(ord, class, supertype_id)
stopifnot(nrow(ord) == 86, !anyNA(ord$class))
stopifnot(ord[, .(k = uniqueN(rleid(subclass_label))), by = subclass_label]$k == 1)  # contiguous

# subclass colours: shades of the class ramp; interleaved so neighbouring subclasses differ
sub <- unique(ord[, .(class, subclass_label)])
sub[, col := {
  k <- .N; sh <- colorRampPalette(RAMP[[as.character(class[1])]])(k)
  sh[order(c(seq(1, k, 2), seq(2, k, 2)))]
}, by = class]
SUB_COL <- setNames(sub$col, sub$subclass_label)

m[, class := factor(CLASS_LAB[class_label], levels = CLASS_LAB)]
m[, subclass := factor(subclass_label, levels = sub$subclass_label)]
m[, cell_type := factor(cell_type, levels = ord$supertype_label)]
stopifnot(!anyNA(m$class), !anyNA(m$subclass), !anyNA(m$cell_type))
m[, star_y := ifelse(beta >= 0, beta + se, beta - se) + sign(beta) * 0.03]

bars <- function(d, mod, ylim, show_x) {
  x <- m[disease == d & modality == mod]
  stopifnot(nrow(x) == 86)
  g <- ggplot(x, aes(x = cell_type, y = beta)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_col(aes(fill = subclass), width = 0.8) +
    geom_errorbar(aes(ymin = beta - se, ymax = beta + se), width = 0, linewidth = 0.25) +
    geom_text(data = x[fdr < 0.05], aes(y = star_y, label = "*"), size = 2.4, vjust = 0.75) +
    facet_grid(. ~ class, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = SUB_COL, name = NULL, drop = FALSE) +
    coord_cartesian(ylim = ylim) +
    labs(x = NULL, y = "Abundance change\n(β ± SE, SD units)", subtitle = MOD_LAB[[mod]]) +
    theme_cns +
    theme(plot.subtitle = element_text(size = fontsize, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          strip.text = element_text(size = 6.5), panel.spacing.x = unit(3, "pt"))
  if (!show_x) g <- g + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  else g <- g + theme(strip.text = element_blank(), strip.background = element_blank())
  g
}

for (d in names(DIS)) {
  x <- m[disease == d]
  ylim <- range(c(x$beta - x$se, x$beta + x$se, x$star_y, 0))
  fig <- (bars(d, "bulk", ylim, FALSE) / bars(d, "snRNA", ylim, TRUE)) +
    plot_layout(guides = "collect") +
    plot_annotation(title = DIS[[d]],
                    theme = theme(plot.title = element_text(size = fontsize + 1, face = "bold"))) &
    theme(legend.position = "bottom", legend.key.size = unit(6, "pt"),
          legend.text = element_text(size = 5.5, margin = margin(l = 1.5, r = 4)),
          legend.margin = margin(0, 0, 0, 0)) &
    # one legend row per class (9 excitatory, 9 inhibitory, 6 non-neuronal subclasses)
    guides(fill = guide_legend(ncol = max(table(sub$class)), byrow = TRUE))
  save_fig(fig, file.path(HERE, "figures", paste0("figure10_disease_effects_", FILE[[d]], "_draft")),
           height_in = 4.8)
}

# numbers for the legend text (not drawn)
cat("\n-- donors per disease x modality (entered cohorts) --\n")
print(co[, .(n_cohorts = .N, n_case = sum(n_case), n_ctrl = sum(n_ctrl)), by = .(disease, modality)])
cat("\n-- FDR < 0.05 supertypes (of 86) --\n")
print(m[, .(n_fdr05 = sum(fdr < 0.05)), by = .(disease, modality)])
