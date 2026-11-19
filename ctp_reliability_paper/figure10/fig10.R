# Figure 10 (DRAFT). Disease effects on the 86 Fig 4 Core Supertypes, measured two ways:
# bulk RNA-seq (production CelMod, out-of-sample donors only) vs snRNA-seq proportions
# (within-compartment CLR). Paired bars = random-effects meta beta +- SE of case vs control,
# outcome z-scored within cohort (SD units); model y ~ case + age + sex + PMI (no RIN, paper
# rule). SCZ/BP adults 18-69, AD all ages. * = BH FDR < 0.05 within disease x modality.
# Inputs: data/fig10_effects_meta.tsv, data/fig10_cohorts.tsv (fig10_data.py).
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig10.R
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(HERE, "..", "common", "theme_cns.R"))

TAX <- "/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv"
DIS_LAB <- c(SCZ = "SCZ", BD = "BP", AD = "AD")
CLASS_LAB <- c("Neuronal: Glutamatergic" = "Excitatory", "Neuronal: GABAergic" = "Inhibitory",
               "Non-neuronal and Non-neural" = "Non-neuronal")
MOD_LAB <- c(bulk = "Bulk RNA-seq (deconvolved)", snRNA = "snRNA-seq (proportion)")
MOD_COL <- c(bulk = GREYS[1], snRNA = GREYS[3])

m <- fread(file.path(HERE, "data", "fig10_effects_meta.tsv"))[model == "primary"]
co <- fread(file.path(HERE, "data", "fig10_cohorts.tsv"))[entered == TRUE]
tax <- fread(TAX)

# x order: class, then taxonomy (subclass, supertype) order
m[, class := factor(CLASS_LAB[class_label], levels = CLASS_LAB)]
stopifnot(!anyNA(m$class))
ord <- tax[supertype_label %in% m$cell_type]
ord[, class := factor(CLASS_LAB[class_label], levels = CLASS_LAB)]
setorder(ord, class, supertype_id)
m[, cell_type := factor(cell_type, levels = ord$supertype_label)]
stopifnot(!anyNA(m$cell_type), uniqueN(m$cell_type) == 86)

# row strips carry the donors entering each modality
n <- co[, .(n_case = sum(n_case), n_ctrl = sum(n_ctrl), k = .N), by = .(disease, modality)]
strip <- n[, .(lab = paste0(DIS_LAB[disease], "\n",
                            paste0(ifelse(modality == "bulk", "bulk ", "sn "), n_case, "/", n_ctrl,
                                   " (", k, ")", collapse = "\n"))), by = disease]
m[, disease_f := factor(strip$lab[match(disease, strip$disease)],
                        levels = strip$lab[match(names(DIS_LAB), strip$disease)])]
m[, modality := factor(modality, levels = names(MOD_LAB))]
m[, star_y := ifelse(beta >= 0, beta + se, beta - se) + sign(beta) * 0.04]

dodge <- position_dodge(width = 0.8)
p <- ggplot(m, aes(x = cell_type, y = beta, fill = modality)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(position = dodge, width = 0.8) +
  geom_errorbar(aes(ymin = beta - se, ymax = beta + se), position = dodge,
                width = 0, linewidth = 0.25, colour = "black") +
  geom_text(data = m[fdr < 0.05], aes(y = star_y, label = "*", group = modality),
            position = dodge, size = 2.4, vjust = 0.75, colour = "black") +
  facet_grid(disease_f ~ class, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = MOD_COL, labels = MOD_LAB, name = NULL) +
  labs(x = NULL, y = "Case vs control abundance change (β ± SE, SD units)") +
  theme_cns +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 6.5),
        panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(6, "pt"),
        legend.position = "top", legend.justification = "left")

save_fig(p, file.path(HERE, "figures", "figure10_disease_effects_draft"), height_in = 6.2)
