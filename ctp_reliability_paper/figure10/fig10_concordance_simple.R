# Figure 10 supplement (DRAFT, simple): do bulk-deconvolved and snRNA-seq composition give the
# same disease effects? One point per Core Supertype (86): case-vs-control beta from the Fig 10
# random-effects meta-analysis, snRNA-seq (x) vs bulk (y). Dashed line = identical effect.
# Annotation: Spearman rho and % of supertypes with the same sign (descriptive; supertypes are
# correlated, so no p value). Full analysis incl. same-donor permutation: fig10_concordance.R.
# Input: data/fig10_effects_meta.tsv. Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig10_concordance_simple.R
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(HERE, "..", "common", "theme_cns.R"))

TAX <- "/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv"
CLASS_LAB <- c("Neuronal: Glutamatergic" = "Excitatory", "Neuronal: GABAergic" = "Inhibitory",
               "Non-neuronal and Non-neural" = "Non-neuronal")
CLS_COL <- c(Excitatory = "#D1491C", Inhibitory = "#2F4FA8", "Non-neuronal" = "#2E8B3C")
DIS_LAB <- c(SCZ = "Schizophrenia", BD = "Bipolar disorder", AD = "Alzheimer's disease")

m <- fread(file.path(HERE, "data", "fig10_effects_meta.tsv"))[model == "primary"]
w <- dcast(m, disease + cell_type ~ modality, value.var = "beta")
stopifnot(!anyNA(w$bulk), !anyNA(w$snRNA), w[, .N, by = disease]$N == 86)
tax <- fread(TAX)
w[, class := factor(CLASS_LAB[tax$class_label[match(cell_type, tax$supertype_label)]], levels = CLASS_LAB)]
st <- w[, .(rho = cor(bulk, snRNA, method = "spearman"), same = mean(sign(bulk) == sign(snRNA))), by = disease]
st[, lab := sprintf("%s\nρ = %.2f", DIS_LAB[disease], rho)]
w[, panel := factor(st$lab[match(disease, st$disease)], levels = st$lab[match(names(DIS_LAB), st$disease)])]
lim <- range(c(w$bulk, w$snRNA))

p <- ggplot(w, aes(snRNA, bulk)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linewidth = 0.3, linetype = "dashed") +
  geom_point(aes(colour = class), size = 1.1, alpha = 0.85) +
  facet_wrap(~panel, nrow = 1) +
  coord_equal(xlim = lim, ylim = lim) +
  scale_colour_manual(values = CLS_COL, name = NULL) +
  labs(x = "snRNA-seq effect (β, SD units)", y = "Bulk-deconvolved effect\n(β, SD units)") +
  theme_cns + theme(legend.position = "bottom")

save_fig(p, file.path(HERE, "figures", "figure10_supp_concordance_simple_draft"), height_in = 2.9)
