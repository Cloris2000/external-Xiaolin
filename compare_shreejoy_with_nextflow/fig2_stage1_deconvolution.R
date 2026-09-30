# Figure 2 -- Stage 1, deconvolution.
#
# Both pipelines measured the same quantity on donors who have both bulk and
# snRNA-seq: how well the deconvolved bulk proportion tracks the directly
# counted one. Mine is 19 broad MGP classes (figure2_celltype_accuracy.tsv);
# his is 33 CelMod supertypes on 1,601 dual-measured people
# (person_pheno/arm_agreement.csv).
#
# The comparison is not like for like, and the direction of the bias matters: a
# finer cell type is harder to deconvolve, so his numbers are achieved on the
# harder task. The figure says so rather than hiding it.
#
# Panel (c) is the point of the figure: agreement at stage 1 against the number
# of genome-wide hits the same trait produced at stage 3.

source(file.path(
  "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow",
  "theme_compare.R"))
suppressPackageStartupMessages(library(ggrepel))

dec  <- read_stage("stage1_deconvolution.tsv")
disc <- read_stage("stage3_discovery.tsv")
sub  <- read_stage("subclass_map.tsv")

dec[, r := as.numeric(r)]
med <- dec[, .(m = median(r)), by = pipeline]

# --- (a) every cell type, both pipelines, one axis --------------------------
dec[, lab := factor(label, levels = dec[order(r), label])]
dec[, pipe_lab := factor(
  c(mine = "Mine: MGP, broad class", his = "His: CelMod, supertype")[pipeline],
  levels = c("Mine: MGP, broad class", "His: CelMod, supertype"))]

pa <- ggplot(dec, aes(r, lab, colour = pipeline)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = r, yend = lab), linewidth = 0.5,
               alpha = 0.55) +
  geom_point(size = 1.9) +
  facet_grid(pipe_lab ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_colour_manual(values = PIPE_COL, guide = "none") +
  labs(x = "correlation, deconvolved bulk vs snRNA-seq", y = NULL,
       title = "a  Agreement between technologies, per cell type") +
  theme_cmp(10) +
  theme(axis.text.y = element_text(size = rel(0.72)),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90),
        panel.grid.major.y = element_line(colour = "grey94"))

# --- (b) the distributions, with medians ------------------------------------
pb <- ggplot(dec, aes(r, pipeline, fill = pipeline, colour = pipeline)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_boxplot(width = 0.42, alpha = 0.22, outlier.shape = NA,
               linewidth = 0.5) +
  geom_jitter(height = 0.13, size = 1.5, alpha = 0.75) +
  geom_text(data = med, aes(x = m, y = pipeline, label = sprintf("%.2f", m)),
            vjust = -1.9, fontface = "bold", size = 3.3, show.legend = FALSE) +
  scale_fill_manual(values = PIPE_COL, guide = "none") +
  scale_colour_manual(values = PIPE_COL, guide = "none") +
  scale_y_discrete(labels = c(mine = "19 broad\nMGP classes",
                              his = "33 SEA-AD\nsupertypes"),
                   limits = c("mine", "his")) +
  labs(x = "correlation, deconvolved bulk vs snRNA-seq", y = NULL,
       title = "b  He gets better agreement at finer resolution",
       subtitle = paste("Finer cell types are harder to deconvolve,",
                        "so this comparison understates the gap")) +
  theme_cmp(10)

# --- (c) stage 1 quality against stage 3 output -----------------------------
mine <- merge(
  dec[pipeline == "mine", .(mgp_trait, r)],
  disc[pipeline == "mine", .(mgp_trait, n_gw = as.integer(n_gw),
                             n_het = as.integer(n_gw_high_het))],
  by = "mgp_trait")

# Does the trait have any counterpart that survived his portability filter?
core_by_mgp <- sub[, .(n_core = sum(n_core)), by = mgp_trait]
mine <- merge(mine, core_by_mgp, by = "mgp_trait", all.x = TRUE)
mine[is.na(n_core), n_core := 0]
mine[, portable := ifelse(n_core > 0, "kept in his CORE set",
                          "dropped as non-portable")]

pc <- ggplot(mine, aes(r, n_gw + 1)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = portable, size = n_het + 1), alpha = 0.85) +
  geom_text_repel(aes(label = mgp_trait), size = 2.7, max.overlaps = 20,
                  segment.colour = "grey70", min.segment.length = 0.2) +
  scale_y_log10(breaks = c(1, 2, 6, 11, 51, 201, 501),
                labels = c(0, 1, 5, 10, 50, 200, 500)) +
  scale_size_area(max_size = 9, breaks = c(1, 51, 201, 351),
                  labels = c(0, 50, 200, 350),
                  name = "of which\nhigh-I2") +
  scale_colour_manual(values = c("dropped as non-portable" = "#C1553B",
                                 "kept in his CORE set" = "#2E6F8E"),
                      name = NULL) +
  labs(x = "stage 1: correlation, deconvolved bulk vs snRNA-seq",
       y = "stage 3: genome-wide hits in my meta-analysis",
       title = "c  My signal is loudest where the phenotype is weakest",
       subtitle = paste0(
         "L5.ET deconvolves at r = -0.01 and yields 473 hits, 351 high-I2. ",
         "SST deconvolves best\n(r = 0.44) and yields none. ",
         "520 of my 687 hits come from traits he dropped as non-portable.")) +
  theme_cmp(10) +
  theme(legend.position = "right", panel.grid.major.x = element_line(
    colour = "grey94"))

p <- (pa | (pb / pc + plot_layout(heights = c(1, 2.1)))) +
  plot_layout(widths = c(1, 1.45)) +
  plot_annotation(
    title = "Stage 1: MGP against CelMod, and what it costs downstream",
    caption = paste(
      "Mine: nextflow/manuscript_figure/figure2_celltype_accuracy.tsv,",
      "mean within-cohort r over 5 validation cohorts (1,299 donors).",
      "\nHis: ctpgwas/results/person_pheno/arm_agreement.csv,",
      "r over 1,601 dual-measured people. Portability from",
      "celltype-composition/refs/taxonomy_DFC_2026.tsv.",
      "\nHit counts from results/.../heterogeneity/meta_heterogeneity_summary.tsv."),
    theme = theme_cmp(12))

save_fig(p, "fig2_stage1_deconvolution.png", width = 15, height = 9.5)
