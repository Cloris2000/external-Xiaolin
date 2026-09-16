# Nature-style Figure 2 focus: Panel A = all cell types; Panel B = large
# exemplar scatters with colored points (not 18 tiny 5-color clouds).
# Does not overwrite the original or the _1 figure.

OUT_DIR <- Sys.getenv(
  "FIG2_OUT_DIR",
  unset = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"
)
PAIRED_FILE <- file.path(OUT_DIR, "combined_bulk_snrna_paired.tsv")
ACC_FILE    <- file.path(OUT_DIR, "figure2_celltype_accuracy.tsv")

# One strong type per class, plus the next-best validated types.
EXEMPLARS <- c("sst", "vip", "lamp5", "microglia", "oligodendrocyte", "astrocyte")

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

make_clean_fn <- function(x) {
  sub("^_+|_+$", "",
      tolower(gsub("[^a-zA-Z0-9]+", "_",
                   gsub("([a-z])([A-Z])", "\\1_\\2", x))))
}

CT_DISPLAY <- c(
  oligodendrocyte = "Oligodendrocyte", opc = "OPC", astrocyte = "Astrocyte",
  microglia = "Microglia", endothelial = "Endothelial", pericyte = "Pericyte",
  vlmc = "VLMC",
  it = "IT", l4_it = "L4 IT", l5_6_np = "L5/6 NP", l5_et = "L5 ET",
  l6_ct = "L6 CT", l5_6_it_car3 = "L5/6 IT Car3", l6b = "L6b",
  pvalb = "PVALB", sst = "SST", vip = "VIP", lamp5 = "LAMP5", pax6 = "PAX6"
)
CT_CLASS <- c(
  oligodendrocyte = "Non-neuronal", opc = "Non-neuronal", astrocyte = "Non-neuronal",
  microglia = "Non-neuronal", endothelial = "Non-neuronal", pericyte = "Non-neuronal",
  vlmc = "Non-neuronal",
  it = "Excitatory", l4_it = "Excitatory", l5_6_np = "Excitatory",
  l5_et = "Excitatory", l6_ct = "Excitatory", l5_6_it_car3 = "Excitatory",
  l6b = "Excitatory",
  pvalb = "Inhibitory", sst = "Inhibitory", vip = "Inhibitory",
  lamp5 = "Inhibitory", pax6 = "Inhibitory"
)
CLASS_LEVELS <- c("Excitatory", "Inhibitory", "Non-neuronal")
COHORT_COLORS <- c(
  ROSMAP  = "#0072B2",
  HBCC    = "#E69F00",
  MSBB    = "#CC79A7",
  Mathys  = "#56B4E9",
  Ruzicka = "#009E73"
)
COHORT_LEVELS <- names(COHORT_COLORS)

base_th <- theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text       = element_text(size = 10),
    axis.title      = element_text(size = 12),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 12, face = "bold")
  )

paired_val <- read_tsv(PAIRED_FILE, show_col_types = FALSE) %>%
  mutate(
    ct_clean = make_clean_fn(cell_type),
    validation_cohort = factor(validation_cohort, levels = COHORT_LEVELS)
  ) %>%
  filter(snrna_proportion < 0.99, !ct_clean %in% "l4_it")

acc_mean_df <- read_tsv(ACC_FILE, show_col_types = FALSE) %>%
  filter(!ct_clean %in% "l4_it")

acc_df <- paired_val %>%
  group_by(validation_cohort, ct_clean) %>%
  summarise(
    pearson_r = suppressWarnings(
      cor(bulk_proportion, snrna_proportion,
          method = "pearson", use = "pairwise.complete.obs")),
    .groups = "drop"
  ) %>%
  filter(!is.na(pearson_r))

ord_tbl <- acc_mean_df %>%
  mutate(class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
  arrange(class_f, desc(mean_r))
ct_order_a <- ord_tbl$ct_clean

mean_a <- acc_mean_df %>%
  mutate(ct_ord  = factor(ct_clean, levels = ct_order_a),
         class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS))
coh_a <- acc_df %>%
  mutate(ct_ord  = factor(ct_clean, levels = ct_order_a),
         class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
  filter(!is.na(ct_ord))

panel_A <- ggplot(mean_a, aes(x = ct_ord, y = mean_r)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey55") +
  geom_errorbar(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r),
                width = 0.25, linewidth = 0.45, color = "grey60") +
  geom_point(data = coh_a,
             aes(x = ct_ord, y = pearson_r, color = validation_cohort),
             size = 1.8, alpha = 0.9, shape = 16, inherit.aes = FALSE,
             position = position_dodge(width = 0.35)) +
  geom_point(size = 3.0, shape = 18, color = "#08519c") +
  facet_grid(cols = vars(class_f), scales = "free_x", space = "free_x") +
  scale_color_manual(values = COHORT_COLORS, breaks = COHORT_LEVELS, guide = "none") +
  scale_x_discrete(labels = function(x) {
    lbl <- CT_DISPLAY[x]; ifelse(is.na(lbl), x, lbl)
  }) +
  labs(x = NULL, y = "Pearson r") +
  base_th +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        axis.text.x        = element_text(size = 10, angle = 45, hjust = 1),
        strip.text         = element_text(size = 11, face = "bold"),
        strip.background   = element_rect(fill = "grey94", color = NA),
        panel.spacing.x    = unit(0.4, "lines"),
        panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4))

ex_order <- acc_mean_df %>%
  filter(ct_clean %in% EXEMPLARS) %>%
  mutate(class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
  arrange(class_f, desc(mean_r)) %>%
  pull(ct_clean)

zdat <- paired_val %>%
  filter(ct_clean %in% EXEMPLARS) %>%
  mutate(
    sn_pct  = snrna_proportion * 100,
    bulk_au = bulk_proportion,
    ct_lab  = factor(CT_DISPLAY[ct_clean], levels = CT_DISPLAY[ex_order])
  ) %>%
  filter(is.finite(sn_pct), is.finite(bulk_au), !is.na(ct_lab))

r_facet <- acc_mean_df %>%
  filter(ct_clean %in% EXEMPLARS) %>%
  mutate(ct_lab = factor(CT_DISPLAY[ct_clean], levels = CT_DISPLAY[ex_order]),
         lbl    = paste0("r = ", sprintf("%.2f", mean_r))) %>%
  select(ct_lab, lbl)
pos_df <- zdat %>%
  group_by(ct_lab) %>%
  summarise(xpos = min(sn_pct, na.rm = TRUE) + 0.02 * diff(range(sn_pct, na.rm = TRUE)),
            ypos = max(bulk_au, na.rm = TRUE), .groups = "drop")
r_facet <- left_join(r_facet, pos_df, by = "ct_lab")

panel_B <- ggplot(zdat, aes(x = sn_pct, y = bulk_au, color = validation_cohort)) +
  geom_point(size = 1.35, alpha = 0.45, stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 0.9) +
  geom_text(data = r_facet,
            aes(x = xpos, y = ypos, label = lbl),
            inherit.aes = FALSE, hjust = 0, vjust = 1,
            size = 3.6, fontface = "bold", color = "grey15") +
  facet_wrap(~ ct_lab, nrow = 2, scales = "free") +
  scale_color_manual(values = COHORT_COLORS, name = "Cohort",
                     breaks = COHORT_LEVELS,
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = "snRNA-seq-derived proportion (%)",
       y = "Bulk-derived proportion (AU)") +
  base_th +
  theme(strip.text       = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "grey94", color = NA),
        panel.spacing    = unit(0.7, "lines"),
        axis.text        = element_text(size = 8.5),
        legend.position  = "bottom",
        panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.4))

fig2 <- plot_grid(
  panel_A, panel_B,
  nrow = 2, rel_heights = c(0.72, 1.15),
  labels = c("A", "B"), label_size = 15, label_fontface = "bold"
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

out_png <- file.path(OUT_DIR, "figure2_celltype_validation_focus_2.png")
out_pdf <- file.path(OUT_DIR, "figure2_celltype_validation_focus_2.pdf")
ggsave(out_png, fig2, width = 15, height = 11.5, dpi = 300, bg = "white")
ggsave(out_pdf, fig2, width = 15, height = 11.5, bg = "white")
cat("Saved", out_png, "\n")
cat("Saved", out_pdf, "\n")
