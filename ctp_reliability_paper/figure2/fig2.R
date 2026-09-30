# Figure 2. Only a limited number of supertypes give a reliable snRNA-seq signal.
#   (a) median ICC per supertype across the six independent matched-donor pairs
#   (b) donor-level CLR, Green vs Mathys (227 shared donors), for the nine
#       supertypes with the highest median ICC (measured in >= 5 pairs)
#
# Layout and styling follow compare_shreejoy_with_nextflow/fig7_combined.R.
# Taxonomy: SEA-AD DFC 2026, reactive (-SEAAD) states included.
# Inputs:
#   117_whitelist_dfc.csv            median/IQR of ICC over the 6 independent pairs
#   data/fig2a_clr_green_mathys.tsv   from fig2a_clr_data.py (reproduces 114)
#   taxonomy_DFC_2026.tsv            compartment and reactive_state
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig2.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
  library(patchwork)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure2"
FIG_DIR <- file.path(HERE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)
WL <- "/scratch/shreejoy/cell_type_bias/wp3/results/117_whitelist_dfc.csv"
TAX <- "/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv"

N_SHOW <- 9
MIN_PAIRS_SHOW <- 5

# ---- style (fig7_combined.R) -------------------------------------------------
NEU <- "#3C5488"
GLIA <- "#E64B35"
fontsize <- 8
FILL <- c("Neuronal" = NEU, "Non-neuronal" = GLIA, "Reactive (-SEAAD)" = GLIA)
COL <- FILL
SHP <- c("Neuronal" = 21, "Non-neuronal" = 21, "Reactive (-SEAAD)" = 23)
PT_SIZE <- 1.15

theme_cns <- theme_classic(base_size = fontsize, base_family = "sans") +
  theme(
    text = element_text(size = fontsize, colour = "black",
                        family = "sans", face = "plain"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.text = element_text(size = fontsize, colour = "black", family = "sans"),
    axis.title = element_text(size = fontsize, colour = "black",
                              family = "sans", face = "plain"),
    strip.background = element_blank(),
    strip.text = element_text(size = fontsize, colour = "black",
                              family = "sans", face = "plain"),
    legend.text = element_text(size = fontsize, colour = "black",
                               family = "sans", face = "plain"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(8, "pt"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    plot.tag = element_text(size = fontsize, family = "sans",
                            face = "bold", colour = "black")
  )

tax <- fread(TAX)[, .(cell_type = supertype_label, compartment, reactive_state)]
key_of <- function(compartment, reactive) {
  factor(fcase(reactive, "Reactive (-SEAAD)",
               compartment == "neu", "Neuronal",
               default = "Non-neuronal"), levels = names(FILL))
}

# ---- (a) ranked median ICC --------------------------------------------------
wl <- fread(WL)
wl <- merge(wl, tax, by = "cell_type", all.x = TRUE)
stopifnot(!anyNA(wl$compartment), all(wl$group == wl$compartment))
wl[, key := key_of(compartment, reactive_state)]
setorder(wl, -median_icc)
wl[, rank := .I]

SHOW <- wl[n_pairs >= MIN_PAIRS_SHOW][1:N_SHOW, cell_type]
cat(sprintf("%d supertypes (%d reactive); median ICC %.3f; %d with ICC >= 0.30\n",
            nrow(wl), wl[reactive_state == TRUE, .N], median(wl$median_icc),
            wl[median_icc >= 0.30, .N]))
cat("shown in b:", paste(SHOW, collapse = ", "), "\n")

lab <- wl[cell_type %in% SHOW][order(rank)]
p_a <- ggplot(wl[order(as.integer(key))], aes(rank, median_icc)) +
  geom_linerange(aes(ymin = q25, ymax = q75),
                 colour = "#B7B7B7", linewidth = 0.28) +
  geom_point(aes(fill = key, colour = key, shape = key),
             size = PT_SIZE, stroke = 0.28) +
  geom_text_repel(
    data = lab,
    aes(label = cell_type, colour = key),
    size = fontsize / ggplot2::.pt, family = "sans", fontface = "plain",
    show.legend = FALSE, seed = 1, min.segment.length = 0,
    segment.size = 0.22, segment.colour = "grey60",
    box.padding = 0.22, point.padding = 0.18, max.overlaps = Inf,
    hjust = 0, direction = "y", nudge_x = 22 - lab$rank,
    ylim = c(0.30, NA)) +
  scale_fill_manual(values = FILL, name = NULL) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  labs(x = "Supertypes", y = "Intraclass correlation (ICC)") +
  guides(
    fill = guide_legend(nrow = 1, override.aes = list(size = PT_SIZE, stroke = 0.4)),
    colour = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)) +
  theme_cns +
  theme(legend.position = "bottom", legend.justification = "center")

# ---- (b) Green x Mathys CLR, same nine types --------------------------------
dt <- fread(file.path(HERE, "data/fig2a_clr_green_mathys.tsv"))
dt <- dt[cell_type %in% SHOW]
stopifnot(setequal(unique(dt$cell_type), SHOW))
dt <- merge(dt, tax, by = "cell_type")
dt[, key := key_of(compartment, reactive_state)]
dt[, cell_type := factor(cell_type, levels = SHOW)]

st <- dt[, {
  ac <- clr_a - mean(clr_a); bc <- clr_b - mean(clr_b)
  list(icc = cov(ac, bc) / var(ac), n = .N)
}, by = cell_type]
st[, lab := sprintf("ICC = %.2f", icc)]
cat("\nGreen x Mathys ICC for the panels in b:\n"); print(st)

lims <- dt[, {
  r <- range(c(clr_a, clr_b), na.rm = TRUE)
  pad <- max(diff(r) * 0.06, 0.15)
  list(clr_a = c(r[1] - pad, r[2] + pad),
       clr_b = c(r[1] - pad, r[2] + pad))
}, by = cell_type]

p_b <- ggplot(dt, aes(clr_a, clr_b)) +
  geom_blank(data = lims, aes(clr_a, clr_b), inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55",
              linewidth = 0.3, linetype = "dashed") +
  geom_point(aes(colour = key), size = 0.9, alpha = 0.7, stroke = 0,
             show.legend = FALSE) +
  geom_text(data = st, aes(label = lab),
            x = -Inf, y = Inf, hjust = -0.08, vjust = 1.25,
            size = fontsize / ggplot2::.pt, family = "sans",
            colour = "black", inherit.aes = FALSE) +
  scale_colour_manual(values = FILL, guide = "none") +
  facet_wrap(~ cell_type, nrow = 3, scales = "free") +
  labs(x = "Centred log-ratio of cell-type\nproportion (Green, ln)",
       y = "Centred log-ratio of cell-type\nproportion (Mathys, ln)") +
  theme_cns +
  theme(legend.position = "none")

p_a <- p_a + theme(plot.margin = margin(2, 4, 2, 4))
p_b <- p_b + theme(plot.margin = margin(2, 4, 2, 4))

fig <- (p_a / p_b) +
  plot_layout(heights = c(1, 2.1), widths = 1, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        legend.justification = "center",
        plot.tag = element_text(size = fontsize, family = "sans",
                                face = "bold", colour = "black"))

f <- file.path(FIG_DIR, "figure2_snrna_reliability")
ggsave(paste0(f, ".png"), fig, width = 7.4, height = 9.8, dpi = 300, bg = "white")
ggsave(paste0(f, ".pdf"), fig, width = 7.4, height = 9.8, bg = "white",
       device = cairo_pdf)
message("wrote ", f, ".{png,pdf}")
