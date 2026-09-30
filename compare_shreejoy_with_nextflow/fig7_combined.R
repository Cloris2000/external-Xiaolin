# Combined Figure 7: (a) supertype ICC ranking  (b) Green vs Mathys CLR
# Same 8 pt sans, Nature colours, and key mapping as the standalone panels.
# Writes under compare_shreejoy_with_nextflow/. Does not edit the wp3 tree.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
  library(patchwork)
})

OUT_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

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

# ---- (a) ranked median ICC -------------------------------------------------
INDEP <- c(
  "Green x Mathys", "Green x Multiome2025", "Mathys x Multiome2025",
  "Multiome2025 x PsychAD", "BrainSCOPE CMC x PsychAD",
  "PsychAD x Ruzicka (MSSM2 x MSSM1)")

wl <- fread("/scratch/shreejoy/cell_type_bias/wp3/results/93_whitelist_v2.csv")
wl[, seaad := as.logical(seaad)]
wl[, reactive := seaad | grepl("SEAAD$|SEAD$", cell_type)]
wl[, single_pair := FALSE]
pc <- fread("/scratch/shreejoy/cell_type_bias/wp3/results/51_pair_percelltype_supertype.csv")
pc[, reactive := grepl("SEAAD$|SEAD$", cell_type)]
missing <- setdiff(unique(pc[reactive == TRUE]$cell_type), unique(wl$cell_type))
extra <- rbindlist(lapply(missing, function(ct) {
  rows <- pc[cell_type == ct]
  indep <- rows[pair %in% INDEP]
  use <- if (nrow(indep)) indep else rows
  icc <- use$icc_person_a
  data.table(
    cell_type = ct, n_pairs = nrow(use),
    median_icc = median(icc),
    q25_icc = if (nrow(use) >= 2) as.numeric(quantile(icc, 0.25)) else use$icc_ci_lo[1],
    q75_icc = if (nrow(use) >= 2) as.numeric(quantile(icc, 0.75)) else use$icc_ci_hi[1],
    group = use$group[1], seaad = TRUE, bulk_r = NA_real_,
    reactive = TRUE, single_pair = TRUE)
}))
wl <- rbind(wl, extra, fill = TRUE)
wl[, class_lab := fifelse(group == "glia", "Non-neuronal", "Neuronal")]
wl <- wl[order(-median_icc)]
wl[, rank := .I]
wl[, key := fcase(
  reactive, "Reactive (-SEAAD)",
  class_lab == "Neuronal", "Neuronal",
  default = "Non-neuronal"
)]
wl[, key := factor(key, levels = names(FILL))]
wl <- wl[order(as.integer(key))]

lab <- copy(wl[!is.na(bulk_r) & median_icc >= 0.30 & bulk_r >= 0.30][order(rank)])
lab[, nudge_x := c(5, 14, 6, 14, 5, 16, 7, 15, 8)]
lab[, nudge_y := c(0.045, 0.055, 0.040, 0.050, 0.035, 0.045, 0.030, 0.040, 0.030)]

p_a <- ggplot(wl, aes(rank, median_icc)) +
  geom_linerange(data = wl[n_pairs >= 2],
                 aes(ymin = q25_icc, ymax = q75_icc),
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
    force = 0.9, force_pull = 0.7,
    nudge_x = lab$nudge_x, nudge_y = lab$nudge_y, direction = "both") +
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

# ---- (b) Green x Mathys CLR, same nine types -------------------------------
SHORTLIST <- c(
  "Astro_6-SEAAD", "OPC_2_2-SEAAD", "OPC_2",
  "Sst_25", "Sst_3", "Micro-PVM_2_1-SEAAD",
  "Sst_23", "Micro-PVM_2", "Sst_2")
dt <- fread(file.path(OUT_DIR, "data/fig7b_clr_scatter.tsv"))
dt <- dt[pair == "Green x Mathys" & cell_type %in% SHORTLIST]
dt[, cell_type := factor(cell_type, levels = SHORTLIST)]
dt[, reactive := grepl("SEAAD$", cell_type)]
dt[, key := fcase(
  reactive, "Reactive (-SEAAD)",
  group == "neu", "Neuronal",
  default = "Non-neuronal"
)]
dt[, key := factor(key, levels = names(FILL))]

st <- dt[, {
  ok <- is.finite(clr_a) & is.finite(clr_b)
  va <- var(clr_a[ok])
  list(icc = cov(clr_a[ok], clr_b[ok]) / va)
}, by = cell_type]
st[, lab := sprintf("ICC = %.2f", icc)]

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

png <- file.path(FIG_DIR, "fig7_reliability_and_clr.png")
pdf <- file.path(FIG_DIR, "fig7_reliability_and_clr.pdf")
ggsave(png, fig, width = 7.4, height = 9.8, dpi = 300, bg = "white")
ggsave(pdf, fig, width = 7.4, height = 9.8, bg = "white")
message("wrote ", png)
message("wrote ", pdf)
