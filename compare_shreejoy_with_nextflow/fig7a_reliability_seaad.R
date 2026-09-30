# Regenerated Figure 7 panel A.
#   - Label the nine-type shortlist: median ICC >= 0.30 and bulk r >= 0.30.
#   - Include the seven -SEAAD types that have an ICC in only one pair
#     (whitelist v2 required >=2 pairs and dropped them).
#   - Draw neurons first so they do not cover glia.
#
# Writes under compare_shreejoy_with_nextflow/. Does not edit the wp3 tree.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
})

OUT_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

INDEP <- c(
  "Green x Mathys", "Green x Multiome2025", "Mathys x Multiome2025",
  "Multiome2025 x PsychAD", "BrainSCOPE CMC x PsychAD",
  "PsychAD x Ruzicka (MSSM2 x MSSM1)")

wl <- fread("/scratch/shreejoy/cell_type_bias/wp3/results/93_whitelist_v2.csv")
wl[, seaad := as.logical(seaad)]
wl[, reactive := seaad | grepl("SEAAD$|SEAD$", cell_type)]
wl[, single_pair := FALSE]
wl[, independent := TRUE]

pc <- fread("/scratch/shreejoy/cell_type_bias/wp3/results/51_pair_percelltype_supertype.csv")
pc[, reactive := grepl("SEAAD$|SEAD$", cell_type)]
seaad_all <- pc[reactive == TRUE]
have <- unique(wl$cell_type)
missing <- setdiff(unique(seaad_all$cell_type), have)

extra <- rbindlist(lapply(missing, function(ct) {
  rows <- seaad_all[cell_type == ct]
  indep <- rows[pair %in% INDEP]
  use <- if (nrow(indep)) indep else rows
  icc <- use$icc_person_a
  data.table(
    cell_type = ct,
    n_pairs = nrow(use),
    median_icc = median(icc),
    min_icc = min(icc),
    q25_icc = if (nrow(use) >= 2) as.numeric(quantile(icc, 0.25))
              else use$icc_ci_lo[1],
    q75_icc = if (nrow(use) >= 2) as.numeric(quantile(icc, 0.75))
              else use$icc_ci_hi[1],
    n_pairs_ge_030 = sum(icc >= 0.30),
    group = use$group[1],
    seaad = TRUE,
    bulk_r = NA_real_,
    bulk_n_pairs = NA_real_,
    in_v1 = FALSE,
    reactive = TRUE,
    single_pair = TRUE,
    independent = all(use$pair %in% INDEP)
  )
}))

wl <- rbind(wl, extra, fill = TRUE)
wl[, class_lab := fifelse(group == "glia", "Non-neuronal", "Neuronal")]
wl <- wl[order(-median_icc)]
wl[, rank := .I]
wl[, draw := fifelse(reactive, 3L, fifelse(group == "glia", 2L, 1L))]
wl <- wl[order(draw)]

# cnsplots Nature palette and ggplot theme (cns.setup_ggplot).
NEU <- "#3C5488"
NEU_MUTED <- "#8491B4"
GLIA <- "#E64B35"
fontsize <- 8

theme_cns <- theme_classic(base_size = fontsize, base_family = "sans") +
  theme(
    text = element_text(size = fontsize, colour = "black",
                        family = "sans", face = "plain"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.text = element_text(size = fontsize, colour = "black"),
    axis.title = element_text(size = fontsize, colour = "black", face = "plain"),
    legend.text = element_text(size = fontsize, colour = "black",
                               family = "sans", face = "plain"),
    legend.title = element_text(size = fontsize, colour = "black",
                                family = "sans", face = "plain", hjust = 0),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

wl[, key := fcase(
  reactive, "Reactive (-SEAAD)",
  class_lab == "Neuronal", "Neuronal",
  default = "Non-neuronal"
)]
wl[, key := factor(key, levels = c(
  "Neuronal", "Non-neuronal", "Reactive (-SEAAD)"))]
wl[, draw := as.integer(key)]
wl <- wl[order(draw)]

FILL <- c(
  "Neuronal" = NEU,
  "Non-neuronal" = GLIA,
  "Reactive (-SEAAD)" = GLIA)
COL <- c(
  "Neuronal" = NEU,
  "Non-neuronal" = GLIA,
  "Reactive (-SEAAD)" = "#111111")
SHP <- c(
  "Neuronal" = 21,
  "Non-neuronal" = 21,
  "Reactive (-SEAAD)" = 23)
PT_SIZE <- 1.15

lab <- copy(wl[!is.na(bulk_r) & median_icc >= 0.30 & bulk_r >= 0.30][order(rank)])
lab[, nudge_x := c(5, 14, 6, 14, 5, 16, 7, 15, 8)]
lab[, nudge_y := c(0.045, 0.055, 0.040, 0.050, 0.035, 0.045, 0.030, 0.040, 0.030)]

p <- ggplot(wl, aes(rank, median_icc)) +
  geom_linerange(data = wl[n_pairs >= 2],
                 aes(ymin = q25_icc, ymax = q75_icc),
                 colour = "#B7B7B7", linewidth = 0.28) +
  geom_point(aes(fill = key, colour = key, shape = key),
             size = PT_SIZE, stroke = 0.28) +
  geom_text_repel(
    data = lab,
    aes(label = cell_type, colour = key),
    size = fontsize / ggplot2::.pt, family = "sans", fontface = "plain",
    show.legend = FALSE,
    seed = 1, min.segment.length = 0, segment.size = 0.22,
    segment.colour = "grey60", box.padding = 0.22, point.padding = 0.18,
    max.overlaps = Inf, force = 0.9, force_pull = 0.7,
    nudge_x = lab$nudge_x, nudge_y = lab$nudge_y,
    direction = "both") +
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
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.key.size = unit(8, "pt"),
    legend.spacing.y = unit(1, "pt")
  )

png <- file.path(FIG_DIR, "fig7a_supertype_reliability_seaad.png")
pdf <- file.path(FIG_DIR, "fig7a_supertype_reliability_seaad.pdf")
ggsave(png, p, width = 7.4, height = 4.5, dpi = 300, bg = "white")
ggsave(pdf, p, width = 7.4, height = 4.5, bg = "white")
message("wrote ", png)
message("wrote ", pdf)

cat("\nlabelled shortlist (ICC >= 0.30 and bulk r >= 0.30):\n")
print(lab[order(rank),
          .(rank, cell_type, n_pairs,
            median_icc = round(median_icc, 3),
            bulk_r = round(bulk_r, 3), group)])
