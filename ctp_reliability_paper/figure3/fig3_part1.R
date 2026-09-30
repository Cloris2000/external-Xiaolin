# Figure 3, part 1 (DRAFT). Tissue blocks sample different cortical depths,
# and that lowers agreement between datasets.
#   a  depth distribution of sampled nuclei, per dataset
#   b  the same donor sampled by two studies: sampled depth in each
#   c  depth gap between a donor's two blocks vs agreement of that donor's
#      inhibitory-neuron composition (depth from excitatory neurons only, so the
#      two measures use different cells)
#
# Taxonomy: SEA-AD DFC 2026. Inputs from fig3_data.py.
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig3_part1.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure3"
DATA <- file.path(HERE, "data")
FIG_DIR <- file.path(HERE, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

# Olah 2020 is microglia-enriched (92% Micro-PVM_2), so its composition says
# nothing about where in the cortex the block was taken.
EXCLUDE <- "Olah_2020"
REF <- "SEAAD_DFC_2026"
# Display names. Only SEA-AD is renamed for now; the others wait on confirmed
# citations (first author + year, or consortium name).
LABEL <- c(SEAAD_DFC_2026 = "SEA-AD DLPFC", Multiome2025_DLPFC = "Multiome 2025",
           BrainSCOPE_CMC = "BrainSCOPE CMC", BrainSCOPE_UCLA = "BrainSCOPE UCLA",
           BrainSCOPE_Multiome = "BrainSCOPE Multiome")
pretty <- function(x) ifelse(x %in% names(LABEL), LABEL[x], gsub("_", " ", x))

NEU <- "#3C5488"
fontsize <- 8
# Study pairs use the cnsplots "Cell" palette; the "Nature" blue/red stay
# reserved for neuronal/non-neuronal, as in Figure 2.
PAIR_COL <- c("Green x Mathys" = "#2F7E8F", "Green x Multiome" = "#E1A22E",
              "Mathys x Multiome" = "#8B6FA8", "Multiome x PsychAD" = "#C84C3A",
              "BrainSCOPE CMC x PsychAD" = "#5F9862", "PsychAD x Ruzicka" = "#4E5A8A")
pair_lab <- function(x) gsub(" x ", " vs ", x)
theme_cns <- theme_classic(base_size = fontsize, base_family = "sans") +
  theme(
    text = element_text(size = fontsize, colour = "black"),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.text = element_text(size = fontsize, colour = "black"),
    axis.title = element_text(size = fontsize, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = fontsize, colour = "black"),
    legend.title = element_text(size = fontsize, colour = "black"),
    legend.text = element_text(size = fontsize, colour = "black"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(8, "pt"),
    plot.tag = element_text(size = fontsize, face = "bold", colour = "black"),
    plot.margin = margin(2, 4, 2, 4)
  )
spear <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
fmt_p <- function(p) ifelse(p < 1e-3, sprintf("%.0e", p), sprintf("%.2g", p))

# ---- a: depth distribution per dataset --------------------------------------
dd <- fread(file.path(DATA, "fig3_donor_depth.tsv"))[!dataset %in% EXCLUDE]
prof <- fread(file.path(DATA, "fig3_dataset_profiles.tsv"))[!dataset %in% EXCLUDE]
ord <- dd[, .(m = median(depth_all), n = .N), by = dataset]
# deepest at the bottom, shallowest above, reference on top
ord <- rbind(ord[dataset != REF][order(-m)], ord[dataset == REF])
ord[, lab := sprintf("%s (n = %d)", pretty(dataset), n)]
prof[, y := match(dataset, ord$dataset)]
prof[, dens := 0.9 * density / max(density), by = dataset]

# draw the top row first so each lower ridge's peak sits in front of the row above
prof[, grp := factor(dataset, levels = rev(ord$dataset))]
p_a <- ggplot(prof, aes(depth, ymin = y, ymax = y + dens * 1.6, group = grp)) +
  geom_ribbon(fill = "grey80", colour = "grey20", linewidth = 0.2) +
  scale_y_continuous(breaks = seq_len(nrow(ord)) + 0.5, labels = ord$lab,
                     limits = c(0.8, nrow(ord) + 1.5), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Cortical depth of sampled nuclei\n(0 = pia, 1 = white matter)",
       y = "Dataset") +
  theme_cns + theme(axis.ticks.y = element_blank())

# ---- b: same donor, two studies ---------------------------------------------
pr <- fread(file.path(DATA, "fig3_pair_donor.tsv"))
pr[, pair := factor(pair_lab(pair), levels = pair_lab(names(PAIR_COL)))]
PCOL <- setNames(PAIR_COL, pair_lab(names(PAIR_COL)))
st_b <- pr[is.finite(depth_exc_a) & is.finite(depth_exc_b),
           .(r = cor(depth_exc_a, depth_exc_b), n = .N), by = pair]
st_b[, lab := sprintf("r = %.2f\nn = %d", r, n)]
print(st_b)

pr_b <- pr[is.finite(depth_exc_a) & is.finite(depth_exc_b)]
p_b <- ggplot(pr_b, aes(depth_exc_a, depth_exc_b)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.3,
              linetype = "dashed") +
  geom_point(aes(colour = pair), size = 0.7, alpha = 0.6, stroke = 0,
             show.legend = FALSE) +
  geom_text(data = st_b, aes(label = lab), x = -Inf, y = Inf, hjust = -0.1,
            vjust = 1.2, size = fontsize / .pt, lineheight = 0.9,
            inherit.aes = FALSE) +
  scale_colour_manual(values = PCOL) +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
  scale_y_continuous(breaks = c(0.2, 0.5, 0.8)) +
  coord_cartesian(xlim = c(0.1, 0.9), ylim = c(0.1, 0.9)) +
  facet_wrap(~ pair, nrow = 2, axes = "all", axis.labels = "margins", labeller = as_labeller(function(x)
    sub("BrainSCOPE CMC vs\n", "BrainSCOPE\nCMC vs\n", sub(" vs ", " vs\n", x)))) +
  labs(x = "Sampled depth, first study", y = "Sampled depth of the\nsame donor, second study") +
  theme_cns

# ---- c: depth gap vs agreement ----------------------------------------------
pe <- pr[is.finite(depth_exc_a) & is.finite(depth_exc_b) & is.finite(agree_inh)]
pe[, gap := abs(depth_exc_a - depth_exc_b)]
pe[, rk := frank(gap) / .N, by = pair]
pool <- spear(pe$rk, pe$agree_inh)
print(pe[, {t <- spear(gap, agree_inh); .(rho = t$estimate, p = t$p.value, n = .N)}, by = pair])

p_c <- ggplot(pe, aes(gap, agree_inh, colour = pair)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(size = 0.6, alpha = 0.45, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
           size = fontsize / .pt,
           label = sprintf("ρ = %.2f, p = %s\nn = %d donors",
                           pool$estimate, fmt_p(pool$p.value), nrow(pe))) +
  scale_colour_manual(values = PCOL, name = "Study pair") +
  # headroom above the data so the legend sits over empty space
  scale_y_continuous(breaks = c(-0.4, 0, 0.4, 0.8), limits = c(NA, 1.15)) +
  labs(x = "Difference in sampled depth between\nthe donor's two samples",
       y = "Agreement of the donor's inhibitory-neuron\ncomposition between studies (r)") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.8, alpha = 1, linewidth = 0.8))) +
  theme_cns +
  # top-right corner, where donors with large depth gaps are few
  theme(legend.position = "inside", legend.position.inside = c(1, 1),
        legend.justification = c(1, 1), legend.key.height = unit(7, "pt"),
        legend.key.width = unit(10, "pt"), legend.spacing.y = unit(1, "pt"),
        legend.background = element_rect(fill = "white", colour = NA))

design <- "
AB
AC
"
fig <- wrap_plots(A = free(p_a, type = "panel", side = "t"), B = p_b, C = p_c,
                  design = design) +
  plot_layout(widths = c(1, 1.1), heights = c(0.85, 1)) +
  plot_annotation(tag_levels = "a")

f <- file.path(FIG_DIR, "figure3_part1_draft")
ggsave(paste0(f, ".png"), fig, width = 7.4, height = 8, dpi = 300, bg = "white")
ggsave(paste0(f, ".pdf"), fig, width = 7.4, height = 8, bg = "white", device = cairo_pdf)
message("wrote ", f, ".{png,pdf}")
