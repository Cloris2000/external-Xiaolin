# Figure 3. Inconsistent sampling depth between tissue blocks lowers agreement
# between datasets, and AD pathology looks like deeper sampling, so depth cannot
# simply be covaried.
#   a  depth distribution of sampled nuclei, per dataset
#   b  the same donor sampled by two studies: sampled depth in each
#   c  depth gap between a donor's two blocks vs agreement of that donor's
#      inhibitory-neuron composition (depth from excitatory neurons only, so the
#      two measures use different cells)
#   d  Green: per-supertype change in composition when a block is sampled deeper
#      (donors without dementia) and per-supertype CERAD effect, both against
#      the supertype's cortical depth
#   e  the depth-disease-effect relationship (Spearman rho across supertypes)
#      before and after adjusting for sampled depth, all datasets and phenotypes
#
# The share of the CERAD gradient depth alone could produce (~60%, 95% CI
# 35-170%, Green) is reported in text from fig3_part2_summary.tsv, not plotted.
#
# Taxonomy: SEA-AD DFC 2026 (reactive states included; 15 reactive states
# without a spatial depth reference are excluded from depth-based panels).
# Inputs: fig3_data.py (a-c), fig3_part2_data.py (d-e).
# Colours: cnsplots "Nature" blue = neuronal (d); cnsplots "Cell" for study
# pairs (b, c) and for model state (e), no colour used with two meanings.
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig3.R

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
LABEL <- c(SEAAD_DFC_2026 = "SEA-AD DLPFC", Multiome2025_DLPFC = "Multiome 2025",
           BrainSCOPE_CMC = "BrainSCOPE CMC", BrainSCOPE_UCLA = "BrainSCOPE UCLA",
           BrainSCOPE_Multiome = "BrainSCOPE Multiome")
pretty <- function(x) ifelse(x %in% names(LABEL), LABEL[x], gsub("_", " ", x))
PHENO_LAB <- c(cerad_ad = "CERAD", cogdx = "Cognitive diagnosis", dementia = "Dementia")

# ---- palette (cnsplots hex values) --------------------------------------------
# Nature: neuronal classes
CLS_COL <- c(Excitatory = "#3C5488", Inhibitory = "#8491B4")
# Cell: study pairs (b, c)
PAIR_COL <- c("Green x Mathys" = "#2F7E8F", "Green x Multiome" = "#E1A22E",
              "Mathys x Multiome" = "#8B6FA8", "Multiome x PsychAD" = "#C84C3A",
              "BrainSCOPE CMC x PsychAD" = "#5F9862", "PsychAD x Ruzicka" = "#4E5A8A")
# Cell: model state (e), the two Cell colours not used by any study pair
MODEL_COL <- c(Unadjusted = "#B85F7A", `Depth-adjusted` = "#7B8C9E")
pair_lab <- function(x) gsub(" x ", " vs ", x)

fontsize <- 8
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

# ==== Part 1: sampling depth differs between blocks ============================

# ---- a: depth distribution per dataset ----------------------------------------
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

# ---- b: same donor, two studies -------------------------------------------------
pr <- fread(file.path(DATA, "fig3_pair_donor.tsv"))
pr[, pair := factor(pair_lab(pair), levels = pair_lab(names(PAIR_COL)))]
PCOL <- setNames(PAIR_COL, pair_lab(names(PAIR_COL)))
pr_b <- pr[is.finite(depth_exc_a) & is.finite(depth_exc_b)]
st_b <- pr_b[, .(r = cor(depth_exc_a, depth_exc_b), n = .N), by = pair]
st_b[, lab := sprintf("r = %.2f\nn = %d", r, n)]
cat("b: same donor, two studies\n"); print(st_b)

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
  facet_wrap(~ pair, nrow = 2, axes = "all", axis.labels = "margins",
             labeller = as_labeller(function(x)
               sub("BrainSCOPE CMC vs\n", "BrainSCOPE\nCMC vs\n", sub(" vs ", " vs\n", x)))) +
  labs(x = "Sampled depth, first study",
       y = "Sampled depth of the\nsame donor, second study") +
  theme_cns

# ---- c: depth gap vs agreement --------------------------------------------------
pe <- pr[is.finite(depth_exc_a) & is.finite(depth_exc_b) & is.finite(agree_inh)]
pe[, gap := abs(depth_exc_a - depth_exc_b)]
pe[, rk := frank(gap) / .N, by = pair]
pool <- spear(pe$rk, pe$agree_inh)
cat("\nc: depth gap vs inhibitory agreement, per pair\n")
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
  theme(legend.position = "inside", legend.position.inside = c(1, 1),
        legend.justification = c(1, 1), legend.key.height = unit(7, "pt"),
        legend.key.width = unit(10, "pt"), legend.spacing.y = unit(1, "pt"),
        legend.background = element_rect(fill = "white", colour = NA))

# ==== Part 2: AD looks like deeper sampling ====================================
eff <- fread(file.path(DATA, "fig3_part2_effects.tsv"))
sm <- fread(file.path(DATA, "fig3_part2_summary.tsv"))

# ---- d: the two patterns, Green x CERAD -----------------------------------------
g <- eff[dataset == "Green" & phenotype == "cerad_ad"]
long <- rbind(
  g[, .(cell_type, cls, depth, value = gamma, what = "Sampled deeper")],
  g[, .(cell_type, cls, depth, value = beta, what = "Higher CERAD")])
long[, what := factor(what, levels = c("Sampled deeper", "Higher CERAD"))]
st_d <- long[, {t <- spear(depth, value)
                .(lab = sprintf("ρ = %.2f\np = %s", t$estimate, fmt_p(t$p.value)))}, by = what]
cat("\nd: Green x CERAD\n"); print(st_d)

p_d <- ggplot(long, aes(depth, value)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey25", fill = "grey85",
              linewidth = 0.4) +
  geom_point(aes(colour = cls), size = 1, alpha = 0.85, stroke = 0) +
  geom_text(data = st_d, aes(label = lab), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            size = fontsize / .pt, lineheight = 0.9, inherit.aes = FALSE) +
  scale_colour_manual(values = CLS_COL) +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
  facet_wrap(~ what, nrow = 1, scales = "free_y", axes = "all") +
  labs(x = "Cortical depth of the supertype (0 = pia, 1 = white matter)",
       y = "Change in supertype CLR (Δ CLR)") +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  theme_cns +
  theme(legend.title = element_blank(),
        legend.position = "inside", legend.position.inside = c(1, 0),
        legend.justification = c(1, 0))

# ---- e: depth-disease-effect rho before/after depth adjustment ------------------
sm[, row := sprintf("%s, %s", dataset, PHENO_LAB[phenotype])]
sm[, main := dataset == "Green" & phenotype == "cerad_ad"]
ord_e <- sm[order(factor(dataset, levels = c("Green", "Mathys", "Multiome 2025", "PsychAD")),
                  factor(phenotype, levels = names(PHENO_LAB)))]$row
sm[, row := factor(row, levels = rev(ord_e))]
cat("\ne: rho before/after depth adjustment\n")
print(sm[, .(row, rho_obs = round(rho_obs, 2), rho_adj = round(rho_adj, 2),
             share = round(100 * share), lo = round(100 * share_lo), hi = round(100 * share_hi))])
cat(sprintf("\nGreen x CERAD: depth shift alone could produce %.0f%% (95%% CI %.0f-%.0f%%) of the\n",
            100 * sm[main == TRUE, share], 100 * sm[main == TRUE, share_lo],
            100 * sm[main == TRUE, share_hi]),
    "observed CERAD gradient (text statement; see fig3_part2_summary.tsv).\n")

e_dat <- melt(sm, id.vars = c("row", "main"),
              measure.vars = c("rho_obs", "rho_adj"), variable.name = "model", value.name = "rho")
e_dat[, model := factor(ifelse(model == "rho_obs", "Unadjusted", "Depth-adjusted"),
                        levels = names(MODEL_COL))]

p_e <- ggplot(e_dat, aes(rho, row)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(aes(group = row), colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = model), size = 1.9, stroke = 0) +
  scale_colour_manual(values = MODEL_COL) +
  scale_y_discrete(limits = levels(sm$row)) +
  labs(x = "Depth-disease-effect relationship\n(ρ, supertype abundance)", y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.2))) +
  theme_cns +
  theme(legend.title = element_blank(), legend.position = "top",
        legend.location = "plot", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0), legend.key.spacing.x = unit(6, "pt"))

# ==== layout ====================================================================
# free() stops patchwork aligning panel areas across the grid; otherwise e's
# long row labels shrink b and c's plot areas, and a's ridge labels shrink d's.
# A 20-column grid lets the bottom row use a different column split from the
# top rows (e needs more width than b/c because of its row labels) without
# nesting, which free() does not survive.
design <- "
AAAAAAAAABBBBBBBBBBB
AAAAAAAAACCCCCCCCCCC
DDDDDDDDDEEEEEEEEEEE
"
fig <- wrap_plots(A = free(p_a, type = "panel"), B = free(p_b, type = "panel"),
                  C = free(p_c, type = "panel"), D = free(p_d, type = "panel"),
                  E = free(p_e, type = "panel"), design = design) +
  plot_layout(heights = c(0.8, 1, 1.05)) +
  plot_annotation(tag_levels = "a")

# 183 mm = Nature double-column width
W_IN <- 183 / 25.4
H_IN <- 9.4
f <- file.path(FIG_DIR, "figure3_sampling_depth")
ggsave(paste0(f, ".pdf"), fig, width = W_IN, height = H_IN, bg = "white", device = cairo_pdf)
ggsave(paste0(f, ".png"), fig, width = W_IN, height = H_IN, dpi = 300, bg = "white")
message("wrote ", f, ".{pdf,png}")
