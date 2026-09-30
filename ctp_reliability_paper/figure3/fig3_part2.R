# Figure 3, part 2 (DRAFT). AD looks like deeper sampling, so depth cannot
# simply be covaried.
#   d  Green: per-supertype change in composition when a block is sampled
#      deeper (donors without dementia) and per-supertype CERAD effect, both
#      against the supertype's cortical depth. d alone carries the qualitative
#      claim (both panels slope the same way against cortical depth).
#   e  the AD gradient before and after adjusting for sampled depth
#
# The share of the gradient depth alone could produce (~60%, 95% CI 35-170%
# for Green x CERAD) is reported in text/table, not as its own panel: two
# overlapping trend lines through ~106 points asks readers to compare
# slopes by eye, which tested poorly. Numbers for that statement stay in
# fig3_part2_data.py's summary table.
#
# Inputs from fig3_part2_data.py (design of the PI's 67/68, on DFC 2026).
# Run: /home/zhoux156/miniforge3/envs/test/bin/Rscript fig3_part2.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure3"
DATA <- file.path(HERE, "data")
FIG_DIR <- file.path(HERE, "figures")

fontsize <- 8
# cnsplots "Nature": blue is neuronal across the paper; two shades for the
# two neuronal classes
CLS_COL <- c(Excitatory = "#3C5488", Inhibitory = "#8491B4")
theme_cns <- theme_classic(base_size = fontsize, base_family = "sans") +
  theme(
    text = element_text(size = fontsize, colour = "black"),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.text = element_text(size = fontsize, colour = "black"),
    axis.title = element_text(size = fontsize, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = fontsize, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = fontsize, colour = "black"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(8, "pt"),
    plot.tag = element_text(size = fontsize, face = "bold", colour = "black"),
    plot.margin = margin(2, 4, 2, 4)
  )
fmt_p <- function(p) ifelse(p < 1e-3, sprintf("%.0e", p), sprintf("%.2g", p))
PHENO_LAB <- c(cerad_ad = "CERAD", cogdx = "Cognitive diagnosis", dementia = "Dementia")

eff <- fread(file.path(DATA, "fig3_part2_effects.tsv"))
sm <- fread(file.path(DATA, "fig3_part2_summary.tsv"))

# ---- d: the two patterns, Green x CERAD ---------------------------------------
g <- eff[dataset == "Green" & phenotype == "cerad_ad"]
long <- rbind(
  g[, .(cell_type, cls, depth, value = gamma, what = "Sampled deeper")],
  g[, .(cell_type, cls, depth, value = beta, what = "Higher CERAD")])
long[, what := factor(what, levels = c("Sampled deeper", "Higher CERAD"))]
st_d <- long[, {t <- suppressWarnings(cor.test(depth, value, method = "spearman", exact = FALSE))
                .(lab = sprintf("ρ = %.2f\np = %s", t$estimate, fmt_p(t$p.value)))}, by = what]
print(st_d)

p_d <- ggplot(long, aes(depth, value)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey25", fill = "grey85",
              linewidth = 0.4) +
  geom_point(aes(colour = cls), size = 1, alpha = 0.85, stroke = 0) +
  geom_text(data = st_d, aes(label = lab), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            size = fontsize / .pt, lineheight = 0.9, inherit.aes = FALSE) +
  scale_colour_manual(values = CLS_COL) +
  facet_wrap(~ what, nrow = 1, scales = "free_y", axes = "all") +
  labs(x = "Cortical depth of the supertype (0 = pia, 1 = white matter)",
       y = "Change in supertype CLR (Δ CLR)") +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  theme_cns +
  theme(legend.position = "inside", legend.position.inside = c(1, 0),
        legend.justification = c(1, 0))

# ---- e: the AD gradient before/after depth adjustment, all datasets/phenotypes --
sm[, row := sprintf("%s, %s", dataset, PHENO_LAB[phenotype])]
sm[, main := dataset == "Green" & phenotype == "cerad_ad"]
sm[, no_tilt := p_obs > 0.05]
ord <- sm[order(factor(dataset, levels = c("Green", "Mathys", "Multiome 2025", "PsychAD")),
                factor(phenotype, levels = names(PHENO_LAB)))]$row
sm[, row := factor(row, levels = rev(ord))]
print(sm[, .(row, share = round(100 * share), lo = round(100 * share_lo),
             hi = round(100 * share_hi), rho_obs = round(rho_obs, 2), rho_adj = round(rho_adj, 2))])
cat(sprintf("\nGreen x CERAD: depth shift alone could produce %.0f%% (95%% CI %.0f-%.0f%%) of the\n",
            100 * sm[main == TRUE, share], 100 * sm[main == TRUE, share_lo],
            100 * sm[main == TRUE, share_hi]),
    "observed CERAD gradient (text statement, not plotted; see fig3_part2_summary.tsv).\n")

f_dat <- melt(sm, id.vars = c("row", "main", "no_tilt"),
              measure.vars = c("rho_obs", "rho_adj"), variable.name = "model", value.name = "rho")
f_dat[, model := factor(ifelse(model == "rho_obs", "Unadjusted", "Depth-adjusted"),
                        levels = c("Unadjusted", "Depth-adjusted"))]
p_e <- ggplot(f_dat, aes(rho, row)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(aes(group = row), colour = "grey55", linewidth = 0.4) +
  geom_point(aes(shape = model), size = 1.8, stroke = 0.4, fill = "white") +
  scale_shape_manual(values = c(Unadjusted = 16, `Depth-adjusted` = 21)) +
  scale_y_discrete(limits = levels(sm$row)) +
  labs(x = "Depth-disease-effect relationship\n(ρ, supertype abundance)", y = NULL) +
  theme_cns +
  theme(legend.position = "top",
        legend.justification = "left", legend.margin = margin(0, 0, 0, 0))

design <- "
D
E
"
fig <- wrap_plots(D = free(p_d, type = "panel"), E = p_e, design = design) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = list(c("d", "e")))

f <- file.path(FIG_DIR, "figure3_part2_draft")
ggsave(paste0(f, ".png"), fig, width = 7.4, height = 6.2, dpi = 300, bg = "white")
ggsave(paste0(f, ".pdf"), fig, width = 7.4, height = 6.2, bg = "white", device = cairo_pdf)
message("wrote ", f, ".{png,pdf}")
