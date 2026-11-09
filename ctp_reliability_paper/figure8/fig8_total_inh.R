# Figure 8 variant: y = total inhibitory share of neurons (one number per
# donor), instead of the mean over inhibitory supertypes. See
# fig8_total_inh_data.py. Same donors and x as the main figure.
#   (a) per-pair scatter: x = cortdepth Wasserstein depth discordance,
#       y = |residual logit(inh/neurons)_A - _B|
#   (b) per pair, Spearman rho (with Fisher-z 95% CI) of depth discordance
#       against the two y definitions side by side, plus the pooled
#       fixed-effects partial R2 of each
# Inputs: data/fig8_total_inh.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig8_total_inh.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure8"
source(file.path(dirname(HERE), "common", "theme_cns.R"))
FIG_DIR <- file.path(HERE, "figures")

PAIRS <- c("Green x Mathys", "Green x Multiome", "Mathys x Multiome",
          "Multiome x PsychAD", "BrainSCOPE CMC x PsychAD")
# Figure 3's pair palette (cnsplots "Cell"): one colour per pair across the paper
PAIR_COL <- c("Green x Mathys" = "#2F7E8F", "Green x Multiome" = "#E1A22E",
             "Mathys x Multiome" = "#8B6FA8", "Multiome x PsychAD" = "#C84C3A",
             "BrainSCOPE CMC x PsychAD" = "#5F9862")
cohorts <- fread(file.path(dirname(HERE), "common", "cohort_names.tsv"))
canon <- function(key) {
  nm <- cohorts$canonical_name[cohorts$source_key == key]
  stopifnot(length(nm) == 1)
  nm
}
SIDE_KEY <- c("Green" = "Green_DFC_person", "Mathys" = "Mathys_DFC_person",
              "Multiome" = "Multiome2025_DLPFC_DFC_person", "PsychAD" = "PsychAD_full_DFC_person",
              "BrainSCOPE CMC" = "BrainSCOPE_CMC_DFC_person")
wrap_name <- function(nm, width = 24) {
  if (nchar(nm) <= width) return(nm)
  cut <- max(gregexpr("_", substr(nm, 1, width))[[1]])
  if (cut < 1) return(nm)
  paste0(substr(nm, 1, cut), "\n", substr(nm, cut + 1, nchar(nm)))
}
STRIP_LAB <- setNames(vapply(strsplit(PAIRS, " x "), function(s)
  paste0(wrap_name(canon(SIDE_KEY[[s[1]]])), "\n\u00d7 ", wrap_name(canon(SIDE_KEY[[s[2]]]))), ""), PAIRS)
ROW_LAB <- setNames(vapply(strsplit(PAIRS, " x "), function(s)
  paste0(canon(SIDE_KEY[[s[1]]]), "\n\u00d7 ", canon(SIDE_KEY[[s[2]]])), ""), PAIRS)

Y_DEFS <- c(disc_inh_total = "Total inhibitory share of neurons",
            disc_inh = "Mean over 71 inhibitory supertypes")
Y_COL <- c("Total inhibitory share of neurons" = "#111111",
           "Mean over 71 inhibitory supertypes" = "#9A9A9A")

d <- fread(file.path(HERE, "data/fig8_total_inh.tsv"))[pass == TRUE]
d[, pair := factor(pair, levels = PAIRS)]
cat(sprintf("donors used: %d\n", nrow(d)))

# ---- pooled fixed-effects, both y definitions ----------------------------
pooled <- function(y) {
  f <- lm(reformulate(c("wasserstein", "pair"), y), d); n0 <- lm(reformulate("pair", y), d)
  cf <- summary(f)$coefficients["wasserstein", ]
  data.table(y = y, beta = cf["Estimate"], p = cf["Pr(>|t|)"],
             partial_r2 = 1 - sum(residuals(f)^2) / sum(residuals(n0)^2))
}
pool <- rbindlist(lapply(names(Y_DEFS), pooled))
cat("\npooled fixed-effects (pair as covariate):\n")
print(pool[, .(y, beta = round(beta, 3), p = fmt_p(p), partial_r2 = sprintf("%.1f%%", 100 * partial_r2))])

# ---- per-pair Spearman with Fisher-z CI, both y ---------------------------
rho_ci <- function(x, y) {
  s <- spear(x, y); r <- unname(s$estimate); n <- length(x)
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  list(n = n, rho = r, lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se), p = s$p.value)
}
stats <- rbindlist(lapply(names(Y_DEFS), function(y)
  d[, c(list(y = y), rho_ci(wasserstein, get(y))), by = pair]))
stats[, ydef := factor(Y_DEFS[y], levels = Y_DEFS)]
cat("\nper pair, Spearman rho of depth discordance vs y:\n")
print(dcast(stats, pair ~ ydef, value.var = "rho")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])

# ---- (a) scatter, total-inhibitory y --------------------------------------
lab <- stats[y == "disc_inh_total", .(pair, lab = sprintf("n=%d\nrho=%+.2f, p=%s", n, rho, fmt_p(p)))]
lab[, `:=`(x = min(d$wasserstein), yy = max(d$disc_inh_total))]
p_a <- ggplot(d, aes(wasserstein, disc_inh_total)) +
  geom_point(aes(colour = pair), size = 1.1, alpha = 0.75, stroke = 0) +
  geom_smooth(aes(colour = pair), method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = lab, aes(x, yy, label = lab), hjust = 0, vjust = 1,
           size = (fontsize - 1) / ggplot2::.pt, colour = "grey30", lineheight = 0.9) +
  scale_colour_manual(values = PAIR_COL, name = NULL) +
  facet_wrap(~ pair, nrow = 1, labeller = as_labeller(STRIP_LAB)) +
  labs(x = "Dissection-depth discordance between the two samples\n(Wasserstein distance of excitatory depth profiles, cortex units)",
       y = "Cross-dataset disagreement in the donor's\ninhibitory share of neurons (|\u0394 residual logit|)",
       subtitle = "One point per donor; y is a single number per donor (inhibitory nuclei / all neuronal nuclei), not a supertype average") +
  theme_cns +
  theme(legend.position = "none",
        strip.text = element_text(size = fontsize - 1.5, lineheight = 0.95),
        plot.subtitle = element_text(size = fontsize - 1, colour = "grey30"))

# ---- (b) rho per pair, two y definitions side by side ----------------------
stats[, row_lab := factor(ROW_LAB[as.character(pair)], levels = rev(unname(ROW_LAB)))]
pool_lab <- pool[, sprintf("%s: pooled partial R\u00b2 = %.1f%% (p = %s)", Y_DEFS[y], 100 * partial_r2, fmt_p(p))]
p_b <- ggplot(stats, aes(rho, row_lab, colour = ydef)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y",
               position = position_dodge(width = 0.5), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  scale_colour_manual(values = Y_COL, name = "y definition") +
  scale_x_continuous(limits = c(-0.35, 1), breaks = seq(-0.25, 1, 0.25)) +
  labs(x = "Spearman rho, depth discordance vs cross-dataset disagreement (95% CI)", y = NULL,
       subtitle = paste(pool_lab, collapse = "\n")) +
  theme_cns +
  theme(axis.text.y = element_text(hjust = 0, lineheight = 0.95),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        plot.subtitle = element_text(size = fontsize - 1, colour = "grey30"))

fig <- (free(p_a, type = "panel") / free(p_b, type = "panel")) +
  plot_layout(heights = c(1.5, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = fontsize, face = "bold", colour = "black"))

save_fig(fig, file.path(FIG_DIR, "figure8_total_inhibitory_check"), height_in = 6.8, width_mm = 240)
