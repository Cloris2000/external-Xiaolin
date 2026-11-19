# Figure 10 supplement (DRAFT): Sst_25 vs Braak stage, Mount Sinai (MSBB) cohorts.
#   a  all AD-continuum donors (Control, MCI, AD) per modality: bulk MSBB_FP, snRNA PsychAD_MSBB
#   b  the donors measured by both
# Boxes = Sst_25 z-scored within modality (CLR scale; bulk = CelMod prediction). p = linear
# Braak trend from z ~ Braak + age + sex + PMI (fig10_braak_data.py). Points with |z| > 3.5
# drawn at the axis limit as open triangles (display only; all donors in the models).
# Inputs: data/fig10_braak_donors.tsv, data/fig10_braak_effects.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig10_braak.R
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(HERE, "..", "common", "theme_cns.R"))

MOD_COL <- c(bulk = "#E69F00", snRNA = "#0072B2")   # as figure10_supp_sst25 c, d
MOD_LAB <- c(bulk = "Bulk RNA-seq (deconvolved), MSBB_FP", snRNA = "snRNA-seq, PsychAD_MSBB")
YCAP <- 3.5
fmtp <- function(p) ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p))

dn <- fread(file.path(HERE, "data", "fig10_braak_donors.tsv"))
ef <- fread(file.path(HERE, "data", "fig10_braak_effects.tsv"))
dn[, mod := factor(modality, levels = names(MOD_COL))]
dn[, stage := factor(braak, levels = 0:6)]
dn[, capped := abs(z) > YCAP][, yplot := pmax(pmin(z, YCAP), -YCAP)]
cat("donors per Braak stage:\n"); print(dcast(dn, set + modality ~ braak, fun.aggregate = length))

panel <- function(s, title) {
  x <- dn[set == s]; e <- ef[set == s & modality %in% names(MOD_COL)]
  lab <- paste0(c(bulk = "Bulk", snRNA = "snRNA")[e$modality], ": β = ", sprintf("%.3f", e$beta),
                " SD per stage, p = ", fmtp(e$p), collapse = "\n")
  d <- ef[set == s & modality == "bulk_minus_snRNA"]
  if (nrow(d)) lab <- paste0(lab, "\nBulk vs snRNA slope: p = ", fmtp(d$p))
  nn <- x[, .(n = uniqueN(person_id)), by = modality]
  ntxt <- if (s == "same_donor") paste0("n = ", nn$n[1]) else
    paste0("bulk n = ", nn[modality == "bulk"]$n, ", snRNA n = ", nn[modality == "snRNA"]$n)
  ggplot(x, aes(x = stage, y = yplot)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60", linetype = "dashed") +
    geom_boxplot(aes(fill = mod), outlier.shape = NA, width = 0.75, linewidth = 0.3, alpha = 0.7,
                 position = position_dodge(0.85)) +
    geom_point(data = x[capped == FALSE], aes(fill = mod), show.legend = FALSE,
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85, seed = 1),
               size = 0.25, alpha = 0.3) +
    geom_point(data = x[capped == TRUE], aes(fill = mod), shape = 2, size = 0.9,
               position = position_dodge(0.85), show.legend = FALSE) +
    annotate("text", x = 0.5, y = YCAP + 0.3, label = lab, hjust = 0, vjust = 0, size = 2.2,
             lineheight = 0.95) +
    scale_fill_manual(values = MOD_COL, labels = MOD_LAB, name = NULL, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = seq(-3, 3, 1)) +
    coord_cartesian(ylim = c(-YCAP, YCAP + 1.6), clip = "off") +
    labs(x = "Braak stage", y = "Sst_25 (z-scored CLR)",
         title = paste0(title, "\n(", ntxt, ")")) +
    theme_cns +
    theme(plot.title = element_text(size = fontsize, face = "bold"), legend.position = "bottom")
}

fig <- (panel("all", "All donors per modality") | panel("same_donor", "Same donors, both modalities")) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
save_fig(fig, file.path(HERE, "figures", "figure10_supp_braak_draft"), height_in = 3.4)
