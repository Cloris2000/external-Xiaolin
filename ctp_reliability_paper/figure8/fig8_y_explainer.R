# Explainer for the Figure 8 y-axis ("% of inhibitory neurons that differ"), for
# discussion with the PI, not a paper panel (2026-11-03). Two real Green x Mathys
# donors (fig8_y_explainer_data.py): one sampled at the same depth in both
# datasets, one at different depths. Each bar is the donor's inhibitory neurons in
# one dataset, split by supertype (thin white lines) and coloured by subclass;
# the % printed for each donor is the Figure 8 y value.
# Input: data/fig8_y_explainer.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig8_y_explainer.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure8"
source(file.path(dirname(HERE), "common", "theme_cns.R"))

# CGE-derived subclasses first, then MGE-derived; no laminar order is implied
SUBCLASS <- c("Lamp5", "Sncg", "Vip", "Pax6", "Lamp5 Lhx6", "Sst", "Sst Chodl",
              "Pvalb", "Chandelier")
SUB_COL <- setNames(c("#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#882255",
                      "#56B4E9", "#999999", "#0072B2", "#009E73"), SUBCLASS)
CASES <- c("Same sampled depth", "Different sampled depth")

x <- fread(file.path(HERE, "data/fig8_y_explainer.tsv"))
stopifnot(all(x$subclass %in% SUBCLASS))
x[, subclass := factor(subclass, levels = SUBCLASS)]
setorder(x, case, dataset, subclass, supertype)
x[, `:=`(xmax = cumsum(share)), by = .(case, dataset)]
x[, xmin := xmax - share]
x[, dataset := factor(dataset, levels = c("Mathys", "Green"))]   # Green on top

lab <- unique(x[, .(case, donor, tv_inh)])
lab[, strip := sprintf("%s (%s): %.0f%% of inhibitory neurons differ", case, donor, tv_inh)]
x <- merge(x, lab[, .(case, strip)], by = "case")
x[, strip := factor(strip, levels = lab[match(CASES, case), strip])]
depth <- unique(x[, .(strip, dataset, upper_frac)])
depth[, lab := sprintf("upper-layer share %.2f", upper_frac)]

p <- ggplot(x) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = as.integer(dataset) - 0.38,
                ymax = as.integer(dataset) + 0.38, fill = subclass),
            colour = "white", linewidth = 0.15) +
  geom_text(data = depth, aes(x = 101.5, y = as.integer(dataset), label = lab),
            hjust = 0, size = (fontsize - 1) / .pt, colour = "grey30") +
  scale_fill_manual(values = SUB_COL, name = NULL) +
  scale_y_continuous(breaks = 1:2, labels = c("Mathys", "Green"), expand = expansion(add = 0.1)) +
  scale_x_continuous(breaks = seq(0, 100, 25), labels = function(v) paste0(v, "%"),
                     expand = expansion(add = c(0, 0))) +
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  facet_wrap(~ strip, ncol = 1) +
  labs(x = "Share of the donor's inhibitory neurons", y = NULL,
       subtitle = paste0("% that differ = half the sum, over inhibitory supertypes, of |share in Green - share in Mathys|",
                         "\n(upper-layer share = L2/3 IT / (L2/3 IT + L6 IT) among the sample's excitatory neurons)")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_cns +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", hjust = 0),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = fontsize - 1, colour = "grey30"),
        plot.margin = margin(4, 100, 4, 4))   # room for the upper-layer share labels

save_fig(p, file.path(HERE, "figures", "figure8_y_explainer"), height_in = 3.4, width_mm = 200)
