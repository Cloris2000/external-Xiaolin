# Shared ggplot theme and palette for the CTP reliability paper (Figs 4-6).
# Figs 2-3 carry inline copies of the same values; keep these in sync.
# Colours are cnsplots hex values. One meaning per colour across the paper:
#   Nature blue/red  = neuronal / non-neuronal (cell class)
#   Cell palette     = study pairs (Fig 3), model state (Fig 3e)
#   ACCENT (Nature teal) = the arm being argued for (consensus / production)
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

fontsize <- 8
NEU <- "#3C5488"; GLIA <- "#E64B35"; ACCENT <- "#00A087"
CLASS_COL <- c("Neuronal" = NEU, "Non-neuronal" = GLIA, "Reactive (-SEAAD)" = GLIA)
CLASS_OUT <- c("Neuronal" = NEU, "Non-neuronal" = GLIA, "Reactive (-SEAAD)" = "#111111")
CLASS_SHP <- c("Neuronal" = 21, "Non-neuronal" = 21, "Reactive (-SEAAD)" = 23)
GREYS <- c("#B7B7B7", "#7B8C9E", "#3D3D3D")

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
class_key <- function(compartment, reactive) {
  factor(fcase(reactive, "Reactive (-SEAAD)", compartment == "neu", "Neuronal",
               default = "Non-neuronal"), levels = names(CLASS_COL))
}
save_fig <- function(fig, path, height_in, width_mm = 183) {
  ggsave(paste0(path, ".pdf"), fig, width = width_mm / 25.4, height = height_in,
         bg = "white", device = cairo_pdf)
  ggsave(paste0(path, ".png"), fig, width = width_mm / 25.4, height = height_in,
         dpi = 300, bg = "white")
  message("wrote ", path, ".{pdf,png}")
}
