# Shared theme and constants for the pipeline comparison figures.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

CMP_DIR  <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow"
DATA_DIR <- file.path(CMP_DIR, "data")
FIG_DIR  <- file.path(CMP_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# One colour per pipeline, used consistently across all six figures.
PIPE_COL <- c(mine = "#C1553B", his = "#2E6F8E")
PIPE_LAB <- c(mine = "Nextflow, 19 MGP classes, meta",
              his  = "ctp-gwas phase 2, 33 supertypes, mega")

# SEA-AD class colours, for the supertype fan-out panels.
CLASS_COL <- c(Excitatory = "#4C72B0", Inhibitory = "#DD8452", Glia = "#55A868",
               `NA` = "grey70")

GW_LINE   <- 5e-8      # conventional genome-wide threshold, used by mine
MEFF_LINE <- 2.94e-09  # Shreejoy's Li-and-Ji corrected threshold (REPORT.md)

theme_cmp <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(face = "bold", size = rel(0.9)),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(colour = "grey30", size = rel(0.88)),
      plot.caption = element_text(colour = "grey40", size = rel(0.72),
                                  hjust = 0),
      legend.key.size = unit(0.9, "lines")
    )
}

save_fig <- function(plot, name, width, height, dpi = 200) {
  path <- file.path(FIG_DIR, name)
  ggsave(path, plot, width = width, height = height, dpi = dpi,
         units = "in", bg = "white")
  message(sprintf("  wrote %s (%.0f x %.0f in)", path, width, height))
}

read_stage <- function(name) fread(file.path(DATA_DIR, name), sep = "\t")
