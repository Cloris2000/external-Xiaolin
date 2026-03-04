#!/usr/bin/env Rscript
# Generate side-by-side Manhattan plots for all 19 cell types:
#   Pipeline meta-analysis  vs  Individual script meta-analysis
#
# Output: manhattan_plots/pipeline_vs_ind/{CELL_TYPE}_pipeline_vs_ind.png
# Usage: Rscript scripts/generate_all_cell_type_manhattan.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(gridExtra)
  library(grid)
})

# ── Paths ──────────────────────────────────────────────────────────────────

PIPELINE_DIR <- "results/meta_analysis"
IND_DIR      <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/METAL/Joint_AMP_AD_meta_with_CMC_NABEC_GTEx"
OUT_DIR      <- "manhattan_plots/pipeline_vs_ind"

# Pipeline uses: {CT}_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP.tbl from results/meta_analysis
# Individual uses: {CT}_meta_analysis1.tbl from /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/METAL/Joint_AMP_AD_meta_with_CMC_NABEC_GTEx

CELL_TYPES <- c("Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3",
                "L5.6.NP", "L5.ET", "L6.CT", "L6b", "LAMP5",
                "Microglia", "OPC", "Oligodendrocyte", "PAX6", "PVALB",
                "Pericyte", "SST", "VIP", "VLMC")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ────────────────────────────────────────────────────────────────

parse_marker <- function(marker) {
  # format: chr1:12345:A:G
  parts <- strsplit(as.character(marker), ":")
  chrom <- as.integer(gsub("chr", "", sapply(parts, `[`, 1)))
  pos   <- as.integer(sapply(parts, `[`, 2))
  data.frame(CHROM = chrom, GENPOS = pos)
}

load_metal <- function(path, label, max_bg = 100000) {
  cat("  [", label, "] Reading:", path, "\n")
  df <- fread(path, data.table = FALSE)
  # parse chromosome and position from MarkerName
  coords <- parse_marker(df$MarkerName)
  df <- cbind(df, coords)
  df <- df %>% filter(!is.na(CHROM), !is.na(`P-value`), `P-value` > 0)
  colnames(df)[colnames(df) == "P-value"] <- "P"

  n_total <- nrow(df)
  n_sug   <- sum(df$P < 1e-5, na.rm = TRUE)
  n_sig   <- sum(df$P < 5e-8, na.rm = TRUE)
  lambda  <- median(qchisq(1 - df$P, 1), na.rm = TRUE) / qchisq(0.5, 1)
  cat("    Variants:", n_total, "| Suggestive:", n_sug, "| GW-sig:", n_sig,
      "| Lambda:", round(lambda, 4), "\n")

  # downsample background
  sig <- df %>% filter(P < 1e-5)
  bg  <- df %>% filter(P >= 1e-5)
  if (nrow(bg) > max_bg) bg <- bg[sample(nrow(bg), max_bg), ]
  df_plot <- bind_rows(sig, bg)

  list(df = df_plot, n = n_total, n_sug = n_sug, n_sig = n_sig, lambda = lambda)
}

make_manhattan <- function(df, title, subtitle = NULL,
                           suggestive = 1e-5, gwsig = 5e-8,
                           point_color = c("#444444", "#AAAAAA")) {
  # compute cumulative positions
  chr_info <- df %>%
    group_by(CHROM) %>%
    summarise(chr_max = max(GENPOS), .groups = "drop") %>%
    arrange(CHROM) %>%
    mutate(offset = cumsum(as.numeric(lag(chr_max, default = 0))))

  df <- df %>%
    left_join(chr_info %>% select(CHROM, offset), by = "CHROM") %>%
    mutate(BPcum = GENPOS + offset)

  axis_df <- df %>%
    group_by(CHROM) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")

  p <- ggplot(df, aes(x = BPcum, y = -log10(P),
                      color = as.factor(CHROM %% 2))) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = -log10(gwsig),      color = "#E41A1C", linetype = "dashed",  linewidth = 0.7) +
    geom_hline(yintercept = -log10(suggestive), color = "#377EB8", linetype = "dotted",  linewidth = 0.5) +
    scale_color_manual(values = setNames(point_color, c("0","1")), guide = "none") +
    scale_x_continuous(labels = axis_df$CHROM, breaks = axis_df$center) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title    = title,
         subtitle = subtitle,
         x        = "Chromosome",
         y        = expression(-log[10](p))) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      axis.line          = element_line(linewidth = 0.4),
      plot.title         = element_text(size = 10, face = "bold"),
      plot.subtitle      = element_text(size = 8, color = "grey40")
    )
  p
}

# ── Main loop ──────────────────────────────────────────────────────────────

for (ct in CELL_TYPES) {
  cat("\n==============================\n")
  cat("Cell type:", ct, "\n")
  cat("==============================\n")

  pipe_file <- file.path(PIPELINE_DIR,
    paste0(ct, "_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP.tbl"))
  ind_file  <- file.path(IND_DIR, paste0(ct, "_meta_analysis1.tbl"))

  if (!file.exists(pipe_file)) {
    cat("  [SKIP] Pipeline file not found:", pipe_file, "\n"); next
  }
  if (!file.exists(ind_file)) {
    cat("  [SKIP] Individual file not found:", ind_file, "\n"); next
  }

  pipe_res <- load_metal(pipe_file, "Pipeline")
  ind_res  <- load_metal(ind_file,  "Individual")

  pipe_sub <- paste0("N=", pipe_res$n, " | Sug=", pipe_res$n_sug,
                     " | GW-sig=", pipe_res$n_sig, " | λ=", round(pipe_res$lambda, 3))
  ind_sub  <- paste0("N=", ind_res$n,  " | Sug=", ind_res$n_sug,
                     " | GW-sig=", ind_res$n_sig,  " | λ=", round(ind_res$lambda, 3))

  p_pipe <- make_manhattan(pipe_res$df,
    title    = paste0(ct, "  —  Pipeline (8 cohorts meta-analysis)"),
    subtitle = pipe_sub,
    point_color = c("#1B4F72", "#85C1E9"))   # blue tones for pipeline

  p_ind  <- make_manhattan(ind_res$df,
    title    = paste0(ct, "  —  Individual Scripts (meta-analysis)"),
    subtitle = ind_sub,
    point_color = c("#6E2F1A", "#F0A87A"))   # orange tones for individual

  # stack vertically: pipeline on top, individual on bottom
  out_png <- file.path(OUT_DIR, paste0(ct, "_pipeline_vs_ind.png"))
  png(out_png, width = 2400, height = 900, res = 150)
  grid.arrange(p_pipe, p_ind, nrow = 2,
               top = textGrob(paste(ct, ": Pipeline vs Individual Scripts"),
                              gp = gpar(fontsize = 13, fontface = "bold")))
  dev.off()
  cat("  Saved:", out_png, "\n")
}

cat("\n\nDone! Plots saved in:", OUT_DIR, "\n")
cat("Files:\n")
cat(paste(list.files(OUT_DIR, pattern = "\\.png$"), collapse = "\n"), "\n")
