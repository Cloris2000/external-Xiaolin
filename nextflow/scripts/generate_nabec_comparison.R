#!/usr/bin/env Rscript
# Generate NABEC pipeline vs individual scripts comparison Manhattan plots
# for SST, VIP, and L5ET cell types

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(topr)
  library(cowplot)
})

PIPELINE_DIR <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/NABEC/regenie_step2"
IND_DIR      <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/NABEC_wgs_step2"
OUT_DIR      <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/manhattan_plots/NABEC_pipeline_vs_ind"
DOWNSAMPLE_N <- 100000
SUGGESTIVE   <- 1e-5

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cell_types <- list(
  list(label = "SST",
       pipeline = file.path(PIPELINE_DIR, "NABEC_SST_step2.regenie.raw_p"),
       ind      = file.path(IND_DIR,      "NABEC_WGS_step2_update_hg19_SST.regenie.raw_p")),
  list(label = "VIP",
       pipeline = file.path(PIPELINE_DIR, "NABEC_VIP_step2.regenie.raw_p"),
       ind      = file.path(IND_DIR,      "NABEC_WGS_step2_update_hg19_VIP.regenie.raw_p")),
  list(label = "L5.ET",
       pipeline = file.path(PIPELINE_DIR, "NABEC_L5.ET_step2.regenie.raw_p"),
       ind      = file.path(IND_DIR,      "NABEC_WGS_step2_update_hg19_L5.ET.regenie.raw_p"))
)

read_gwas <- function(path, label) {
  cat(paste0("  Reading ", label, ": ", basename(path), "\n"))
  dt <- fread(path, header = TRUE, showProgress = FALSE) %>%
    as.data.frame() %>%
    rename(POS = GENPOS) %>%
    mutate(
      CHROM = as.numeric(gsub("chr", "", CHROM)),
      POS   = as.numeric(POS),
      P     = as.numeric(P)
    ) %>%
    filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1) %>%
    select(ID, CHROM, POS, P)

  n_total  <- nrow(dt)
  lambda   <- median(qchisq(1 - dt$P, 1), na.rm = TRUE) / qchisq(0.5, 1)
  n_sig    <- sum(dt$P < 5e-8, na.rm = TRUE)
  n_sugg   <- sum(dt$P < SUGGESTIVE, na.rm = TRUE)
  cat(paste0("    Variants: ", n_total,
             "  |  Lambda: ", round(lambda, 3),
             "  |  GW-sig: ", n_sig,
             "  |  Suggestive: ", n_sugg, "\n"))

  # Downsample below suggestive threshold for faster plotting
  important <- dt %>% filter(P < SUGGESTIVE)
  other     <- dt %>% filter(P >= SUGGESTIVE)
  if (nrow(other) > DOWNSAMPLE_N) {
    set.seed(42)
    other <- other %>% sample_n(DOWNSAMPLE_N)
  }
  plot_dt <- bind_rows(important, other) %>% arrange(CHROM, POS)

  list(data = plot_dt, lambda = lambda, n = n_total, n_sig = n_sig, n_sugg = n_sugg)
}

make_manhattan <- function(gwas, title_str) {
  tryCatch(
    manhattanExtra(
      df                 = gwas$data,
      genome_wide_thresh = 5e-8,
      suggestive_thresh  = SUGGESTIVE,
      annotate           = SUGGESTIVE,
      build              = 37,
      label_size         = 3.5,
      title = paste0(title_str,
                     "  [n=", gwas$n,
                     ", λ=", round(gwas$lambda, 3),
                     ", GW-sig=", gwas$n_sig,
                     ", Sugg=", gwas$n_sugg, "]")
    ),
    error = function(e) {
      cat(paste0("  Warning: topr failed (", e$message, "), using basic ggplot\n"))
      ggplot(gwas$data, aes(x = POS, y = -log10(P))) +
        geom_point(alpha = 0.4, size = 0.3) +
        geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
        geom_hline(yintercept = -log10(1e-5),  color = "blue", linetype = "dashed") +
        facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_blank(), strip.text = element_text(size = 7)) +
        labs(title = title_str, x = "Chromosome", y = expression(-log[10](p)))
    }
  )
}

for (ct in cell_types) {
  label     <- ct$label
  label_safe <- gsub("\\.", "_", label)
  out_file  <- file.path(OUT_DIR, paste0("NABEC_", label_safe, "_pipeline_vs_ind.png"))

  cat(paste0("\n=== ", label, " ===\n"))

  pipeline_gwas <- read_gwas(ct$pipeline, "Pipeline")
  ind_gwas      <- read_gwas(ct$ind,      "Ind. Scripts")

  p_pipeline <- make_manhattan(pipeline_gwas,
                               paste0(label, " — NABEC Pipeline (", pipeline_gwas$n, " variants)"))
  p_ind      <- make_manhattan(ind_gwas,
                               paste0(label, " — NABEC Ind. Scripts (", ind_gwas$n, " variants)"))

  title_grob <- ggdraw() +
    draw_label(
      paste0("NABEC ", label, ": Pipeline vs Individual Scripts\n",
             "Pipeline: 210 samples (updated pipeline, Mar 2026)  |  ",
             "Ind. Scripts: 214 samples (individual analysis)"),
      fontface = "bold", size = 12, hjust = 0.5
    )

  combined <- plot_grid(title_grob, p_pipeline, p_ind,
                        ncol = 1, rel_heights = c(0.08, 1, 1))

  ggsave(out_file, plot = combined, width = 16, height = 11, dpi = 200)
  cat(paste0("  Saved: ", out_file, "\n"))
}

cat("\nAll plots done!\n")
