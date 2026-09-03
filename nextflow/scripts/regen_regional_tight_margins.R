#!/usr/bin/env Rscript
# Regenerate the 3 focus regional coloc plots with tight margins for Figure 4.
# Removes the title banner panel and zero-pads all outer margins so there is
# no whitespace gap above/below the content when embedded in a cowplot grid.
#
# Targets:
#   VIP × MDD MDD2025          (chr7:12,284,430)
#   L5.6.IT.Car3 × BD bip2024  (chr12:2,324,042)
#   Microglia × BD bip2024     (chr6:164,862,615)
#
# Output: results/coloc/full/plots/regional/  (overwrites the originals)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(topr)
})

ROOT     <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
LOCI_DIR <- file.path(ROOT, "results/coloc/loci_full")
DIS_DIR  <- file.path(ROOT, "results/coloc/disease_gwas")
OUT_DIR  <- file.path(ROOT, "results/coloc/full/plots/regional")

# ── helper functions (adapted from scripts/coloc/05_regional_coloc_plots.R) ──

is_strand_ambiguous <- function(ref, alt) {
  paste(toupper(ref), toupper(alt)) %in% c("A T", "T A", "C G", "G C")
}

harmonize_snps <- function(d1, d2) {
  m <- merge(d1, d2, by = "snp", suffixes = c("_ct", "_dis"))
  if (nrow(m) == 0) return(character(0))
  m <- m[!is_strand_ambiguous(m$ref_ct, m$alt_ct), , drop = FALSE]
  if (nrow(m) == 0) return(character(0))
  flipped <- toupper(m$ref_ct) == toupper(m$alt_dis) &
             toupper(m$alt_ct) == toupper(m$ref_dis)
  matched <- (toupper(m$ref_ct) == toupper(m$ref_dis) &
              toupper(m$alt_ct) == toupper(m$alt_dis)) | flipped
  m$snp[matched]
}

neglog10p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  p[p <= 0] <- .Machine$double.xmin
  -log10(p)
}

# ── tight-margin GWAS panel ───────────────────────────────────────────────────
# show_legend only on the first (CTP) panel to avoid duplication
make_gwas_panel <- function(plot_data, panel_label, lead_pos_mb,
                            lead_points, show_legend, show_x_axis) {
  p <- ggplot(plot_data, aes(x = pos / 1e6, y = neglog10p)) +
    geom_vline(xintercept = lead_pos_mb, linetype = "dashed",
               linewidth = 0.3, color = "grey50") +
    geom_point(aes(color = in_coloc), alpha = 0.65, size = 1.2) +
    scale_color_manual(
      values = c("FALSE" = "grey55", "TRUE" = "#c0392b"),
      labels = c("Other SNPs", "Harmonized SNPs"),
      name   = NULL
    ) +
    labs(
      title = panel_label,
      x     = if (show_x_axis) sprintf("Chromosome %s position (Mb)",
                                       unique(plot_data$chr_label)) else NULL,
      y     = expression(-log[10](p))
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      legend.position  = if (show_legend) "right" else "none",
      axis.title.x     = if (show_x_axis) element_text() else element_blank(),
      axis.text.x      = if (show_x_axis) element_text() else element_blank(),
      # ── tight margins: zero top/bottom so cowplot grid has no gaps ──
      plot.margin      = margin(t = 2, r = 4, b = 0, l = 4)
    )
  if (nrow(lead_points) >= 1)
    p <- p + geom_point(data = lead_points,
                        aes(x = pos / 1e6, y = neglog10p),
                        inherit.aes = FALSE,
                        color = "#1f78b4", size = 2.8, shape = 18)
  p
}

# ── gene track via topr::regionplot ──────────────────────────────────────────
prepare_topr_df <- function(ct_plot, lead_snp) {
  df <- ct_plot %>%
    mutate(CHROM = as.integer(chr), POS = as.integer(pos),
           P = as.numeric(p), ID = snp,
           Effect = if ("beta" %in% names(.)) as.numeric(beta) else 0) %>%
    filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1)
  if (nrow(df) == 0) return(NULL)
  df_pos <- df %>% filter(is.na(Effect) | Effect >= 0) %>%
    select(CHROM, POS, P, ID, Effect)
  df_neg <- df %>% filter(!is.na(Effect), Effect < 0)  %>%
    select(CHROM, POS, P, ID, Effect)
  if (nrow(df_neg) == 0) return(list(df_pos))
  if (lead_snp %in% df_pos$ID) list(df_pos, df_neg) else list(df_neg, df_pos)
}

make_gene_track <- function(ct_plot, lead_snp, region_size, genome_build = 37) {
  datasets <- prepare_topr_df(ct_plot, lead_snp)
  if (is.null(datasets)) return(NULL)
  plots <- tryCatch(
    regionplot(datasets, variant = lead_snp, region_size = region_size,
               build = genome_build, show_overview = FALSE,
               extract_plots = TRUE, title = NULL, legend_name = NULL),
    error = function(e) { warning("Gene track error: ", e$message); NULL }
  )
  if (is.null(plots)) return(NULL)
  plots$gene_plot +
    theme(plot.margin = margin(t = 0, r = 4, b = 4, l = 4),
          axis.title.x = element_text(size = 10))
}

# ── main plot function ────────────────────────────────────────────────────────
plot_regional_tight <- function(cell_type, locus_id, disease,
                                pp_h4, window_kb = 500) {
  cat(sprintf("\nPlotting %s x %s ...\n", locus_id, disease))

  # --- load locus metadata ---
  loci_file <- file.path(LOCI_DIR, paste0(cell_type, "_loci.tsv"))
  loci      <- fread(loci_file, sep = "\t", data.table = FALSE)
  locus     <- loci[loci$locus_id == locus_id, , drop = FALSE][1, ]

  win_start  <- max(1L, as.integer(locus$lead_pos) - window_kb * 1000L)
  win_end    <- as.integer(locus$lead_pos) + window_kb * 1000L
  locus_chr  <- as.character(locus$chr)
  lead_pos_mb <- locus$lead_pos / 1e6
  region_size <- win_end - win_start

  # --- load cell-type GWAS (from pre-computed locus_data.tsv.gz) ---
  ct_all <- fread(file.path(LOCI_DIR, paste0(cell_type, "_locus_data.tsv.gz")),
                  sep = "\t", data.table = FALSE) %>%
    mutate(chr = as.character(chr), pos = as.integer(pos),
           p = as.numeric(p), ref = toupper(ref), alt = toupper(alt))

  ct_plot <- ct_all %>%
    filter(pos >= win_start, pos <= win_end) %>%
    mutate(neglog10p = neglog10p(p), chr_label = locus_chr)

  # --- load disease GWAS ---
  dis_file <- file.path(DIS_DIR, disease, paste0(disease, "_hg19.tsv"))
  dis_all  <- fread(dis_file, sep = "\t",
                    select = c("snp", "chr", "pos", "ref", "alt", "p"),
                    showProgress = FALSE)
  dis_all[, chr := as.character(chr)][, pos := as.integer(pos)]
  dis_plot <- dis_all[chr == locus_chr & pos >= win_start & pos <= win_end] %>%
    as.data.frame() %>%
    mutate(p = as.numeric(p), neglog10p = neglog10p(p),
           chr_label = locus_chr)

  if (nrow(dis_plot) == 0)
    stop("No disease SNPs in window for ", locus_id, " x ", disease)

  # --- harmonize ---
  shared_snps <- harmonize_snps(
    ct_all  %>% select(snp, ref, alt),
    dis_plot %>% select(snp, ref, alt)
  )
  ct_plot  <- ct_plot  %>% mutate(in_coloc = snp %in% shared_snps,
                                   lead     = snp == locus$lead_snp)
  dis_plot <- dis_plot %>% mutate(in_coloc = snp %in% shared_snps,
                                   lead     = pos == as.integer(locus$lead_pos))

  lead_ct  <- ct_plot  %>% filter(lead) %>% slice(1)
  lead_dis <- dis_plot %>% arrange(p)   %>% slice(1)

  disease_label <- sub("_", " ", disease)

  # --- build panels ---
  ct_panel  <- make_gwas_panel(ct_plot,
                               sprintf("%s cell-type GWAS", cell_type),
                               lead_pos_mb, lead_ct,
                               show_legend = TRUE, show_x_axis = FALSE)
  dis_panel <- make_gwas_panel(dis_plot,
                               sprintf("%s disease GWAS", disease_label),
                               lead_pos_mb, lead_dis,
                               show_legend = FALSE, show_x_axis = FALSE)

  gene_panel <- make_gene_track(ct_plot, locus$lead_snp, region_size)

  # --- subtitle line (PP.H4 + locus info) embedded in CTP panel title ---
  # (replaces the separate title banner, which caused the top gap)
  subtitle_txt <- sprintf(
    "%s  |  lead %s:%s  |  PP.H4 = %.3f  |  %d harmonized SNPs",
    locus_id, locus_chr, format(locus$lead_pos, big.mark = ","),
    pp_h4, length(shared_snps)
  )
  ct_panel <- ct_panel +
    labs(subtitle = subtitle_txt) +
    theme(plot.subtitle = element_text(size = 8, colour = "grey40"))

  # --- assemble body (NO separate title banner → eliminates top gap) ---
  if (!is.null(gene_panel)) {
    dis_panel <- dis_panel +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    body <- plot_grid(ct_panel, dis_panel, gene_panel,
                      ncol = 1, rel_heights = c(3, 3, 1.6),
                      align = "v", axis = "lr")
  } else {
    dis_panel <- dis_panel +
      labs(x = sprintf("Chromosome %s position (Mb)", locus_chr))
    body <- plot_grid(ct_panel, dis_panel,
                      ncol = 1, rel_heights = c(1, 1),
                      align = "v", axis = "lr")
  }

  # Tight outer margin on the assembled body
  body <- body + theme(plot.margin = margin(0, 0, 0, 0))

  # --- save ---
  out_file <- file.path(
    OUT_DIR,
    sprintf("%s_%s_%s_regional.png",
            cell_type,
            sub("^[^_]+_", "", locus_id),
            disease)
  )
  ggsave(out_file, body, width = 9, height = 10, dpi = 200, bg = "white")
  cat(sprintf("  Saved: %s\n", basename(out_file)))
  invisible(out_file)
}

# ── regenerate the 3 focus plots ─────────────────────────────────────────────
hits <- list(
  list(cell_type = "VIP",
       locus_id  = "VIP_chr7_12284430",
       disease   = "MDD_MDD2025",
       pp_h4     = 1.000),
  list(cell_type = "L5.6.IT.Car3",
       locus_id  = "L5.6.IT.Car3_chr12_2324042",
       disease   = "BD_bip2024",
       pp_h4     = 0.877),
  list(cell_type = "Microglia",
       locus_id  = "Microglia_chr6_164862615",
       disease   = "BD_bip2024",
       pp_h4     = 0.616)
)

for (h in hits) {
  tryCatch(
    do.call(plot_regional_tight, h),
    error = function(e) warning("Failed: ", h$locus_id, " — ", e$message)
  )
}

cat("\nDone. Regional plots regenerated with tight margins.\n")
