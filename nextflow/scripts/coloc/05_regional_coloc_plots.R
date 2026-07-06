#!/usr/bin/env Rscript
# Regional locus plots for significant colocalization hits.
# Two-panel -log10(p) plots: cell-type GWAS (top) vs disease GWAS (bottom).
#
# Usage:
#   Rscript scripts/coloc/05_regional_coloc_plots.R \
#     --sig_file   results/coloc/coloc_significant.tsv \
#     --loci_dir   results/coloc/loci \
#     --disease_dir results/coloc/disease_gwas \
#     --output_dir results/coloc/plots/regional

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(topr)
})

option_list <- list(
  make_option("--sig_file",    type = "character", default = "results/coloc/coloc_significant.tsv"),
  make_option("--loci_dir",    type = "character", default = "results/coloc/loci"),
  make_option("--disease_dir", type = "character", default = "results/coloc/disease_gwas"),
  make_option("--meta_dir",    type = "character", default = NULL,
              help = "15-cohort meta-analysis dir; used when locus_data.tsv.gz is missing"),
  make_option("--output_dir",  type = "character", default = "results/coloc/plots/regional"),
  make_option("--pp_h4_sig",   type = "double",    default = 0.5),
  make_option("--method",      type = "character", default = "coloc.abf"),
  make_option("--window_kb",   type = "integer",   default = 500,
              help = "Plot window +/- kb around lead SNP [default: 500]"),
  make_option("--genome_build", type = "integer", default = 37,
              help = "Genome build for topr gene track (37 = hg19) [default: 37]"),
  make_option("--skip_genes",  action = "store_true", default = FALSE,
              help = "Skip topr gene annotation track")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

is_strand_ambiguous <- function(ref, alt) {
  paste(toupper(ref), toupper(alt)) %in% c("A T", "T A", "C G", "G C")
}

harmonize_snps <- function(d1, d2) {
  m <- merge(d1, d2, by = "snp", suffixes = c("_ct", "_dis"))
  if (nrow(m) == 0) return(character(0))

  strand_amb <- is_strand_ambiguous(m$ref_ct, m$alt_ct)
  m <- m[!strand_amb, , drop = FALSE]
  if (nrow(m) == 0) return(character(0))

  flipped <- (toupper(m$ref_ct) == toupper(m$alt_dis)) &
             (toupper(m$alt_ct) == toupper(m$ref_dis))
  matched <- (toupper(m$ref_ct) == toupper(m$ref_dis) &
              toupper(m$alt_ct) == toupper(m$alt_dis)) | flipped
  m$snp[matched]
}

neglog10p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  p[p <= 0] <- .Machine$double.xmin
  -log10(p)
}

load_sig_hits <- function(sig_file, pp_h4_sig, method) {
  if (!file.exists(sig_file)) {
    stop("Significant results file not found: ", sig_file)
  }
  hits <- fread(sig_file, sep = "\t", data.table = FALSE)
  if (nrow(hits) == 0) {
    stop("No rows in ", sig_file)
  }
  hits %>%
    mutate(PP.H4 = as.numeric(PP.H4)) %>%
    filter(PP.H4 >= pp_h4_sig, coloc_method == method) %>%
    arrange(desc(PP.H4))
}

find_disease_file <- function(disease_dir, disease_label) {
  matches <- list.files(
    disease_dir,
    pattern = paste0("^", disease_label, "_hg19\\.tsv$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(matches) == 0) {
    stop("Disease GWAS not found for ", disease_label, " under ", disease_dir)
  }
  matches[1]
}

find_meta_tbl <- function(meta_dir, cell_type) {
  matches <- list.files(
    meta_dir,
    pattern = paste0("^", gsub("([.])", "\\\\.", cell_type), "_meta_analysis_.*\\.tbl$"),
    full.names = TRUE
  )
  matches <- matches[!grepl("[0-9]\\.tbl$", basename(matches))]
  if (length(matches) == 0) {
    stop("Meta-analysis .tbl not found for ", cell_type, " under ", meta_dir)
  }
  matches[1]
}

read_meta_window <- function(tbl_file, locus_chr, win_start, win_end) {
  cmd <- sprintf(
    paste0(
      "awk -F'\\t' 'NR==1 {print; next} ",
      "{ split($1,a,\":\"); c=a[1]; gsub(/^chr/,\"\",c); pos=a[2]+0; ",
      "if(c==\"%s\" && pos>=%d && pos<=%d) print }' %s"
    ),
    locus_chr, win_start, win_end, shQuote(tbl_file)
  )
  dt <- fread(cmd = cmd, sep = "\t", showProgress = FALSE)
  if (nrow(dt) == 0) return(dt)
  setnames(dt, c("MarkerName", "Effect", "StdErr", "P-value"),
           c("snp", "beta", "se", "p"), skip_absent = TRUE)
  dt[, p := as.numeric(p)]
  dt[, beta := as.numeric(beta)]
  dt <- dt[!is.na(p) & p > 0 & p <= 1]
  dt[, c("chr", "pos", "ref", "alt") := tstrsplit(snp, ":", fixed = TRUE)]
  dt[, chr := sub("^chr", "", chr, ignore.case = TRUE)]
  dt[, pos := as.integer(pos)]
  dt[, ref := toupper(ref)]
  dt[, alt := toupper(alt)]
  dt[!is.na(pos)]
}

load_celltype_window <- function(cell_type, loci_dir, meta_dir, locus_chr, win_start, win_end) {
  locus_data_file <- file.path(loci_dir, paste0(cell_type, "_locus_data.tsv.gz"))
  if (file.exists(locus_data_file)) {
    ct_all <- fread(locus_data_file, sep = "\t", data.table = FALSE)
    return(ct_all %>%
      mutate(
        chr = as.character(chr),
        pos = as.integer(pos),
        p = as.numeric(p),
        ref = toupper(ref),
        alt = toupper(alt)
      ))
  }
  if (is.null(meta_dir) || !dir.exists(meta_dir)) {
    stop("Missing locus data for ", cell_type,
         " and --meta_dir not set (or directory not found)")
  }
  tbl_file <- find_meta_tbl(meta_dir, cell_type)
  cat(sprintf("  Loading cell-type GWAS from %s\n", basename(tbl_file)))
  read_meta_window(tbl_file, locus_chr, win_start, win_end) %>%
    as.data.frame()
}

load_disease_window <- function(dis_file, disease_label, locus_chr, win_start, win_end, cache) {
  if (is.null(cache[[dis_file]])) {
    cat(sprintf("  Loading disease file: %s\n", basename(dis_file)))
    dis_all <- fread(
      dis_file,
      sep = "\t",
      select = c("snp", "chr", "pos", "ref", "alt", "p"),
      showProgress = FALSE
    )
    dis_all[, chr := as.character(chr)]
    dis_all[, pos := as.integer(pos)]
    cache[[dis_file]] <- dis_all
  }

  dis_all <- cache[[dis_file]]
  dis_all[chr == locus_chr & pos >= win_start & pos <= win_end] %>%
    as.data.frame() %>%
    mutate(
      p = as.numeric(p),
      neglog10p = neglog10p(p),
      panel = paste0(disease_label, " GWAS")
    )
}

make_gwas_panel <- function(plot_data, panel_label, lead_pos_mb, lead_points,
                            show_legend = FALSE, show_x_axis = FALSE) {
  p <- ggplot(plot_data, aes(x = pos / 1e6, y = neglog10p)) +
    geom_vline(xintercept = lead_pos_mb, linetype = "dashed",
               linewidth = 0.3, color = "grey50") +
    geom_point(aes(color = in_coloc), alpha = 0.65, size = 1.2) +
    scale_color_manual(
      values = c("FALSE" = "grey55", "TRUE" = "#c0392b"),
      labels = c("Other SNPs", "Harmonized SNPs"),
      name = NULL
    ) +
    labs(
      title = panel_label,
      x = if (show_x_axis) {
        sprintf("Chromosome %s position (Mb)", unique(plot_data$chr_label))
      } else {
        NULL
      },
      y = expression(-log[10](p))
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = if (show_legend) "bottom" else "none",
      axis.title.x = if (show_x_axis) element_text() else element_blank(),
      axis.text.x = if (show_x_axis) element_text() else element_blank()
    )

  if (nrow(lead_points) >= 1) {
    p <- p + geom_point(
      data = lead_points,
      aes(x = pos / 1e6, y = neglog10p),
      inherit.aes = FALSE,
      color = "#1f78b4", size = 2.8, shape = 18
    )
  }
  p
}

prepare_topr_datasets <- function(ct_plot, lead_snp) {
  df <- ct_plot %>%
    mutate(
      CHROM = as.integer(chr),
      POS = as.integer(pos),
      P = as.numeric(p),
      ID = snp,
      Effect = if ("beta" %in% names(.)) as.numeric(beta) else 0
    ) %>%
    filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1)

  if (nrow(df) == 0) return(NULL)

  df_pos <- df %>% filter(is.na(Effect) | Effect >= 0) %>% select(CHROM, POS, P, ID, Effect)
  df_neg <- df %>% filter(!is.na(Effect), Effect < 0) %>% select(CHROM, POS, P, ID, Effect)

  if (nrow(df_neg) == 0) {
    return(list(df_pos))
  }
  if (lead_snp %in% df_pos$ID) {
    list(df_pos, df_neg)
  } else {
    list(df_neg, df_pos)
  }
}

make_gene_track <- function(ct_plot, lead_snp, region_size, genome_build) {
  datasets <- prepare_topr_datasets(ct_plot, lead_snp)
  if (is.null(datasets)) return(NULL)

  plots <- regionplot(
    datasets,
    variant = lead_snp,
    region_size = region_size,
    build = genome_build,
    show_overview = FALSE,
    extract_plots = TRUE,
    title = NULL,
    legend_name = NULL
  )
  plots$gene_plot +
    theme(
      plot.margin = margin(0, 5.5, 5.5, 5.5),
      axis.title.x = element_text(size = 10)
    )
}

plot_regional_hit <- function(hit, loci_dir, disease_dir, meta_dir, output_dir,
                              window_kb, genome_build, skip_genes, cache) {
  cell_type <- hit$cell_type
  locus_id  <- hit$locus_id
  disease   <- hit$disease
  pp_h4     <- hit$PP.H4

  loci_file <- file.path(loci_dir, paste0(cell_type, "_loci.tsv"))
  if (!file.exists(loci_file)) {
    warning("Missing loci file for ", cell_type, " — skipping ", locus_id)
    return(invisible(NULL))
  }

  loci <- fread(loci_file, sep = "\t", data.table = FALSE)
  locus <- loci[loci$locus_id == locus_id, , drop = FALSE]
  if (nrow(locus) == 0) {
    warning("Locus not found in ", loci_file, ": ", locus_id)
    return(invisible(NULL))
  }
  locus <- locus[1, ]

  win_start <- max(1L, as.integer(locus$lead_pos) - window_kb * 1000L)
  win_end   <- as.integer(locus$lead_pos) + window_kb * 1000L
  locus_chr <- as.character(locus$chr)

  ct_all <- load_celltype_window(cell_type, loci_dir, meta_dir, locus_chr, win_start, win_end) %>%
    mutate(
      neglog10p = neglog10p(p),
      panel = paste0(cell_type, " GWAS")
    )

  ct_plot <- ct_all %>%
    filter(pos >= win_start, pos <= win_end) %>%
    mutate(chr_label = locus_chr)

  dis_file <- find_disease_file(disease_dir, disease)
  dis_plot <- load_disease_window(dis_file, disease, locus_chr, win_start, win_end, cache)
  if (nrow(dis_plot) == 0) {
    warning("No disease SNPs in window for ", locus_id, " x ", disease)
    return(invisible(NULL))
  }

  shared_snps <- harmonize_snps(
    ct_all %>% select(snp, ref, alt),
    dis_plot %>% select(snp, ref, alt)
  )

  ct_plot <- ct_plot %>%
    mutate(
      in_coloc = snp %in% shared_snps,
      lead = snp == locus$lead_snp
    )
  dis_plot <- dis_plot %>%
    mutate(
      in_coloc = snp %in% shared_snps,
      lead = pos == as.integer(locus$lead_pos),
      chr_label = locus_chr
    )

  lead_ct <- ct_plot %>% filter(lead) %>% slice(1)
  lead_dis <- dis_plot %>% arrange(p) %>% slice(1)
  lead_pos_mb <- locus$lead_pos / 1e6
  region_size <- (win_end - win_start)

  disease_label <- sub("_", " ", disease)
  ct_panel <- make_gwas_panel(
    ct_plot,
    panel_label = sprintf("%s cell-type GWAS", cell_type),
    lead_pos_mb = lead_pos_mb,
    lead_points = lead_ct,
    show_legend = TRUE,
    show_x_axis = FALSE
  )
  dis_panel <- make_gwas_panel(
    dis_plot,
    panel_label = sprintf("%s disease GWAS", disease_label),
    lead_pos_mb = lead_pos_mb,
    lead_points = lead_dis,
    show_legend = FALSE,
    show_x_axis = FALSE
  )

  title <- ggdraw() +
    draw_label(
      sprintf("%s × %s colocalization", cell_type, disease_label),
      fontface = "bold", size = 14, x = 0.01, hjust = 0
    ) +
    draw_label(
      sprintf(
        "%s | lead %s:%s | PP.H4 = %.3f | %d harmonized SNPs",
        locus_id, locus_chr, format(locus$lead_pos, big.mark = ","),
        pp_h4, length(shared_snps)
      ),
      fontface = "plain", size = 10, x = 0.01, y = 0.2, hjust = 0
    )

  if (!skip_genes) {
    gene_panel <- tryCatch(
      make_gene_track(ct_plot, locus$lead_snp, region_size, genome_build),
      error = function(e) {
        warning("Gene track failed for ", locus_id, ": ", conditionMessage(e))
        NULL
      }
    )
  } else {
    gene_panel <- NULL
  }

  if (!is.null(gene_panel)) {
    dis_panel <- dis_panel +
      labs(x = NULL) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    body <- plot_grid(
      ct_panel, dis_panel, gene_panel,
      ncol = 1,
      rel_heights = c(3, 3, 1.6),
      align = "v",
      axis = "lr"
    )
  } else {
    dis_panel <- dis_panel +
      labs(x = sprintf("Chromosome %s position (Mb)", locus_chr))
    body <- plot_grid(
      ct_panel, dis_panel,
      ncol = 1,
      rel_heights = c(1, 1),
      align = "v",
      axis = "lr"
    )
  }

  combined <- plot_grid(title, body, ncol = 1, rel_heights = c(0.08, 1))

  out_file <- file.path(
    output_dir,
    sprintf("%s_%s_%s_regional.png", cell_type, sub("^[^_]+_", "", locus_id), disease)
  )
  ggsave(out_file, combined, width = 10, height = 8, dpi = 200)
  cat(sprintf("  Written %s\n", out_file))
  invisible(out_file)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
hits <- load_sig_hits(opt$sig_file, opt$pp_h4_sig, opt$method)
if (nrow(hits) == 0) {
  stop(sprintf("No hits with PP.H4 >= %.2f and method %s", opt$pp_h4_sig, opt$method))
}

cat(sprintf("Generating regional plots for %d significant hits...\n", nrow(hits)))
dis_cache <- list()

for (i in seq_len(nrow(hits))) {
  hit <- hits[i, ]
  cat(sprintf("[%d/%d] %s × %s\n", i, nrow(hits), hit$locus_id, hit$disease))
  tryCatch(
    plot_regional_hit(
      hit = hit,
      loci_dir = opt$loci_dir,
      disease_dir = opt$disease_dir,
      meta_dir = opt$meta_dir,
      output_dir = opt$output_dir,
      window_kb = opt$window_kb,
      genome_build = opt$genome_build,
      skip_genes = opt$skip_genes,
      cache = dis_cache
    ),
    error = function(e) {
      warning("Failed ", hit$locus_id, " × ", hit$disease, ": ", conditionMessage(e))
    }
  )
}

cat("\nDone.\n")
