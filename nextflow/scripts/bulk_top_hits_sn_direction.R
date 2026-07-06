#!/usr/bin/env Rscript
# Identify ALL suggestive-threshold bulk meta-analysis hits (p < sugg_thresh)
# across all 19 cell types, prune to independent loci (greedy distance-based
# clumping within each chromosome), look up each lead SNP in the matched
# snRNAseq meta-analysis, and report effect-direction concordance.
#
# Key outputs (in <outdir>/):
#   bulk_suggestive_hits_sn_direction.tsv   – all suggestive SNPs + sn lookup
#   independent_loci_concordant.tsv         – pruned, concordant loci only
#                                             (ready for colocalization / LocusZoom)
#   bulk_top_hits_direction_scatter.png
#   bulk_top_hits_direction_bar.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

# ---- options ----------------------------------------------------------------
option_list <- list(
  make_option("--bulk_dir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_13cohorts_design_matrix"),
  make_option("--sn_dir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_sn_rosmap_msbb_hbcc_alias15"),
  make_option("--outdir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity_design_matrix/top_hits"),
  make_option("--sugg_thresh", type = "double", default = 1e-5),
  make_option("--clump_window_kb", type = "integer", default = 500,
    help = "Distance window (kb) for greedy LD pruning [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

SUGG        <- opt$sugg_thresh
CLUMP_WIN   <- opt$clump_window_kb * 1000L   # convert to bp

# ---- cell-type mapping ------------------------------------------------------
# All 19 bulk cell types → nearest matched sn cell type.
crosswalk <- data.frame(
  bulk  = c("Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3",
            "L5.6.NP", "L5.ET", "L6b", "L6.CT",
            "LAMP5", "Microglia", "Oligodendrocyte", "OPC",
            "PAX6", "Pericyte", "PVALB", "SST", "VIP", "VLMC"),
  sn    = c("ASTRO", "ENDO", "L23IT", "L23IT", "L23IT",
            "L56NP", "L5ET", "L6B", "L6CT",
            "LAMP5LHX6", "PVM", "OLIGO", "OPC",
            "L23IT", "ENDO", "PVALB", "SST", "VIP", "VLMC"),
  # flag the 4 bulk types that have no direct sn counterpart
  direct_match = c(TRUE, TRUE, TRUE, FALSE, FALSE,
                   TRUE, TRUE, TRUE, TRUE,
                   TRUE, FALSE, TRUE, TRUE,
                   FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# ---- helper: locate .tbl file -----------------------------------------------
find_tbl <- function(dir, cell_type) {
  pat <- paste0("^", gsub("([.])", "\\\\.", cell_type),
                "_meta_analysis_.*\\.tbl$")
  f <- list.files(dir, pattern = pat, full.names = TRUE)
  f <- f[!grepl("1\\.tbl$", basename(f))]
  if (length(f) == 0) return(NA_character_)
  f[1]
}

# ---- helper: read .tbl, keep SNP-only rows, add Z-score --------------------
read_meta <- function(path, label) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  message("  reading ", basename(path))
  dt <- fread(path, sep = "\t", data.table = TRUE, verbose = FALSE)
  setnames(dt, tolower(gsub("[^A-Za-z0-9]", "_", names(dt))))
  col_map <- list(
    marker  = "markername",
    allele1 = "allele1",
    allele2 = "allele2",
    effect  = "effect",
    stderr  = "stderr",
    pvalue  = "p_value"
  )
  for (new_nm in names(col_map)) {
    old_nm <- col_map[[new_nm]]
    if (old_nm %in% names(dt) && !(new_nm %in% names(dt)))
      setnames(dt, old_nm, new_nm)
  }
  needed <- c("marker", "allele1", "allele2", "effect", "stderr", "pvalue")
  if (!all(needed %in% names(dt))) {
    message("    WARNING: missing columns – skipping. Present: ",
            paste(names(dt), collapse = ", "))
    return(NULL)
  }
  dt <- dt[, needed, with = FALSE]
  dt <- dt[nchar(allele1) == 1 & nchar(allele2) == 1]
  parts <- tstrsplit(dt$marker, ":", fixed = TRUE)
  dt[, chrom   := parts[[1]]]
  dt[, pos     := as.integer(parts[[2]])]
  dt[, a1      := toupper(allele1)]
  dt[, a2      := toupper(allele2)]
  dt[, effect  := as.numeric(effect)]
  dt[, stderr  := as.numeric(stderr)]
  dt[, pvalue  := as.numeric(pvalue)]
  dt[, z       := effect / stderr]
  dt[, cell_type := label]
  dt[]
}

# ---- helper: greedy distance-based clumping --------------------------------
# Within each chromosome, sort by p-value ascending, then iteratively select
# the most significant SNP as a locus representative and exclude all SNPs
# within CLUMP_WIN bp on either side. Returns indices of selected lead SNPs.
clump_loci <- function(dt, window_bp) {
  dt <- copy(dt)
  dt[, row_idx := .I]
  setorder(dt, chrom, pvalue)
  selected <- integer(0)
  excluded <- logical(nrow(dt))
  for (k in seq_len(nrow(dt))) {
    if (excluded[k]) next
    lead <- dt[k]
    selected <- c(selected, lead$row_idx)
    # exclude all SNPs on the same chrom within the window
    mask <- dt$chrom == lead$chrom &
            abs(dt$pos - lead$pos) <= window_bp
    excluded[mask] <- TRUE
  }
  selected
}

# ---- helper: look up one SNP in sn_dt with allele harmonization ------------
lookup_sn <- function(bulk_row, sn_dt) {
  hit <- sn_dt[chrom == bulk_row$chrom & pos == bulk_row$pos]
  if (nrow(hit) == 0)
    return(data.table(sn_found = FALSE, sn_effect = NA_real_,
                      sn_stderr = NA_real_, sn_p = NA_real_,
                      sn_z = NA_real_, concordant = NA))
  h <- hit[1]
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  ba1 <- bulk_row$a1; ba2 <- bulk_row$a2
  sa1 <- h$a1;        sa2 <- h$a2
  if      (ba1 == sa1 && ba2 == sa2)                         flip <- FALSE
  else if (ba1 == sa2 && ba2 == sa1)                         flip <- TRUE
  else if (!is.na(comp[ba1]) && comp[ba1] == sa1 && comp[ba2] == sa2) flip <- FALSE
  else if (!is.na(comp[ba1]) && comp[ba1] == sa2 && comp[ba2] == sa1) flip <- TRUE
  else
    return(data.table(sn_found = FALSE, sn_effect = NA_real_,
                      sn_stderr = NA_real_, sn_p = NA_real_,
                      sn_z = NA_real_, concordant = NA))
  sn_eff <- if (flip) -h$effect else h$effect
  sn_z   <- if (flip) -h$z      else h$z
  data.table(sn_found   = TRUE,
             sn_effect  = sn_eff,
             sn_stderr  = h$stderr,
             sn_p       = h$pvalue,
             sn_z       = sn_z,
             concordant = (bulk_row$effect > 0) == (sn_eff > 0))
}

# ---- main loop --------------------------------------------------------------
all_results  <- vector("list", nrow(crosswalk))

for (i in seq_len(nrow(crosswalk))) {
  bulk_ct <- crosswalk$bulk[i]
  sn_ct   <- crosswalk$sn[i]
  direct  <- crosswalk$direct_match[i]
  message("\nProcessing bulk: ", bulk_ct, "  →  sn: ", sn_ct,
          if (!direct) "  [proxy mapping]" else "")

  bulk_path <- find_tbl(opt$bulk_dir, bulk_ct)
  sn_path   <- find_tbl(opt$sn_dir,   sn_ct)

  bulk_dt <- read_meta(bulk_path, bulk_ct)
  if (is.null(bulk_dt)) { message("  skipping – bulk not readable"); next }

  # all suggestive hits
  sugg_bulk <- bulk_dt[pvalue < SUGG][order(pvalue)]
  message("  suggestive hits (before clumping): ", nrow(sugg_bulk))
  if (nrow(sugg_bulk) == 0) { next }

  # greedy clump to independent loci
  lead_idx  <- clump_loci(sugg_bulk, CLUMP_WIN)
  lead_bulk <- sugg_bulk[lead_idx]
  message("  independent loci after clumping (", CLUMP_WIN/1000, " kb): ",
          nrow(lead_bulk))

  sn_dt <- read_meta(sn_path, sn_ct)

  sn_lookups <- rbindlist(lapply(seq_len(nrow(lead_bulk)), function(j) {
    lookup_sn(lead_bulk[j], sn_dt)
  }))

  result <- cbind(
    data.table(
      bulk_cell_type = bulk_ct,
      sn_cell_type   = sn_ct,
      direct_match   = direct,
      marker         = lead_bulk$marker,
      chrom          = lead_bulk$chrom,
      pos            = lead_bulk$pos,
      bulk_a1        = lead_bulk$a1,
      bulk_a2        = lead_bulk$a2,
      bulk_effect    = lead_bulk$effect,
      bulk_stderr    = lead_bulk$stderr,
      bulk_p         = lead_bulk$pvalue,
      bulk_z         = lead_bulk$z
    ),
    sn_lookups
  )

  all_results[[i]] <- result
  rm(bulk_dt, sn_dt, sugg_bulk, lead_bulk, sn_lookups, result)
  gc(verbose = FALSE)
}

results <- rbindlist(all_results, use.names = TRUE, fill = TRUE)

# ---- concordance label ------------------------------------------------------
results[, concordance_label := fcase(
  sn_found == FALSE,                     "not found in sn",
  is.na(concordant),                     "allele mismatch",
  concordant == TRUE,                    "same direction",
  concordant == FALSE,                   "opposite direction",
  default = "unknown"
)]
results[, locus_label := paste0(sub("^chr", "", chrom), ":", format(pos, big.mark = ",", scientific = FALSE))]

message("\n=== Overall direction concordance ===")
print(results[, .N, by = concordance_label])

# ---- save full table --------------------------------------------------------
out_full <- file.path(opt$outdir, "bulk_suggestive_hits_sn_direction.tsv")
fwrite(results, out_full, sep = "\t")
message("Saved full table: ", out_full)

# ---- save concordant independent loci list ----------------------------------
# These are the candidates for colocalization / LocusZoom.
# Include locus window coordinates (lead ± 500 kb) for easy subsetting.
concordant_loci <- results[concordant == TRUE][order(bulk_cell_type, bulk_p)]
concordant_loci[, locus_start := pmax(0L, pos - 500000L)]
concordant_loci[, locus_end   := pos + 500000L]
concordant_loci[, chrom_num   := as.integer(sub("^chr", "", chrom))]

# LocusZoom / coloc manifest columns
coloc_manifest <- concordant_loci[, .(
  bulk_cell_type, sn_cell_type, direct_match,
  chrom, pos, locus_start, locus_end,
  marker, bulk_a1, bulk_a2,
  bulk_effect, bulk_stderr, bulk_p, bulk_z,
  sn_effect,  sn_stderr,  sn_p,  sn_z,
  locus_label
)]
setorder(coloc_manifest, bulk_cell_type, bulk_p)

out_concordant <- file.path(opt$outdir, "independent_loci_concordant.tsv")
fwrite(coloc_manifest, out_concordant, sep = "\t")
message("Saved concordant independent loci: ", out_concordant)
message("  Total concordant independent loci: ", nrow(coloc_manifest))
message("  Across ", uniqueN(coloc_manifest$bulk_cell_type), " bulk cell types")

# ---- per-cell-type summary --------------------------------------------------
summary_dt <- results[, .(
  n_suggestive_all  = .N,
  n_found_sn        = sum(sn_found, na.rm = TRUE),
  n_concordant      = sum(concordant == TRUE,  na.rm = TRUE),
  n_discordant      = sum(concordant == FALSE, na.rm = TRUE),
  n_not_in_sn       = sum(!sn_found, na.rm = TRUE)
), by = bulk_cell_type]
summary_dt[, pct_concordant := round(100 * n_concordant / pmax(n_found_sn, 1), 1)]

out_summary <- file.path(opt$outdir, "cell_type_concordance_summary.tsv")
fwrite(summary_dt, out_summary, sep = "\t")
message("\n=== Per-cell-type summary ===")
print(summary_dt)

# ============================================================
# PLOTS
# ============================================================

# ---- 1. Scatter: bulk -log10p vs sn -log10p, coloured by concordance --------
plot_dt <- results[sn_found == TRUE]
ct_order <- summary_dt[order(-pct_concordant)]$bulk_cell_type
plot_dt[, bulk_cell_type := factor(bulk_cell_type, levels = ct_order)]
plot_dt[, neg_log10_bulk_p := -log10(bulk_p)]
plot_dt[, neg_log10_sn_p   := -log10(sn_p)]
plot_dt[, direction := ifelse(concordant, "same direction", "opposite direction")]

p1 <- ggplot(plot_dt,
       aes(x = neg_log10_bulk_p, y = neg_log10_sn_p,
           colour = direction, shape = direction)) +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", colour = "grey65", linewidth = 0.4) +
  geom_hline(yintercept = -log10(SUGG),  linetype = "dotted", colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.80, size = 2.0) +
  geom_text(data = plot_dt[neg_log10_bulk_p > -log10(SUGG) | neg_log10_sn_p > -log10(SUGG)],
            aes(label = locus_label), size = 1.9, hjust = -0.12, vjust = 0.5,
            show.legend = FALSE, check_overlap = TRUE) +
  facet_wrap(~ bulk_cell_type, ncol = 4) +
  scale_colour_manual(values = c("same direction" = "#1a7abf",
                                 "opposite direction" = "#d93b3b")) +
  scale_shape_manual(values  = c("same direction" = 16,
                                 "opposite direction" = 4)) +
  labs(
    title    = "Bulk meta-analysis independent suggestive loci — effect direction in snRNA-seq GWAS",
    subtitle = paste0("Each point = one independent locus (greedy clumping ±", opt$clump_window_kb,
                      " kb) | dashed = sn p<0.05 | dotted = sn p<", SUGG),
    x = expression(-log[10](p[bulk])),
    y = expression(-log[10](p[sn])),
    colour = "Effect direction\n(vs bulk)",
    shape  = "Effect direction\n(vs bulk)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "#e8edf2"),
        strip.text       = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())

ggsave(file.path(opt$outdir, "bulk_top_hits_direction_scatter.png"),
       p1, width = 14, height = 16, dpi = 150)
message("Saved scatter plot.")

# ---- 2. Stacked bar: concordance counts per cell type -----------------------
bar_dt <- melt(
  summary_dt[, .(bulk_cell_type,
                 `same direction`     = n_concordant,
                 `opposite direction` = n_discordant,
                 `not found in sn`    = n_not_in_sn)],
  id.vars      = "bulk_cell_type",
  variable.name = "category",
  value.name    = "count"
)
bar_dt[, bulk_cell_type := factor(bulk_cell_type,
  levels = summary_dt[order(-pct_concordant)]$bulk_cell_type)]

# label position: centre of the blue "same direction" segment
# only show label when the segment is tall enough to fit text (>= 1 unit)
summary_dt[, pct_label    := ifelse(n_concordant >= 1,
                                    paste0(pct_concordant, "%"), "")]
summary_dt[, label_y      := n_concordant / 2]   # mid-point of blue segment

p2 <- ggplot(bar_dt, aes(x = bulk_cell_type, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  # % label centred inside the blue segment, white text
  geom_text(data = summary_dt[pct_label != ""],
            aes(x = bulk_cell_type,
                y = label_y,
                label = pct_label,
                fill  = NULL),
            colour = "white", fontface = "bold", size = 3.2, vjust = 0.5) +
  # total N label just above the full bar
  geom_text(data = summary_dt,
            aes(x = bulk_cell_type,
                y = n_suggestive_all + 0.3,
                label = paste0("n=", n_suggestive_all),
                fill  = NULL),
            colour = "grey30", size = 2.8, vjust = 0) +
  scale_fill_manual(
    values = c(
      "same direction"     = "#1a7abf",
      "opposite direction" = "#d93b3b",
      "not found in sn"    = "#aaaaaa"
    ),
    # fix legend order: same → opposite → not found
    breaks = c("same direction", "opposite direction", "not found in sn")
  ) +
  labs(
    title    = "Direction concordance of bulk independent suggestive loci in snRNA-seq GWAS",
    subtitle = paste0("All bulk hits p < ", SUGG,
                      ", clumped to independent loci (±", opt$clump_window_kb,
                      " kb) | % = concordant among SNPs found in sn"),
    x    = "Bulk cell type",
    y    = "# independent loci",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(opt$outdir, "bulk_top_hits_direction_bar.png"),
       p2, width = 12, height = 6, dpi = 150)
message("Saved bar plot.")

# ---- 3. Forest-style dot plot for concordant loci only ----------------------
# One panel per cell type, x = bulk effect, y = sn effect, line shows SE.
if (nrow(concordant_loci) > 0) {
  forest_dt <- copy(concordant_loci)
  forest_dt[, bulk_cell_type := factor(bulk_cell_type,
    levels = summary_dt[order(-pct_concordant)]$bulk_cell_type)]
  forest_dt[, short_label := paste0(sub("^chr", "", chrom), ":", pos)]

  p3 <- ggplot(forest_dt,
         aes(x = bulk_effect, y = sn_effect, colour = bulk_cell_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_errorbar(aes(ymin = sn_effect - sn_stderr,
                      ymax = sn_effect + sn_stderr),
                  width = 0, alpha = 0.5) +
    geom_errorbarh(aes(xmin = bulk_effect - bulk_stderr,
                       xmax = bulk_effect + bulk_stderr),
                   height = 0, alpha = 0.5) +
    geom_point(size = 2.5) +
    geom_text(aes(label = short_label), size = 2.2, hjust = -0.15,
              show.legend = FALSE, check_overlap = TRUE) +
    labs(
      title    = "Effect sizes: bulk vs sn (concordant independent loci only)",
      subtitle = "Error bars = ±1 SE | dashed line = perfect concordance (slope 1)",
      x      = "Bulk effect size (beta)",
      y      = "sn effect size (beta)",
      colour = "Bulk cell type"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank())

  ggsave(file.path(opt$outdir, "concordant_loci_effect_scatter.png"),
         p3, width = 10, height = 8, dpi = 150)
  message("Saved concordant-loci effect scatter plot.")
}

message("\n=== Done ===")
message("Outputs in: ", opt$outdir)
message("Key file for coloc / LocusZoom: ", out_concordant)
