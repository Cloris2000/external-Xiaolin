#!/usr/bin/env Rscript
# =============================================================================
# TEST FIGURE — Lead-SNP multi-cell-type Manhattan (GWS + suggestive)
#
# Argument: among CTP GWAS loci, which are cell-type-specific vs shared,
#           and how does the suggestive landscape sit beneath GWS hits?
#
# Design:
#   - Empty chromosome Manhattan backbone (no full GWAS cloud)
#   - Suggestive leads (5e-8 ≤ p < 1e-5): small, translucent points
#   - GWS leads (p < 5e-8): larger, opaque points drawn on top
#   - Colour = cell type; shape = specific (circle) vs shared (triangle)
#   - Shared loci get small x-jitter so stacked cell types remain visible
#   - Region annotations only for GWS (keeps the plot readable)
#
# Inputs:
#   figure3_gws_loci_deduplicated.tsv
#   figure3_panelE_top_loci_matrix.tsv
#   figure3_data_checkpoint.rds  (for suggestive lead derivation; cached after)
#
# Outputs:
#   test_multict_gws_manhattan.png / .pdf / .svg
#   test_multict_suggestive_leads.tsv   (cached suggestive+GWS clumped leads)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

OUT_DIR <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"
DATA_CHECKPOINT <- file.path(OUT_DIR, "figure3_data_checkpoint.rds")
SUGG_CACHE <- file.path(OUT_DIR, "test_multict_suggestive_leads.tsv")

PVAL_GWS  <- 5e-8
PVAL_SUGG <- 1e-5
CLUMP_BP  <- 5e5L
COL_GWS   <- "#D55E00"
COL_SUGG  <- "#3A6EA5"

# Nearest-gene annotation (same helper pattern as figure3)
nearest_gene_labels <- function(chr, pos, p) {
  fallback <- paste0("chr", chr, ":", format(pos, scientific = FALSE, trim = TRUE))
  out <- tryCatch({
    ann <- topr::annotate_with_nearest_gene(
      data.frame(CHROM = as.character(chr), POS = as.integer(pos), P = as.numeric(p))
    )
    g <- if ("Gene_Symbol" %in% names(ann)) as.character(ann$Gene_Symbol) else NULL
    if (is.null(g) || length(g) != length(fallback)) fallback else g
  }, error = function(e) {
    cat("  (nearest-gene lookup failed:", conditionMessage(e), "- using chr:pos)\n")
    fallback
  })
  out[is.na(out) | out == ""] <- fallback[is.na(out) | out == ""]
  out
}

# Variant display ID: prefer rsID if present, else chr:pos (drop alleles for brevity)
variant_label <- function(lead_snp, chr, pos) {
  snp <- as.character(lead_snp)
  pos_lab <- paste0("chr", chr, ":", format(as.integer(pos), scientific = FALSE, trim = TRUE))
  ifelse(
    grepl("^rs[0-9]+$", snp, ignore.case = TRUE),
    snp,
    ifelse(grepl("^chr[0-9XYM]+:[0-9]+", snp, ignore.case = TRUE),
           sub("^((chr)?[0-9XYM]+:[0-9]+).*", "\\1", snp, ignore.case = TRUE),
           ifelse(!is.na(snp) & nzchar(snp), snp, pos_lab))
  )
}

# Greedy ±window clumping (mirrors figure3_gwas_discovery.R)
greedy_clump <- function(dt, p_col = "p", pos_col = "pos", window = CLUMP_BP) {
  keep <- logical(nrow(dt))
  used <- rep(FALSE, nrow(dt))
  pos  <- dt[[pos_col]]
  for (i in seq_len(nrow(dt))) {
    if (!used[i]) {
      keep[i] <- TRUE
      diffs <- abs(pos - pos[i])
      used[diffs <= window] <- TRUE
    }
  }
  which(keep)
}

derive_loci <- function(gwas_dt, p_thresh, window = CLUMP_BP) {
  sig <- gwas_dt[p < p_thresh]
  if (nrow(sig) == 0L) return(data.table())
  setorder(sig, p)
  result <- sig[, {
    idx <- greedy_clump(.SD, window = window)
    .SD[idx, .(lead_snp = snp, lead_pos = pos,
               lead_p = p, lead_beta = beta)]
  }, by = .(cell_type, chr)]
  result[, locus_id := paste0(cell_type, "_", chr, "_", lead_pos)]
  result
}

# Chromosome lengths (GRCh38, approx; for cumulative x)
CHR_LEN <- c(
  248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
  114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
  58617616, 64444167, 46709983, 50818468
)
names(CHR_LEN) <- as.character(1:22)

chr_offset <- c(0, cumsum(as.numeric(CHR_LEN)))[1:22]
names(chr_offset) <- as.character(1:22)

# -----------------------------------------------------------------------------
# Suggestive leads: load cache or derive from figure3 checkpoint
# -----------------------------------------------------------------------------
load_suggestive_leads <- function() {
  if (file.exists(SUGG_CACHE)) {
    cat("Loading cached suggestive leads:", SUGG_CACHE, "\n")
    return(fread(SUGG_CACHE))
  }
  if (!file.exists(DATA_CHECKPOINT)) {
    stop("Need either ", SUGG_CACHE, " or ", DATA_CHECKPOINT)
  }
  cat("Deriving suggestive leads from checkpoint (slow, one-time)...\n")
  chk <- readRDS(DATA_CHECKPOINT)
  gwas_sig <- chk$gwas_sig
  # Normalise column names used by derive_loci
  if (!"snp" %in% names(gwas_sig) && "MarkerName" %in% names(gwas_sig))
    setnames(gwas_sig, "MarkerName", "snp")
  if (!"pos" %in% names(gwas_sig) && "BP" %in% names(gwas_sig))
    setnames(gwas_sig, "BP", "pos")
  if (!"p" %in% names(gwas_sig)) {
    pcol <- intersect(c("P-value", "P", "pval", "p_value"), names(gwas_sig))
    if (length(pcol)) setnames(gwas_sig, pcol[1], "p")
  }
  loci_sugg <- derive_loci(gwas_sig, p_thresh = PVAL_SUGG)
  loci_sugg[, cell_type := as.character(cell_type)]
  loci_sugg[, chr := as.integer(chr)]
  loci_sugg[, lead_pos := as.integer(lead_pos)]
  loci_sugg[, lead_p := as.numeric(lead_p)]
  fwrite(loci_sugg, SUGG_CACHE, sep = "\t")
  cat("  Cached", nrow(loci_sugg), "suggestive leads ->", SUGG_CACHE, "\n")
  loci_sugg
}

loci_sugg_all <- load_suggestive_leads()
# Suggestive-only layer (exclude GWS; those come from the panelE/dedup path)
sugg <- loci_sugg_all[lead_p >= PVAL_GWS & lead_p < PVAL_SUGG]
sugg[, `:=`(
  source = "suggestive_clump",
  sig_tier = "Suggestive"
)]

cat("Suggestive-only leads:", nrow(sugg),
    "across", uniqueN(sugg$cell_type), "cell types\n")

# -----------------------------------------------------------------------------
# GWS layer: prefer per-CT p from panel E; backfill from dedup table
# -----------------------------------------------------------------------------
dedup <- fread(file.path(OUT_DIR, "figure3_gws_loci_deduplicated.tsv"))
panelE <- fread(file.path(OUT_DIR, "figure3_panelE_top_loci_matrix.tsv"))

dedup_exp <- dedup[, {
  cts <- trimws(unlist(strsplit(cell_types, ";")))
  .(cell_type = cts,
    chr = as.integer(chr),
    lead_pos = as.integer(lead_pos),
    lead_p = as.numeric(lead_p),
    n_cell_types = as.integer(n_cell_types),
    region_id = as.integer(region_id))
}, by = seq_len(nrow(dedup))]
dedup_exp[, seq_len := NULL]

gE <- panelE[is_gws == TRUE]
gE[, region_bin := paste0(lead_chr, "_", floor(as.numeric(lead_pos) / CLUMP_BP))]
gE <- gE[, .SD[which.min(as.numeric(p_min))], by = .(cell_type, region_bin)]
gE <- gE[, .(
  cell_type,
  chr = as.integer(lead_chr),
  lead_pos = as.integer(lead_pos),
  lead_p = as.numeric(p_min),
  region_bin,
  source = "panelE"
)]

dedup_exp[, region_bin := paste0(chr, "_", floor(lead_pos / CLUMP_BP))]
gE[, key := paste(cell_type, region_bin, sep = "|")]
dedup_exp[, key := paste(cell_type, region_bin, sep = "|")]
missing <- dedup_exp[!key %in% gE$key]
if (nrow(missing)) {
  miss_rows <- missing[, .(
    cell_type, chr, lead_pos, lead_p, region_bin, source = "dedup"
  )]
  gE <- rbind(
    gE[, .(cell_type, chr, lead_pos, lead_p, region_bin, source)],
    miss_rows,
    use.names = TRUE
  )
} else {
  gE <- gE[, .(cell_type, chr, lead_pos, lead_p, region_bin, source)]
}
gE[, sig_tier := "GWS"]

# -----------------------------------------------------------------------------
# Combine. Sharing is tier-aware so GWS labels stay GWS-only:
#   - GWS points: shared if ≥2 cell types are GWS in the 500 kb bin
#   - Suggestive: shared if ≥2 cell types (any tier) hit the same bin
# -----------------------------------------------------------------------------
sugg[, region_bin := paste0(chr, "_", floor(lead_pos / CLUMP_BP))]
pts <- rbind(
  gE[, .(cell_type, chr, lead_pos, lead_p, region_bin, source, sig_tier)],
  sugg[, .(cell_type, chr, lead_pos, lead_p, region_bin, source, sig_tier)],
  use.names = TRUE
)

n_ct_gws  <- pts[sig_tier == "GWS",
                 .(n_gws = uniqueN(cell_type)), by = region_bin]
n_ct_any  <- pts[, .(n_any = uniqueN(cell_type)), by = region_bin]
pts <- merge(pts, n_ct_gws, by = "region_bin", all.x = TRUE)
pts <- merge(pts, n_ct_any, by = "region_bin", all.x = TRUE)
pts[is.na(n_gws), n_gws := 0L]
pts[, n_cell_types := fifelse(sig_tier == "GWS", n_gws, n_any)]
pts[, specificity := fifelse(n_cell_types >= 2L, "Shared", "Specific")]
pts[, nlp := pmin(-log10(lead_p), 20)]
pts[, sig_tier := factor(sig_tier, levels = c("Suggestive", "GWS"))]

# Cumulative genomic position + within-region jitter for shared stacks
pts[, x_base := chr_offset[as.character(chr)] + lead_pos]
setorder(pts, chr, lead_pos, sig_tier, cell_type)
# Dodge within (region × tier) so suggestive cloud doesn't shove GWS sideways
pts[, idx_in_region := seq_len(.N), by = .(region_bin, sig_tier)]
pts[, n_dodge := uniqueN(cell_type), by = .(region_bin, sig_tier)]
pts[, x := x_base + (idx_in_region - (n_dodge + 1) / 2) * 2.2e6]

# -----------------------------------------------------------------------------
# GWS-only region annotations (n_cell_types = GWS count)
# -----------------------------------------------------------------------------
gws_pts <- pts[sig_tier == "GWS"]
ann <- gws_pts[, {
  i <- which.max(nlp)
  .(
    x = mean(x_base),
    nlp = max(nlp),
    n_cell_types = n_gws[1],
    specificity = fifelse(n_gws[1] >= 2L, "Shared", "Specific"),
    chr = chr[i],
    lead_pos = as.integer(lead_pos[i]),
    lead_p = lead_p[i]
  )
}, by = region_bin]

dedup[, region_bin := paste0(chr, "_", floor(as.integer(lead_pos) / CLUMP_BP))]
ann <- merge(
  ann,
  dedup[, .(region_bin, lead_snp, dedup_pos = as.integer(lead_pos))],
  by = "region_bin",
  all.x = TRUE
)
ann[!is.na(dedup_pos), lead_pos := dedup_pos]

cat("Annotating nearest genes with topr...\n")
ann[, gene := nearest_gene_labels(chr, lead_pos, lead_p)]
ann[, var_id := variant_label(lead_snp, chr, lead_pos)]
ann[, label := paste0(
  gene, "\n", var_id, "\n(",
  n_cell_types, " CT",
  fifelse(n_cell_types > 1L, "s, shared", ", specific"),
  ")"
)]

chr_mids <- data.table(
  chr = 1:22,
  mid = chr_offset + CHR_LEN / 2
)
chr_rects <- data.table(
  chr = 1:22,
  xmin = chr_offset,
  xmax = chr_offset + CHR_LEN,
  fill = rep(c("a", "b"), length.out = 22)
)

# Colour palette (Okabe–Ito extended; colourblind-safe)
ct_levels <- sort(unique(pts$cell_type))
pal <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#999999", "#332288", "#88CCEE", "#117733", "#DDCC77",
  "#AA4499", "#44AA99", "#882255", "#661100", "#6699CC", "#888888",
  "#BBBBBB"
)
ct_cols <- setNames(pal[seq_along(ct_levels)], ct_levels)

n_sugg <- nrow(pts[sig_tier == "Suggestive"])
n_gws  <- nrow(pts[sig_tier == "GWS"])
cat("Plotting", n_sugg, "suggestive +", n_gws, "GWS leads across",
    uniqueN(pts$cell_type), "cell types;",
    uniqueN(pts$region_bin), "genomic regions (",
    sum(ann$n_cell_types >= 2), "GWS-shared,",
    sum(ann$n_cell_types == 1), "GWS-specific labelled)\n")

# Split layers so suggestive sits underneath and GWS stays crisp
pts_sugg <- pts[sig_tier == "Suggestive"]
pts_gws  <- pts[sig_tier == "GWS"]

p <- ggplot() +
  geom_rect(
    data = chr_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.35, colour = NA, show.legend = FALSE
  ) +
  scale_fill_manual(values = c(a = "grey96", b = "grey90"), guide = "none") +
  # Threshold lines (suggestive under GWS)
  geom_hline(yintercept = -log10(PVAL_SUGG), linetype = "dotted",
             colour = COL_SUGG, linewidth = 0.35) +
  geom_hline(yintercept = -log10(PVAL_GWS), linetype = "dashed",
             colour = COL_GWS, linewidth = 0.45) +
  # Invisible layer: one point per cell type so the colour legend is complete
  # and opaque (real layers hide colour to avoid alpha/bleed issues)
  geom_point(
    data = data.table(x = -Inf, y = -Inf, cell_type = ct_levels),
    aes(x = x, y = y, colour = cell_type),
    size = 2.5, alpha = 1, shape = 16,
    show.legend = c(colour = TRUE)
  ) +
  # Suggestive: small + translucent
  geom_point(
    data = pts_sugg,
    aes(x = x, y = nlp, colour = cell_type, shape = specificity),
    size = 1.15, alpha = 0.38, stroke = 0.35,
    show.legend = c(colour = FALSE, shape = TRUE)
  ) +
  # GWS: larger + opaque
  geom_point(
    data = pts_gws,
    aes(x = x, y = nlp, colour = cell_type, shape = specificity),
    size = 3.2, alpha = 0.95, stroke = 0.7,
    show.legend = c(colour = FALSE, shape = TRUE)
  ) +
  scale_shape_manual(
    name = "Locus type",
    values = c(Specific = 16, Shared = 17)
  ) +
  scale_colour_manual(
    name = "Cell type",
    values = ct_cols,
    breaks = ct_levels,
    limits = ct_levels,
    drop = FALSE
  ) +
  # Dummy points for significance-tier legend only (no fixed colour —
  # that bleeds into the cell-type legend and turns keys black)
  geom_point(
    data = data.table(
      x = -Inf, y = -Inf,
      tier = factor(c("Suggestive (p < 1e-5)", "Genome-wide (p < 5e-8)"),
                    levels = c("Suggestive (p < 1e-5)", "Genome-wide (p < 5e-8)"))
    ),
    aes(x = x, y = y, size = tier, alpha = tier),
    shape = 16,
    show.legend = c(size = TRUE, alpha = TRUE, colour = FALSE, shape = FALSE)
  ) +
  scale_size_manual(
    name = "Significance",
    values = c("Suggestive (p < 1e-5)" = 1.2, "Genome-wide (p < 5e-8)" = 3.2)
  ) +
  scale_alpha_manual(
    name = "Significance",
    values = c("Suggestive (p < 1e-5)" = 0.4, "Genome-wide (p < 5e-8)" = 0.95)
  ) +
  ggrepel::geom_label_repel(
    data = ann,
    aes(x = x, y = nlp, label = label),
    size = 2.5, lineheight = 0.9,
    box.padding = 0.35, point.padding = 0.2,
    min.segment.length = 0, segment.size = 0.3,
    max.overlaps = Inf, seed = 1,
    fill = alpha("white", 0.88), label.size = 0.2,
    inherit.aes = FALSE
  ) +
  # Threshold labels (right margin)
  annotate(
    "text", x = Inf, y = -log10(PVAL_GWS),
    label = "  GWS 5e-8", hjust = 1, vjust = -0.4,
    size = 2.4, colour = COL_GWS
  ) +
  annotate(
    "text", x = Inf, y = -log10(PVAL_SUGG),
    label = "  sugg. 1e-5", hjust = 1, vjust = -0.4,
    size = 2.4, colour = COL_SUGG
  ) +
  scale_x_continuous(
    breaks = chr_mids$mid,
    labels = as.character(chr_mids$chr),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name = expression(-log[10](italic(P))),
    limits = c(0, 22),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Chromosome") +
  theme_classic(base_size = 11) +
  theme(
    plot.caption = element_text(size = 7, colour = "grey40", hjust = 0),
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    legend.key.size = unit(0.32, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7.2),
    legend.spacing.y = unit(0.08, "cm"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
  ) +
  guides(
    colour = guide_legend(
      order = 1, ncol = 1,
      override.aes = list(size = 2.5, alpha = 1, shape = 16, stroke = 0.5)
    ),
    size = guide_legend(
      order = 2,
      override.aes = list(shape = 16, colour = "grey25", alpha = 1)
    ),
    alpha = "none",
    shape = guide_legend(
      order = 3,
      override.aes = list(size = 3, colour = "grey20", alpha = 1)
    )
  )

out_png <- file.path(OUT_DIR, "test_multict_gws_manhattan.png")
out_pdf <- file.path(OUT_DIR, "test_multict_gws_manhattan.pdf")
out_svg <- file.path(OUT_DIR, "test_multict_gws_manhattan.svg")
# Slightly taller to give suggestive cloud room without crowding labels
ggsave(out_png, p, width = 11.5, height = 6.2, dpi = 200, bg = "white")
ggsave(out_pdf, p, width = 11.5, height = 6.2, bg = "white")
ggsave(out_svg, p, width = 11.5, height = 6.2, bg = "white")

fwrite(pts[order(sig_tier, chr, lead_pos, cell_type),
           .(sig_tier, cell_type, chr, lead_pos, lead_p, nlp,
             n_cell_types, specificity, source)],
       file.path(OUT_DIR, "test_multict_gws_manhattan_points.tsv"), sep = "\t")
fwrite(ann[order(chr, lead_pos),
           .(region_bin, chr, lead_pos, lead_p, gene, var_id, lead_snp,
             n_cell_types, specificity, label)],
       file.path(OUT_DIR, "test_multict_gws_manhattan_labels.tsv"), sep = "\t")

cat("Wrote:\n  ", out_png, "\n  ", out_pdf, "\n  ", out_svg, "\n",
    "  ", file.path(OUT_DIR, "test_multict_gws_manhattan_points.tsv"), "\n",
    "  ", SUGG_CACHE, "\n", sep = "")
