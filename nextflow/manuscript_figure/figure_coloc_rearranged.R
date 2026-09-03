#!/usr/bin/env Rscript
# ============================================================
# Colocalization multi-panel figure (rearranged layout)
# ============================================================
# Layout:
#   A  full-width  Cell type × disease landscape + marginals
#   B  full-width  Cell type × locus–disease architecture
#   C|D|E          Regional examples (VIP–TMEM106B–MDD,
#                  L5.6 IT Car3–CACNA1C–BD, Microglia–RP11-347L18.1–BD)
#
# Unit of a coloc event = independent locus–disease comparison
# (not individual SNPs). Same genomic locus × different diseases
# count as separate events but one unique locus.
#
# Outputs (manuscript_figure/):
#   figure_coloc_rearranged.{pdf,png}
#   panel_A_landscape.tsv
#   panel_B_locus_disease_events.tsv
#   coloc_summary_by_cell_type.tsv
#   panel_{A,B,C,D,E}_*.pdf
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(scales)
  library(grid)
  library(topr)
})

# ── paths ────────────────────────────────────────────────────────────────────
ROOT        <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
COLOC_FILE  <- file.path(ROOT, "results/coloc/full/coloc_all_results.tsv")
LOCI_DIR    <- file.path(ROOT, "results/coloc/loci_full")
DISEASE_DIR <- file.path(ROOT, "results/coloc/disease_gwas")
OUT_DIR     <- file.path(ROOT, "manuscript_figure")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PP_HIT   <- 0.5
PP_STRONG <- 0.8
# Classic LocusZoom flank: ±250 kb around the CTP lead.
WINDOW_KB <- 250L
# Coloc inputs are hg19 (disease files *_hg19.tsv; 05_regional_coloc_plots.R
# defaults genome_build=37). Gene tracks must use build 37 — build 38 mis-annotates
# the Microglia chr6 lead (C6orf118 instead of RP11-347L18.1).
GENOME_BUILD <- 37L
LABEL_SIZE <- 15

# Disease order and label map (primary updated disease GWAS only).
DISEASE_ORDER <- c("MDD", "BD", "SCZ", "AD", "LBD", "PD")
DISEASE_MAP <- c(
  MDD_MDD2025        = "MDD",
  BD_bip2024         = "BD",
  SCZ_Trubetskoy2022 = "SCZ",
  AD_Bellenguez2022  = "AD",
  LBD_Chia2021       = "LBD",
  PD_Nalls2019       = "PD"
)
DISEASES_IN_FIGURE <- DISEASE_ORDER

# Fixed publication cell-type order (Non-neuronal → Inhibitory → Excitatory)
CT_ORDER_RAW <- c(
  "Microglia", "Endothelial", "Oligodendrocyte", "Pericyte", "VLMC",
  "Astrocyte", "OPC",
  "LAMP5", "VIP", "PAX6", "PVALB", "SST",
  "IT", "L5.6.IT.Car3", "L5.6.NP", "L5.ET", "L6.CT", "L4.IT", "L6b"
)
CT_CLASS <- c(
  Microglia = "Non-neuronal", Endothelial = "Non-neuronal",
  Oligodendrocyte = "Non-neuronal", Pericyte = "Non-neuronal",
  VLMC = "Non-neuronal", Astrocyte = "Non-neuronal", OPC = "Non-neuronal",
  LAMP5 = "Inhibitory", VIP = "Inhibitory", PAX6 = "Inhibitory",
  PVALB = "Inhibitory", SST = "Inhibitory",
  IT = "Excitatory", `L5.6.IT.Car3` = "Excitatory", `L5.6.NP` = "Excitatory",
  `L5.ET` = "Excitatory", `L6.CT` = "Excitatory", `L4.IT` = "Excitatory",
  L6b = "Excitatory"
)
FOCUS_CTS <- c("VIP", "L5.6.IT.Car3", "Microglia")
FOCUS_COL <- "#c0392b"

# Fixed biotype colours shared by all gene tracks + the C–E legend.
# topr otherwise assigns ggplot hue colours *per panel*, so the same biotype
# (e.g. protein_coding) can be red in C, green in D, and the hardcoded legend
# (teal/lncRNA/purple) does not match what is drawn.
BIOTYPE_COLORS <- c(
  protein_coding  = "#0072B2",  # blue
  lincRNA        = "#D55E00",  # vermillion
  lncRNA         = "#D55E00",
  antisense      = "#E69F00",  # orange
  sense_intronic = "#009E73",  # green
  pseudogene     = "#CC79A7",  # reddish purple
  processed_pseudogene = "#CC79A7"
)
BIOTYPE_LEGEND_ORDER <- c(
  "protein_coding", "lincRNA", "lncRNA", "antisense",
  "sense_intronic", "pseudogene", "processed_pseudogene"
)

# Clean display labels (no underscores; periods → spaces in subtype names)
ct_display <- function(ct) {
  map <- c(
    "L5.6.IT.Car3" = "L5.6 IT Car3",
    "L5.6.NP" = "L5.6 NP",
    "L5.ET" = "L5 ET",
    "L6.CT" = "L6 CT",
    "L4.IT" = "L4 IT"
  )
  ifelse(ct %in% names(map), map[ct], ct)
}

# Curated gene / region labels for significant coloc loci (hg19 / GRCh37)
# Fallback: "chrN locus"
assign_locus_label <- function(chr, pos) {
  chr <- as.character(chr)
  pos <- as.numeric(pos)
  dplyr::case_when(
    chr == "7"  & pos >= 1.20e7 & pos <= 1.30e7 ~ "TMEM106B",
    chr == "12" & pos >= 1.8e6  & pos <= 3.0e6  ~ "CACNA1C",
    # Microglia BD lead 6:164862615 — nearest gene on hg19 is RP11-347L18.1
    # (lincRNA ~164.77 Mb). Not PRKN (~2 Mb upstream) and not hg38 C6orf118.
    chr == "6"  & pos >= 1.64e8 & pos <= 1.66e8 ~ "RP11-347L18.1 region",
    chr == "6"  & pos >= 1.61e8 & pos <  1.64e8 ~ "PRKN region",
    chr == "11" & pos >= 1.4e7  & pos <= 1.8e7  ~ "chr11 locus",
    chr == "18" & pos >= 6.0e7  & pos <= 6.3e7  ~ "chr18 locus",
    chr == "4"  & pos >= 1.07e8 & pos <= 1.11e8 ~ "chr4 locus",
    TRUE ~ paste0("chr", chr, " locus")
  )
}

# Shared PP.H4 fill scale (0–1), matching the original coloc figure
pp_h4_fill_scale <- function(name = "Max PP.H4") {
  scale_fill_gradientn(
    name   = name,
    colours = c("grey97", "#FEE8C8", "#FDBB84", "#E34A33", "#8B0000"),
    values  = rescale(c(0, 0.2, 0.5, 0.8, 1.0)),
    limits  = c(0, 1),
    breaks  = c(0, 0.25, 0.5, 0.75, 1.0),
    na.value = "grey90",
    oob = squish
  )
}

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8, hjust = 0),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      complete = TRUE
    )
}

# ============================================================
# DATA PREPARATION
# ============================================================
cat("=== Loading coloc results ===\n")
coloc_raw <- fread(COLOC_FILE)
coloc <- copy(coloc_raw)
coloc[, disease_clean := DISEASE_MAP[disease]]
coloc <- coloc[!is.na(disease_clean)]
coloc <- coloc[disease_clean %in% DISEASES_IN_FIGURE]
coloc[, PP.H4 := as.numeric(PP.H4)]
# PD is mostly untested (grey tiles); kept in figure at user request.
n_pd <- nrow(coloc[disease_clean == "PD"])
cat("NOTE: PD included;", n_pd,
    "CTP×PD tests with enough SNPs (most cell types untested → grey).\n")

# Parse genomic coordinates from locus_id: {cell_type}_chr{N}_{pos}
coloc[, genomic_tag := sub("^.+_(chr[0-9XY]+_[0-9]+)$", "\\1", locus_id)]
coloc[, chr := sub("^chr", "", sub("_[0-9]+$", "", genomic_tag))]
coloc[, lead_pos := as.integer(sub("^chr[0-9XY]+_", "", genomic_tag))]
coloc[, locus_label := assign_locus_label(chr, lead_pos)]
coloc[, ct_class := CT_CLASS[cell_type]]
coloc[is.na(ct_class), ct_class := "Other"]

unmatched_ct <- setdiff(unique(coloc$cell_type), CT_ORDER_RAW)
unmatched_dis <- setdiff(unique(coloc_raw$disease), names(DISEASE_MAP))
if (length(unmatched_ct)) {
  warning("Unmatched cell-type labels (not in CT_ORDER_RAW): ",
          paste(unmatched_ct, collapse = ", "))
}
if (length(unmatched_dis)) {
  cat("NOTE: Diseases excluded from figure (not in primary disease set):\n  ",
      paste(unmatched_dis, collapse = ", "), "\n")
}

# Factor levels: display order top→bottom = CT_ORDER_RAW; ggplot y = rev
ct_levels_plot <- rev(CT_ORDER_RAW)          # bottom → top on y
ct_label_levels <- rev(ct_display(CT_ORDER_RAW))
coloc[, cell_type_f := factor(cell_type, levels = ct_levels_plot)]
coloc[, cell_label  := factor(ct_display(cell_type), levels = ct_label_levels)]
coloc[, disease_clean := factor(disease_clean, levels = DISEASE_ORDER)]

# ── Independent locus identity ───────────────────────────────────────────────
# Cluster lead positions within 1 Mb on the same chromosome into one
# independent locus. Multiple diseases at the same locus remain separate
# *events* but share one locus identity for "unique loci" counts.
cluster_independent_loci <- function(dt, window = 1e6) {
  x <- copy(dt)
  x[, indep_locus_id := NA_character_]
  for (ch in unique(x$chr)) {
    idx <- which(x$chr == ch)
    ord <- idx[order(x$lead_pos[idx])]
    if (!length(ord)) next
    cluster_id <- 1L
    cluster_start <- x$lead_pos[ord[1]]
    for (i in ord) {
      if (x$lead_pos[i] - cluster_start > window) {
        cluster_id <- cluster_id + 1L
        cluster_start <- x$lead_pos[i]
      }
      x$indep_locus_id[i] <- sprintf("chr%s_c%d", ch, cluster_id)
    }
  }
  x
}

coloc <- cluster_independent_loci(coloc)

# Significant events (PP.H4 ≥ 0.5)
sig <- coloc[PP.H4 >= PP_HIT]
if (anyDuplicated(sig, by = c("cell_type", "disease_clean", "indep_locus_id"))) {
  stop("Duplicated cell type–disease–locus events after clustering — inspect data.")
}
cat(sprintf("Significant events (PP.H4 >= %.1f): %d\n", PP_HIT, nrow(sig)))

# ── Panel A landscape table ──────────────────────────────────────────────────
# Max PP.H4 per cell type × disease; complete grid for missing/untested
tested_pairs <- unique(coloc[, .(cell_type, disease_clean)])
tested_pairs[, tested := TRUE]

heat <- coloc[, .(
  max_h4 = max(PP.H4, na.rm = TRUE),
  n_tests = .N,
  n_indep_loci_tested = uniqueN(indep_locus_id)
), by = .(cell_type, disease_clean, ct_class)]

grid_A <- CJ(
  cell_type = CT_ORDER_RAW,
  disease_clean = factor(DISEASE_ORDER, levels = DISEASE_ORDER),
  unique = TRUE
)
grid_A[, ct_class := CT_CLASS[cell_type]]
grid_A <- merge(grid_A, heat,
                by = c("cell_type", "disease_clean", "ct_class"), all.x = TRUE)
grid_A[, tested := !is.na(n_tests) & n_tests > 0]
grid_A[, fill_status := fifelse(
  !tested, "untested",
  fifelse(is.finite(max_h4) & max_h4 >= PP_HIT, "hit", "below")
)]
# Plot fill: only hits use continuous PP.H4; others NA (handled via overlays)
grid_A[, fill_h4 := fifelse(fill_status == "hit", max_h4, NA_real_)]
grid_A[, cell_type_f := factor(cell_type, levels = ct_levels_plot)]
grid_A[, cell_label  := factor(ct_display(cell_type), levels = ct_label_levels)]
grid_A[, disease_clean := factor(disease_clean, levels = DISEASE_ORDER)]
grid_A[, strong := fill_status == "hit" & max_h4 >= PP_STRONG]

# ── Marginal counts ──────────────────────────────────────────────────────────
# Right bars: independent CTP-GWAS loci tested per cell type (from loci files
# when available; otherwise unique locus_id in coloc results).
ctp_loci_counts <- rbindlist(lapply(CT_ORDER_RAW, function(ct) {
  f <- file.path(LOCI_DIR, paste0(ct, "_loci.tsv"))
  if (file.exists(f)) {
    n <- nrow(fread(f))
    data.table(cell_type = ct, n_ctp_loci = n, source = "loci_file")
  } else {
    data.table(
      cell_type = ct,
      n_ctp_loci = uniqueN(coloc[cell_type == ct]$locus_id),
      source = "coloc_results"
    )
  }
}))

# Top bars: independent *disease GWAS* loci eligible for coloc.
# This pipeline tests CTP-GWAS loci against full disease summary statistics
# in each window — it does not start from a disease-locus list. True
# "disease loci tested" counts are therefore unavailable.
disease_loci_available <- FALSE
cat("\nWARNING: Eligible independent disease-GWAS-locus counts cannot be ",
    "derived from the current CTP-locus-driven coloc tables. ",
    "Panel A will omit the top marginal bar ('Disease loci tested') ",
    "while preserving the heatmap, right CTP-loci bar, and Events/Loci ",
    "annotations.\n\n", sep = "")

# Events / unique loci per cell type (PP.H4 ≥ 0.5)
# Events = cell type–disease–locus triples
# Unique loci = distinct independent genomic loci (across diseases)
ct_summary <- rbindlist(lapply(CT_ORDER_RAW, function(ct) {
  s <- sig[cell_type == ct]
  data.table(
    cell_type = ct,
    cell_label = ct_display(ct),
    ct_class = CT_CLASS[[ct]],
    n_events = nrow(s),
    n_unique_loci = uniqueN(s$indep_locus_id),
    diseases = if (nrow(s)) paste(sort(unique(as.character(s$disease_clean))),
                                  collapse = ";") else "",
    locus_labels = if (nrow(s)) paste(sort(unique(s$locus_label)),
                                      collapse = ";") else ""
  )
}))

# Validation checks (expected from analysis narrative)
expected_checks <- list(
  list(ct = "L5.6.IT.Car3", events = 3L, loci = 1L),
  list(ct = "LAMP5",        events = 2L, loci = 1L),
  list(ct = "Endothelial",  events = 2L, loci = 2L),
  list(ct = "IT",           events = 2L, loci = 2L),
  list(ct = "Microglia",    events = 2L, loci = 2L),
  list(ct = "Oligodendrocyte", events = 2L, loci = 2L),
  list(ct = "VIP",          events = 1L, loci = 1L),
  list(ct = "SST",          events = 0L, loci = 0L)
)
cat("=== Validation: events / unique loci ===\n")
for (chk in expected_checks) {
  row <- ct_summary[cell_type == chk$ct]
  ok <- row$n_events == chk$events && row$n_unique_loci == chk$loci
  cat(sprintf("  %s: events=%d (expect %d), loci=%d (expect %d) %s\n",
              chk$ct, row$n_events, chk$events,
              row$n_unique_loci, chk$loci,
              if (ok) "OK" else "MISMATCH"))
  if (!ok) warning("Validation mismatch for ", chk$ct)
}
print(ct_summary[, .(cell_type, n_events, n_unique_loci, diseases)])

# Per cell-type × disease: number of independent loci with hits (for checks)
cat("\n=== Independent loci per cell type × disease (PP.H4 >= 0.5) ===\n")
print(sig[, .(n_indep_loci = uniqueN(indep_locus_id),
              loci = paste(sort(unique(locus_label)), collapse = ",")),
          by = .(cell_type, disease_clean)][order(cell_type, disease_clean)])

# ── Panel B events table ─────────────────────────────────────────────────────
# One row per significant cell type–locus–disease event (no aggregation)
panel_B <- copy(sig)
panel_B[, cell_label := ct_display(cell_type)]
# Column key = disease × locus label; order loci within disease by max PP.H4
locus_rank <- panel_B[, .(max_pp = max(PP.H4)),
                      by = .(disease_clean, locus_label, indep_locus_id)]
locus_rank <- locus_rank[order(match(disease_clean, DISEASE_ORDER), -max_pp)]
locus_rank[, col_key := paste(disease_clean, locus_label, sep = "|")]
locus_rank[, col_key := factor(col_key, levels = unique(col_key))]
panel_B <- merge(panel_B, locus_rank,
                 by = c("disease_clean", "locus_label", "indep_locus_id"))
panel_B[, cell_type_f := factor(cell_type, levels = ct_levels_plot)]
panel_B[, cell_label_f := factor(ct_display(cell_type), levels = ct_label_levels)]
panel_B[, strong := PP.H4 >= PP_STRONG]

# ── Effect-direction harmonization for Panel B events ────────────────────────
# For each significant coloc event, pick the representative shared variant
# (smallest CTP meta p among allele-matched, strand-unambiguous SNPs in the
# locus window) and record whether the CTP-increasing allele increases or
# decreases disease risk.
#
# Allele conventions (verified against the pipeline scripts):
#   CTP locus data : `beta` is the effect of METAL `Allele1` (NOT always the
#                    ALT of the snp ID) -> effect allele = toupper(Allele1).
#   Disease files  : 02_download_disease_gwas.sh writes effect_allele into the
#                    `ref` column and other_allele into `alt`, so `beta` is
#                    the effect of the `ref` column allele.
cat("Computing effect directions (CTP-increasing allele vs disease risk)...\n")

.dir_events <- unique(panel_B[, .(cell_type, locus_id, disease, disease_clean,
                                  chr, lead_pos)])

# Locus windows from the per-cell-type lead tables
.lead_windows <- rbindlist(lapply(unique(.dir_events$cell_type), function(ct) {
  f <- file.path(LOCI_DIR, paste0(ct, "_loci.tsv"))
  if (!file.exists(f)) return(NULL)
  fread(f)[, .(cell_type, locus_id, window_start, window_end)]
}), use.names = TRUE)
.dir_events <- merge(.dir_events, .lead_windows,
                     by = c("cell_type", "locus_id"), all.x = TRUE)

# CTP meta sumstats within each event window
.ct_window_dat <- rbindlist(lapply(unique(.dir_events$cell_type), function(ct) {
  f <- file.path(LOCI_DIR, paste0(ct, "_locus_data.tsv.gz"))
  d <- fread(f)
  evs <- .dir_events[cell_type == ct]
  rbindlist(lapply(seq_len(nrow(evs)), function(i) {
    e <- evs[i]
    sub <- d[as.character(chr) == as.character(e$chr) &
             pos >= e$window_start & pos <= e$window_end,
             .(snp, pos, ref = toupper(ref), alt = toupper(alt),
               ea_ct = toupper(Allele1), oa_ct = toupper(Allele2),
               beta_ct = beta, p_ct = p)]
    if (nrow(sub) == 0) return(NULL)
    cbind(cell_type = ct, locus_id = e$locus_id, sub)
  }), use.names = TRUE)
}), use.names = TRUE)

# Disease sumstats: one awk pass per disease over all its event windows
.dis_window_dat <- rbindlist(lapply(unique(.dir_events$disease), function(dk) {
  f <- file.path(DISEASE_DIR, dk, paste0(dk, "_hg19.tsv"))
  if (!file.exists(f)) return(NULL)
  wins <- unique(.dir_events[disease == dk,
                             .(chr, window_start, window_end)])
  cond <- paste(sprintf('($2=="%s" && $3>=%d && $3<=%d)',
                        wins$chr, wins$window_start, wins$window_end),
                collapse = " || ")
  hdr <- names(fread(f, nrows = 0L))
  d <- tryCatch(
    fread(cmd = sprintf("awk -F'\t' 'NR>1 && (%s)' %s", cond, shQuote(f)),
          header = FALSE, col.names = hdr),
    error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  # `ref` column holds the disease effect allele (see header comment above)
  d[, .(disease = dk, chr = as.character(chr), pos = as.integer(pos),
        ea_dis = toupper(ref), oa_dis = toupper(alt),
        beta_dis = as.numeric(beta), p_dis = as.numeric(p))]
}), use.names = TRUE)

.compute_direction <- function(ev) {
  ctd <- .ct_window_dat[cell_type == ev$cell_type & locus_id == ev$locus_id]
  dsd <- .dis_window_dat[disease == ev$disease &
                         chr == as.character(ev$chr) &
                         pos >= ev$window_start & pos <= ev$window_end]
  if (nrow(ctd) == 0 || nrow(dsd) == 0) return(NULL)
  m <- merge(ctd, dsd, by = "pos", allow.cartesian = TRUE)
  if (nrow(m) == 0) return(NULL)
  # allele-set match + drop strand-ambiguous pairs
  m <- m[((ea_ct == ea_dis & oa_ct == oa_dis) |
          (ea_ct == oa_dis & oa_ct == ea_dis))]
  m <- m[!(paste(ea_ct, oa_ct) %in% c("A T", "T A", "C G", "G C"))]
  if (nrow(m) == 0) return(NULL)
  m[, beta_dis_aligned := fifelse(ea_dis == ea_ct, beta_dis, -beta_dis)]
  m <- m[is.finite(beta_ct) & is.finite(beta_dis_aligned)]
  if (nrow(m) == 0) return(NULL)
  best <- m[which.min(p_ct)]
  # orient to the CTP-increasing allele
  flip <- best$beta_ct < 0
  data.table(
    cell_type    = ev$cell_type,
    locus_id     = ev$locus_id,
    disease      = ev$disease,
    rep_snp      = best$snp,
    ctp_inc_allele = if (flip) best$oa_ct else best$ea_ct,
    beta_ctp_inc = abs(best$beta_ct),
    beta_disease_per_ctp_inc = if (flip) -best$beta_dis_aligned
                               else best$beta_dis_aligned,
    p_ct  = best$p_ct,
    p_dis = best$p_dis
  )
}

dir_tbl <- rbindlist(lapply(seq_len(nrow(.dir_events)), function(i)
  .compute_direction(.dir_events[i])), use.names = TRUE)
dir_tbl[, risk_up := beta_disease_per_ctp_inc > 0]

panel_B <- merge(panel_B,
                 dir_tbl[, .(cell_type, locus_id, disease, rep_snp,
                             ctp_inc_allele, beta_ctp_inc,
                             beta_disease_per_ctp_inc, p_dis, risk_up)],
                 by = c("cell_type", "locus_id", "disease"), all.x = TRUE)
panel_B[, dir_glyph := fifelse(is.na(risk_up), "",
                        fifelse(risk_up, "\u25B2", "\u25BC"))]

cat("Direction annotation coverage:",
    sum(!is.na(panel_B$risk_up)), "/", nrow(panel_B), "events\n")
print(panel_B[!is.na(risk_up),
              .(cell_type, disease_clean, locus_label, rep_snp,
                ctp_inc_allele, beta_ctp_inc,
                beta_disease_per_ctp_inc, p_dis, risk_up)])

# Export tables
panel_A_export <- grid_A[, .(
  cell_type, cell_label = as.character(cell_label), ct_class,
  disease = as.character(disease_clean),
  max_PP_H4 = max_h4, n_tests, tested, fill_status, strong
)][order(match(cell_type, CT_ORDER_RAW), match(disease, DISEASE_ORDER))]
fwrite(panel_A_export, file.path(OUT_DIR, "panel_A_landscape.tsv"), sep = "\t")

panel_B_export <- panel_B[, .(
  cell_type, cell_label, ct_class, disease = as.character(disease_clean),
  locus_id, indep_locus_id, locus_label, chr, lead_pos,
  PP.H4, n_snps, strong, col_key = as.character(col_key),
  rep_snp, ctp_inc_allele, beta_ctp_inc, beta_disease_per_ctp_inc,
  p_disease_rep_snp = p_dis, risk_up
)][order(match(disease, DISEASE_ORDER), -PP.H4, cell_type)]
fwrite(panel_B_export, file.path(OUT_DIR, "panel_B_locus_disease_events.tsv"),
       sep = "\t")

summary_export <- merge(ct_summary, ctp_loci_counts, by = "cell_type")
summary_export <- summary_export[order(match(cell_type, CT_ORDER_RAW))]
fwrite(summary_export, file.path(OUT_DIR, "coloc_summary_by_cell_type.tsv"),
       sep = "\t")
cat("Wrote panel_A_landscape.tsv, panel_B_locus_disease_events.tsv, ",
    "coloc_summary_by_cell_type.tsv\n")

# ============================================================
# PANEL A — landscape heatmap + marginals
# ============================================================
cat("Building Panel A...\n")

# Separators + class labels centred on middle cell type of each block.
# y levels bottom→top: Exc (1–7), Inh (8–12), NN (13–19).
n_exc <- 7L
n_inh <- 5L
n_nn  <- 7L
sep_y <- c(n_exc + 0.5, n_exc + n_inh + 0.5)
y_expand_A <- expansion(mult = c(0.02, 0.02))

grid_A[, plot_h4 := fifelse(tested, pmax(0, pmin(1, max_h4)), NA_real_)]

y_colours <- ifelse(levels(grid_A$cell_label) %in% ct_display(FOCUS_CTS),
                    FOCUS_COL, "grey20")
y_faces <- ifelse(levels(grid_A$cell_label) %in% ct_display(FOCUS_CTS),
                  "bold", "plain")

# Anchor each class name on the middle cell-type row of that block
group_lab_dt <- data.table(
  label = c("Excitatory", "Inhibitory", "Non-neuronal"),
  cell_label = factor(
    c(ct_label_levels[as.integer((1 + n_exc) / 2)],
      ct_label_levels[n_exc + as.integer((1 + n_inh) / 2)],
      ct_label_levels[n_exc + n_inh + as.integer((1 + n_nn) / 2)]),
    levels = ct_label_levels
  )
)

p_heat <- ggplot(grid_A, aes(x = disease_clean, y = cell_label)) +
  geom_tile(data = grid_A[tested == FALSE],
            fill = "grey78", colour = "white", linewidth = 0.4) +
  geom_tile(data = grid_A[tested == TRUE],
            aes(fill = plot_h4), colour = "white", linewidth = 0.4) +
  geom_text(data = grid_A[strong == TRUE],
            aes(label = "\u2605"), colour = "white", size = 2.8,
            fontface = "bold", vjust = 0.45) +
  pp_h4_fill_scale("Max PP.H4") +
  geom_hline(yintercept = sep_y, colour = "grey40", linewidth = 0.4) +
  scale_x_discrete(position = "top", drop = FALSE, limits = DISEASE_ORDER) +
  scale_y_discrete(drop = FALSE, expand = y_expand_A) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_pub(9) +
  theme(
    axis.text.x = element_text(size = 8, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 7.5, colour = y_colours, face = y_faces),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.28, "cm"),
    plot.margin = margin(8, 4, 4, 4)
  )

# Left strip: keep a real (invisible) y-axis so cowplot align="h" matches the
# heatmap panel geometry — theme_void strips drift vertically.
p_groups <- ggplot(
  data.table(cell_label = factor(ct_label_levels, levels = ct_label_levels)),
  aes(x = 1, y = cell_label)
) +
  geom_blank() +
  geom_hline(yintercept = sep_y, colour = "grey70", linewidth = 0.3) +
  geom_text(
    data = group_lab_dt,
    aes(x = 1, y = cell_label, label = label),
    angle = 90, size = 2.55, fontface = "italic", colour = "grey35",
    hjust = 0.5, vjust = 0.5
  ) +
  scale_y_discrete(drop = FALSE, expand = y_expand_A) +
  # Invisible top x labels ≈ heatmap disease headers (same vertical panel box)
  scale_x_continuous(
    limits = c(0.5, 1.5), breaks = 1, labels = "", position = "top"
  ) +
  coord_cartesian(clip = "off") +
  theme_pub(9) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, colour = NA),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(8, 0, 4, 2)
  )

# Right bars: # CTP-GWAS loci that entered coloc testing per cell type
ctp_plot_dt <- merge(
  data.table(cell_type = CT_ORDER_RAW,
             cell_label = factor(ct_display(CT_ORDER_RAW),
                                 levels = ct_label_levels)),
  ctp_loci_counts, by = "cell_type"
)
p_ctp <- ggplot(ctp_plot_dt, aes(x = n_ctp_loci, y = cell_label)) +
  geom_col(fill = "grey45", width = 0.7) +
  geom_hline(yintercept = sep_y, colour = "grey70", linewidth = 0.3) +
  scale_y_discrete(drop = FALSE, expand = y_expand_A) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "CTP loci tested", y = NULL) +
  theme_pub(8) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(size = 7),
    plot.margin = margin(8, 2, 4, 0)
  )

# Narrow Events / Loci annotation columns
ann_dt <- merge(
  data.table(cell_type = CT_ORDER_RAW,
             cell_label = factor(ct_display(CT_ORDER_RAW),
                                 levels = ct_label_levels)),
  ct_summary[, .(cell_type, n_events, n_unique_loci)],
  by = "cell_type"
)
make_ann_col <- function(dt, value_col, title) {
  ggplot(dt, aes(x = 1, y = cell_label)) +
    geom_text(aes(label = .data[[value_col]]), size = 2.4,
              colour = "grey15", fontface = "plain") +
    geom_hline(yintercept = sep_y, colour = "grey80", linewidth = 0.25) +
    scale_y_discrete(drop = FALSE, expand = y_expand_A) +
    scale_x_continuous(limits = c(0.5, 1.5), breaks = 1, labels = title,
                       position = "top") +
    labs(x = NULL, y = NULL) +
    theme_pub(7) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_text(size = 7, face = "bold", colour = "grey25"),
      plot.margin = margin(8, 1, 4, 0)
    )
}
p_ev  <- make_ann_col(ann_dt, "n_events", "Events")
p_loc <- make_ann_col(ann_dt, "n_unique_loci", "Loci")

heat_legend <- cowplot::get_legend(
  p_heat +
    labs(fill = "Max PP.H4") +
    theme(
      legend.position = "right",
      legend.key.height = unit(0.55, "cm"),
      legend.key.width = unit(0.28, "cm")
    )
)
p_heat_nl <- p_heat + theme(legend.position = "none")
# Align data panels first (legend grob breaks cowplot axis alignment)
panel_A_body <- plot_grid(
  p_groups, p_heat_nl, p_ctp, p_ev, p_loc,
  nrow = 1,
  rel_widths = c(0.40, 3.6, 1.15, 0.42, 0.42),
  align = "h", axis = "tb"
)
panel_A_core <- plot_grid(
  panel_A_body, heat_legend,
  nrow = 1, rel_widths = c(1, 0.14)
)

# ============================================================
# PANEL B — locus–disease architecture
# ============================================================
cat("Building Panel B...\n")

# Disease strip annotation data
col_meta <- unique(locus_rank[, .(col_key, disease_clean, locus_label)])
col_meta[, col_key := factor(col_key, levels = levels(locus_rank$col_key))]
# Vertical separators between disease groups
dis_bounds <- col_meta[, .(n = .N), by = disease_clean]
dis_bounds[, end := cumsum(n)]
dis_bounds[, start := end - n + 1]
sep_x <- dis_bounds$end[-nrow(dis_bounds)] + 0.5

y_expand_B <- expansion(add = 0.55)

p_arch <- ggplot(panel_B, aes(x = col_key, y = cell_label_f)) +
  geom_tile(
    data = CJ(
      col_key = factor(levels(panel_B$col_key), levels = levels(panel_B$col_key)),
      cell_label_f = factor(ct_label_levels, levels = ct_label_levels)
    ),
    fill = NA, colour = "grey94", linewidth = 0.15
  ) +
  geom_point(aes(fill = PP.H4),
             shape = 21, size = 2.6, colour = "grey40", stroke = 0.3) +
  geom_point(data = panel_B[strong == TRUE],
             aes(fill = PP.H4),
             shape = 21, size = 2.6, colour = "grey10", stroke = 0.8) +
  geom_text(data = panel_B[dir_glyph != ""],
            aes(label = dir_glyph,
                colour = fifelse(PP.H4 >= 0.75, "white", "grey15")),
            size = 1.4, vjust = 0.42, show.legend = FALSE) +
  scale_colour_identity() +
  pp_h4_fill_scale("PP.H4") +
  geom_vline(xintercept = sep_x, colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = sep_y, colour = "grey55", linewidth = 0.35) +
  scale_x_discrete(
    labels = function(x) sub("^[^|]+\\|", "", x),
    expand = expansion(add = 0.55)
  ) +
  scale_y_discrete(drop = FALSE, expand = y_expand_B) +
  labs(x = NULL, y = NULL,
       caption = paste0("\u25B2 / \u25BC  CTP-increasing allele increases / ",
                        "decreases disease risk (representative shared variant)")) +
  coord_cartesian(clip = "off") +
  theme_pub(8) +
  theme(
    axis.text.x = element_text(size = 6.5, angle = 40, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7, colour = y_colours, face = y_faces,
                               lineheight = 1.05),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.caption = element_text(size = 6, colour = "grey30", hjust = 0,
                                margin = margin(t = 4)),
    plot.margin = margin(10, 8, 8, 4)
  )

# Disease group labels above (generous height to avoid cropping)
dis_lab <- dis_bounds[, .(
  disease_clean,
  x = (start + end) / 2
)]
p_dis_strip <- ggplot(dis_lab, aes(x = x, y = 1, label = disease_clean)) +
  geom_text(fontface = "bold", size = 3.2, vjust = 0.5) +
  scale_x_continuous(
    limits = c(0.5, nrow(col_meta) + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.2, 1.8)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(6, 8, 0, 4))

# Class labels share Panel B discrete y; blank spacer = disease-strip height.
# Keep invisible axes (not theme_void) so align="h" locks panel geometry.
p_groups_B <- ggplot(
  data.table(cell_label = factor(ct_label_levels, levels = ct_label_levels)),
  aes(x = 1, y = cell_label)
) +
  geom_blank() +
  geom_hline(yintercept = sep_y, colour = "grey70", linewidth = 0.3) +
  geom_text(
    data = group_lab_dt,
    aes(x = 1, y = cell_label, label = label),
    angle = 90, size = 2.4, fontface = "italic", colour = "grey35",
    hjust = 0.5, vjust = 0.5
  ) +
  scale_y_discrete(drop = FALSE, expand = y_expand_B) +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = 1, labels = "") +
  coord_cartesian(clip = "off") +
  theme_pub(8) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    # Match arch x-axis label band height (angled locus names)
    axis.text.x = element_text(size = 6.5, angle = 40, colour = NA, hjust = 1),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(10, 0, 8, 2)
  )

p_arch_nl <- p_arch + theme(legend.position = "none")
dis_row <- plot_grid(NULL, p_dis_strip, nrow = 1, rel_widths = c(0.10, 1))
arch_row <- plot_grid(p_groups_B, p_arch_nl, nrow = 1, rel_widths = c(0.10, 1),
                      align = "h", axis = "tb")
panel_B_core <- plot_grid(dis_row, arch_row, ncol = 1, rel_heights = c(0.16, 1))

# ============================================================
# PANELS C–E — regional colocalization (reuse 05 logic)
# ============================================================
cat("Building regional panels C–E...\n")

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

find_disease_file <- function(disease_label) {
  matches <- list.files(
    DISEASE_DIR,
    pattern = paste0("^", disease_label, "_hg19\\.tsv$"),
    recursive = TRUE, full.names = TRUE
  )
  if (!length(matches)) stop("Disease GWAS not found: ", disease_label)
  matches[1]
}

load_disease_window <- function(dis_file, locus_chr, win_start, win_end) {
  cmd <- sprintf(
    paste0(
      "awk -F'\\t' 'NR==1 {print; next} ",
      "($2==\"%s\" || $2==%s) && ($3+0)>=%d && ($3+0)<=%d {print}' %s"
    ),
    locus_chr, locus_chr, win_start, win_end, shQuote(dis_file)
  )
  dt <- fread(cmd = cmd, sep = "\t", showProgress = FALSE)
  if (!nrow(dt)) return(dt)
  dt[, `:=`(
    chr = as.character(chr),
    pos = as.integer(pos),
    p = as.numeric(p),
    ref = toupper(ref),
    alt = toupper(alt),
    neglog10p = neglog10p(p)
  )]
  dt
}

load_ct_window <- function(cell_type, locus_chr, win_start, win_end) {
  f <- file.path(LOCI_DIR, paste0(cell_type, "_locus_data.tsv.gz"))
  if (!file.exists(f)) stop("Missing locus data: ", f)
  cat(sprintf("  Loading %s locus data...\n", cell_type))
  dt <- fread(f, sep = "\t", showProgress = FALSE,
              select = c("chr", "pos", "snp", "p", "ref", "alt", "beta", "locus_id"))
  dt <- dt[as.character(chr) == as.character(locus_chr) &
             pos >= win_start & pos <= win_end]
  dt[, `:=`(
    chr = as.character(chr),
    pos = as.integer(pos),
    p = as.numeric(p),
    ref = toupper(ref),
    alt = toupper(alt),
    beta = as.numeric(beta),
    neglog10p = neglog10p(p)
  )]
  dt
}

make_gwas_panel <- function(plot_data, lead_pos_mb, lead_points,
                            show_x_axis = FALSE, chr_label = NULL) {
  # No subplot titles/headers — panel letter + shared SNP legend only.
  p <- ggplot(plot_data, aes(x = pos / 1e6, y = neglog10p)) +
    geom_vline(xintercept = lead_pos_mb, linetype = "dotted",
               linewidth = 0.35, colour = "grey40") +
    geom_point(aes(colour = in_coloc), alpha = 0.65, size = 0.9) +
    scale_colour_manual(
      values = c("FALSE" = "grey55", "TRUE" = "#c0392b"),
      labels = c("Other SNPs", "Harmonized SNPs"),
      name = NULL,
      drop = FALSE
    ) +
    labs(
      x = if (show_x_axis) {
        sprintf("Chromosome %s position (Mb)", chr_label)
      } else {
        NULL
      },
      y = expression(-log[10](p))
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_blank(),
      legend.position = "none",
      axis.title.x = if (show_x_axis) element_text(size = 7) else element_blank(),
      axis.text.x = if (show_x_axis) element_text(size = 6.5) else element_blank(),
      axis.title.y = element_text(size = 7),
      axis.text.y = element_text(size = 6.5),
      plot.margin = margin(2, 4, 1, 4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.2)
    )
  # data.table empty[1] yields a 1-row NA template — drop those.
  if (is.data.frame(lead_points) && nrow(lead_points) >= 1 &&
      !is.na(lead_points$pos[1]) && !is.na(lead_points$neglog10p[1])) {
    p <- p + geom_point(
      data = lead_points,
      aes(x = pos / 1e6, y = neglog10p),
      inherit.aes = FALSE,
      colour = "#1f78b4", size = 2.4, shape = 18
    )
  }
  p
}

prepare_topr_datasets <- function(ct_plot, lead_snp) {
  df <- as.data.frame(ct_plot) %>%
    mutate(
      CHROM = as.integer(chr),
      POS = as.integer(pos),
      P = as.numeric(p),
      ID = snp,
      Effect = if ("beta" %in% names(.)) as.numeric(beta) else 0
    ) %>%
    filter(!is.na(CHROM), !is.na(POS), !is.na(P), P > 0, P <= 1)
  if (!nrow(df)) return(NULL)
  df_pos <- df %>% filter(is.na(Effect) | Effect >= 0) %>%
    select(CHROM, POS, P, ID, Effect)
  df_neg <- df %>% filter(!is.na(Effect), Effect < 0) %>%
    select(CHROM, POS, P, ID, Effect)
  if (!nrow(df_neg)) return(list(df_pos))
  if (lead_snp %in% df_pos$ID) list(df_pos, df_neg) else list(df_neg, df_pos)
}

make_empty_gene_track <- function(locus_chr, lead_pos, region_size, note) {
  x_min <- (lead_pos - region_size / 2) / 1e6
  x_max <- (lead_pos + region_size / 2) / 1e6
  ggplot() +
    annotate("text", x = (x_min + x_max) / 2, y = 0.55, label = note,
             size = 2.7, colour = "grey30") +
    scale_x_continuous(
      name = sprintf("Chromosome %s position (Mb)", locus_chr),
      limits = c(x_min, x_max)
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_bw(base_size = 8) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 4, 2, 4)
    )
}

make_gene_track <- function(ct_plot, lead_snp, region_size, locus_chr, lead_pos,
                            empty_note = NULL) {
  if (!is.null(empty_note)) {
    return(make_empty_gene_track(locus_chr, lead_pos, region_size, empty_note))
  }

  datasets <- prepare_topr_datasets(ct_plot, lead_snp)
  if (is.null(datasets)) return(NULL)
  ds1 <- datasets[[1]]
  if (!lead_snp %in% ds1$ID) {
    lead_row <- as.data.frame(ct_plot) %>%
      filter(snp == lead_snp) %>%
      mutate(
        CHROM = as.integer(chr), POS = as.integer(pos),
        P = as.numeric(p), ID = snp,
        Effect = if ("beta" %in% names(.)) as.numeric(beta) else 0
      ) %>%
      select(CHROM, POS, P, ID, Effect)
    if (nrow(lead_row)) {
      datasets[[1]] <- bind_rows(lead_row, ds1)
    }
  }

  # Always keep biotype colouring (protein_coding / lncRNA / pseudogene).
  # Do NOT set protein_coding_only=TRUE — that forces every gene the same colour.
  #
  # Use region= (not variant=) so the gene track matches the association window
  # exactly (±WINDOW_KB). topr's variant= path can hard-cap the zoom.
  half <- as.integer(region_size / 2)
  gene_start <- as.integer(lead_pos) - half
  gene_end   <- as.integer(lead_pos) + half
  region_str <- sprintf("%s:%d-%d", locus_chr, gene_start, gene_end)
  cat(sprintf("  Gene track: region=%s build=%d\n", region_str, GENOME_BUILD))
  plots <- tryCatch(
    regionplot(
      datasets,
      region = region_str,
      build = as.integer(GENOME_BUILD),
      show_overview = FALSE,
      extract_plots = TRUE,
      title = NULL,
      legend_name = "biotype",
      show_gene_legend = TRUE,
      protein_coding_only = FALSE,
      gene_label_size = 1.8,
      angle = 0,
      max_genes = 30,
      vline = as.integer(lead_pos)
    ),
    error = function(e) NULL
  )

  if (is.null(plots) || is.null(plots$gene_plot)) {
    return(make_empty_gene_track(
      locus_chr, lead_pos, region_size,
      "No annotated genes in window"
    ))
  }

  present_biotypes <- tryCatch({
    sc <- ggplot_build(plots$gene_plot)$plot$scales$get_scales("fill")
    if (is.null(sc)) character(0) else as.character(sc$get_limits())
  }, error = function(e) character(0))
  # Keep only biotypes we have colours for; fall back to identity otherwise.
  fill_vals <- BIOTYPE_COLORS[intersect(names(BIOTYPE_COLORS), present_biotypes)]
  if (!length(fill_vals)) fill_vals <- BIOTYPE_COLORS

  gp <- plots$gene_plot +
    scale_fill_manual(values = fill_vals, drop = TRUE) +
    theme(
      plot.margin = margin(0, 4, 2, 4),
      axis.title.x = element_text(size = 7),
      axis.text = element_text(size = 5.5),
      legend.position = "none"
    )
  # Force smaller gene-name text (topr defaults are too large → overlaps)
  for (i in seq_along(gp$layers)) {
    geom_name <- class(gp$layers[[i]]$geom)[1]
    if (geom_name %in% c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel")) {
      gp$layers[[i]]$aes_params$size <- 1.6
    }
  }
  attr(gp, "biotypes") <- present_biotypes
  gp
}

build_regional_panel <- function(cell_type, locus_id, disease_raw,
                                 panel_letter, empty_gene_note = NULL,
                                 window_kb = WINDOW_KB,
                                 wide_gene_track = FALSE,
                                 gene_region_start = NULL,
                                 gene_region_end = NULL) {
  # Bind args to local names so data.table does not treat `locus_id == locus_id`
  # as a tautology against the column of the same name.
  target_locus   <- locus_id
  target_disease <- disease_raw
  target_ct      <- cell_type
  win_kb         <- as.integer(window_kb)

  hit <- coloc[locus_id == target_locus & disease == target_disease]
  if (!nrow(hit)) stop("Missing coloc hit: ", target_locus, " × ", target_disease)
  hit <- hit[1]
  pp_h4 <- hit$PP.H4

  loci <- fread(file.path(LOCI_DIR, paste0(target_ct, "_loci.tsv")))
  locus <- loci[locus_id == target_locus]
  if (!nrow(locus)) stop("Locus not in loci file: ", target_locus)
  locus <- locus[1]
  win_start <- max(1L, as.integer(locus$lead_pos) - win_kb * 1000L)
  win_end   <- as.integer(locus$lead_pos) + win_kb * 1000L
  locus_chr <- as.character(locus$chr)

  cat(sprintf("  Regional %s: %s × %s (CTP lead %s:%s, window ±%d kb, PP.H4=%.3f)\n",
              panel_letter, target_ct, target_disease,
              locus_chr, format(as.integer(locus$lead_pos), big.mark = ","),
              win_kb, pp_h4))
  ct_plot <- load_ct_window(target_ct, locus_chr, win_start, win_end)
  dis_file <- find_disease_file(target_disease)
  cat(sprintf("  Loading disease window %s...\n", target_disease))
  dis_plot <- load_disease_window(dis_file, locus_chr, win_start, win_end)
  if (!nrow(dis_plot)) stop("No disease SNPs in window for ", target_locus)

  shared <- harmonize_snps(
    as.data.frame(ct_plot[, .(snp, ref, alt)]),
    as.data.frame(dis_plot[, .(snp, ref, alt)])
  )
  ct_plot[, in_coloc := snp %in% shared]
  ct_plot[, lead := snp == locus$lead_snp]
  dis_plot[, in_coloc := snp %in% shared]
  dis_plot[, lead := pos == as.integer(locus$lead_pos)]

  # Match scripts/coloc/05_regional_coloc_plots.R marker convention:
  #   CTP diamond  = CTP locus lead SNP
  #   Disease diamond = strongest disease SNP in the plotted window
  #   Dotted line  = CTP lead position (shared reference)
  lead_ct <- ct_plot[snp == locus$lead_snp]
  if (!nrow(lead_ct)) lead_ct <- ct_plot[0]
  lead_dis <- dis_plot[order(p)]
  if (nrow(lead_dis)) lead_dis <- lead_dis[1]
  if (nrow(lead_ct) && nrow(lead_dis)) {
    cat(sprintf(
      "    Markers: CTP lead %s (−log10p=%.2f); disease top %s (−log10p=%.2f; %+d bp from CTP lead)\n",
      lead_ct$snp[1], lead_ct$neglog10p[1],
      lead_dis$snp[1], lead_dis$neglog10p[1],
      as.integer(lead_dis$pos[1]) - as.integer(locus$lead_pos)
    ))
  }
  lead_pos_mb <- locus$lead_pos / 1e6
  region_size <- win_end - win_start

  ct_lead_df <- if (nrow(lead_ct)) as.data.frame(lead_ct[1]) else data.frame()
  dis_lead_df <- if (nrow(lead_dis)) as.data.frame(lead_dis[1]) else data.frame()
  ct_panel <- make_gwas_panel(
    as.data.frame(ct_plot),
    lead_pos_mb = lead_pos_mb,
    lead_points = ct_lead_df,
    show_x_axis = FALSE,
    chr_label = locus_chr
  )
  dis_panel <- make_gwas_panel(
    as.data.frame(dis_plot),
    lead_pos_mb = lead_pos_mb,
    lead_points = dis_lead_df,
    show_x_axis = FALSE,
    chr_label = locus_chr
  )

  gene_panel <- tryCatch(
    make_gene_track(
      ct_plot, locus$lead_snp, region_size, locus_chr,
      as.integer(locus$lead_pos),
      empty_note = empty_gene_note
    ),
    error = function(e) {
      warning("Gene track failed for ", target_locus, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(gene_panel)) {
    body <- plot_grid(ct_panel, dis_panel, gene_panel, ncol = 1,
                      rel_heights = c(3, 3, 1.5), align = "v", axis = "lr")
  } else {
    dis_panel <- dis_panel +
      labs(x = sprintf("Chromosome %s position (Mb)", locus_chr)) +
      theme(axis.title.x = element_text(size = 7),
            axis.text.x = element_text(size = 6.5))
    body <- plot_grid(ct_panel, dis_panel, ncol = 1,
                      rel_heights = c(1, 1), align = "v", axis = "lr")
  }

  # Panel letter only — no title / subtitle / PP.H4 header
  labeled <- ggdraw(body) +
    draw_label(panel_letter, x = 0.01, y = 0.995, hjust = 0, vjust = 1,
               fontface = "bold", size = LABEL_SIZE)
  attr(labeled, "pp_h4") <- pp_h4
  attr(labeled, "n_harmonized") <- length(shared)
  attr(labeled, "gene_panel") <- gene_panel
  labeled
}

# Regional examples — classic LocusZoom ±250 kb, gene track hg19 (build 37).
# Panel E nearest gene on hg19: RP11-347L18.1 (~93 kb upstream of lead).
panel_C <- build_regional_panel(
  cell_type = "VIP",
  locus_id = "VIP_chr7_12284430",
  disease_raw = "MDD_MDD2025",
  panel_letter = "C"
)
panel_D <- build_regional_panel(
  cell_type = "L5.6.IT.Car3",
  locus_id = "L5.6.IT.Car3_chr12_2324042",
  disease_raw = "BD_bip2024",
  panel_letter = "D"
)
panel_E <- build_regional_panel(
  cell_type = "Microglia",
  locus_id = "Microglia_chr6_164862615",
  disease_raw = "BD_bip2024",
  panel_letter = "E"
)

# ============================================================
# ASSEMBLE + SAVE
# ============================================================
cat("Assembling full figure...\n")

panel_A_labeled <- ggdraw(panel_A_core) +
  draw_label("A", x = 0.005, y = 0.995, hjust = 0, vjust = 1,
             fontface = "bold", size = LABEL_SIZE)
# Extra top padding so Panel B disease labels / top row are not clipped
panel_B_labeled <- plot_grid(
  NULL,
  ggdraw(panel_B_core) +
    draw_label("B", x = 0.005, y = 0.99, hjust = 0, vjust = 1,
               fontface = "bold", size = LABEL_SIZE),
  ncol = 1, rel_heights = c(0.06, 1)
)

# Give Panel B more vertical space so 19 cell-type rows do not overlap
row_AB <- plot_grid(
  panel_A_labeled, panel_B_labeled,
  ncol = 1,
  rel_heights = c(0.22, 0.30) / (0.22 + 0.30)
)

# Shared SNP colour legend
snp_legend <- plot_grid(
  ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
    geom_point(colour = "grey55", size = 2.8) + theme_void(),
  ggdraw() + draw_label("Other SNPs", size = 8, hjust = 0, x = 0.05),
  ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
    geom_point(colour = "#c0392b", size = 2.8) + theme_void(),
  ggdraw() + draw_label("Harmonized SNPs", size = 8, hjust = 0, x = 0.05),
  nrow = 1, rel_widths = c(0.05, 0.18, 0.05, 0.22)
)

# Biotype legend: only biotypes actually drawn in C–E, colours = BIOTYPE_COLORS
biotypes_used <- unique(unlist(lapply(
  list(panel_C, panel_D, panel_E),
  function(p) {
    gp <- attr(p, "gene_panel")
    if (is.null(gp)) return(character(0))
    attr(gp, "biotypes")
  }
)))
biotypes_used <- intersect(BIOTYPE_LEGEND_ORDER, biotypes_used)
# Collapse alias labels for the legend key (lncRNA/lincRNA share a colour)
legend_biotypes <- biotypes_used
if ("lincRNA" %in% legend_biotypes && "lncRNA" %in% legend_biotypes) {
  legend_biotypes <- setdiff(legend_biotypes, "lncRNA")
}
if ("processed_pseudogene" %in% legend_biotypes && "pseudogene" %in% legend_biotypes) {
  legend_biotypes <- setdiff(legend_biotypes, "processed_pseudogene")
}
biotype_legend_parts <- lapply(legend_biotypes, function(bt) {
  list(
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_point(shape = 15, colour = BIOTYPE_COLORS[[bt]], size = 3.5) +
      theme_void(),
    ggdraw() + draw_label(bt, size = 7.5, hjust = 0, x = 0.02)
  )
})
biotype_legend_parts <- unlist(biotype_legend_parts, recursive = FALSE)
biotype_legend <- plot_grid(
  plotlist = biotype_legend_parts, nrow = 1,
  rel_widths = rep(c(0.035, 0.14), length(legend_biotypes))
)

cde_legends <- plot_grid(snp_legend, biotype_legend, nrow = 1, rel_widths = c(0.7, 1.3))

row_CDE_plots <- plot_grid(panel_C, panel_D, panel_E, ncol = 3,
                           rel_widths = c(1, 1, 1))
row_CDE <- plot_grid(row_CDE_plots, cde_legends, ncol = 1,
                     rel_heights = c(1, 0.07))

fig <- plot_grid(row_AB, row_CDE, ncol = 1,
                 rel_heights = c(0.52, 0.48))

# Page width ≈ 185 mm; taller page so Panel B rows are not stacked
FIG_W <- 185 / 25.4
FIG_H <- 340 / 25.4

out_pdf <- file.path(OUT_DIR, "figure_coloc_rearranged.pdf")
out_png <- file.path(OUT_DIR, "figure_coloc_rearranged.png")

cairo_pdf(out_pdf, width = FIG_W, height = FIG_H)
print(fig)
dev.off()

png(out_png, width = FIG_W, height = FIG_H, units = "in",
    res = 600, type = "cairo")
print(fig)
dev.off()

# Standalone panels
ggsave(file.path(OUT_DIR, "panel_A_landscape.pdf"), panel_A_labeled,
       width = FIG_W, height = FIG_W * 0.55, device = cairo_pdf)
panel_B_standalone <- plot_grid(
  panel_B_labeled,
  cowplot::get_legend(
    p_arch + theme(legend.position = "right",
                   legend.key.height = unit(0.5, "cm"))
  ),
  nrow = 1, rel_widths = c(1, 0.12)
)
ggsave(file.path(OUT_DIR, "panel_B_architecture.pdf"), panel_B_standalone,
       width = FIG_W, height = FIG_W * 0.42, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "panel_C_TMEM106B_MDD.pdf"), panel_C,
       width = FIG_W / 3 * 1.15, height = 6.5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "panel_D_CACNA1C_BD.pdf"), panel_D,
       width = FIG_W / 3 * 1.15, height = 6.5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "panel_E_RP11-347L18.1_BD.pdf"), panel_E,
       width = FIG_W / 3 * 1.15, height = 6.5, device = cairo_pdf)

cat("\nSaved:\n")
cat("  ", out_pdf, "\n")
cat("  ", out_png, "\n")
cat("  panel_A_landscape.pdf / panel_B_architecture.pdf\n")
cat("  panel_C_TMEM106B_MDD.pdf / panel_D_CACNA1C_BD.pdf / panel_E_RP11-347L18.1_BD.pdf\n")
cat("  panel_A_landscape.tsv / panel_B_locus_disease_events.tsv / coloc_summary_by_cell_type.tsv\n")
cat("\nDone.\n")
