#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--sn_dir", type = "character",
              default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_sn_rosmap_msbb_hbcc_alias15"),
  make_option("--bulk_dir", type = "character",
              default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_13cohorts"),
  make_option("--outdir", type = "character",
              default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity"),
  make_option("--ld_dir", type = "character",
              default = "/external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr"),
  make_option("--gw_thresh", type = "double", default = 5e-8),
  make_option("--sugg_thresh", type = "double", default = 1e-5),
  make_option("--locus_window", type = "integer", default = 1000000),
  make_option("--coloc_window", type = "integer", default = 500000),
  make_option("--max_scatter_points", type = "integer", default = 200000),
  make_option("--prepare_ldsc_inputs", type = "logical", default = TRUE)
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "harmonized_pairs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "coloc_ready_loci"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "ldsc", "sumstats"), recursive = TRUE, showWarnings = FALSE)

crosswalk <- data.frame(
  sn_cell_type = c("ASTRO", "ENDO", "L23IT", "L56NP", "L5ET", "L6B", "L6CT",
                   "LAMP5LHX6", "OLIGO", "OPC", "PVALB", "PVM", "SST", "VIP", "VLMC"),
  bulk_cell_type = c("Astrocyte", "Endothelial", "IT", "L5.6.NP", "L5.ET", "L6b", "L6.CT",
                     "LAMP5", "Oligodendrocyte", "OPC", "PVALB", "Pericyte", "SST", "VIP", "VLMC"),
  stringsAsFactors = FALSE
)

fwrite(crosswalk, file.path(opt$outdir, "cell_type_crosswalk.tsv"), sep = "\t")

message("Validating input METAL tables...")

find_meta_file <- function(dir, cell_type) {
  pattern <- paste0("^", gsub("([.])", "\\\\\\1", cell_type), "_meta_analysis_.*\\.tbl$")
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  files <- files[!grepl("1\\.tbl$", basename(files))]
  if (length(files) == 0) return(NA_character_)
  files[1]
}

crosswalk$sn_meta_tbl <- vapply(crosswalk$sn_cell_type, find_meta_file, character(1), dir = opt$sn_dir)
crosswalk$bulk_meta_tbl <- vapply(crosswalk$bulk_cell_type, find_meta_file, character(1), dir = opt$bulk_dir)
crosswalk$sn_found <- file.exists(crosswalk$sn_meta_tbl)
crosswalk$bulk_found <- file.exists(crosswalk$bulk_meta_tbl)

fwrite(crosswalk, file.path(opt$outdir, "input_file_validation.tsv"), sep = "\t")

if (!all(crosswalk$sn_found & crosswalk$bulk_found)) {
  missing_rows <- crosswalk[!(crosswalk$sn_found & crosswalk$bulk_found), ]
  print(missing_rows)
  stop("Missing one or more input .tbl files; see input_file_validation.tsv")
}

standard_cols <- c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value",
                   "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal")

parse_marker <- function(dt) {
  marker_parts <- tstrsplit(dt$MarkerName, ":", fixed = TRUE)
  dt[, CHR := gsub("^chr", "", marker_parts[[1]], ignore.case = TRUE)]
  dt[, POS := suppressWarnings(as.integer(marker_parts[[2]]))]
  dt[, REF := toupper(marker_parts[[3]])]
  dt[, ALT := toupper(marker_parts[[4]])]
  dt[, CHR_NUM := suppressWarnings(as.integer(CHR))]
  dt
}

read_meta <- function(path) {
  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  keep <- intersect(standard_cols, header)
  dt <- fread(path, select = keep, showProgress = FALSE)
  setnames(dt, old = c("Effect", "StdErr", "P-value"),
           new = c("BETA", "SE", "P"), skip_absent = TRUE)
  dt[, `:=`(
    Allele1 = toupper(Allele1),
    Allele2 = toupper(Allele2),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P)
  )]
  dt <- dt[!is.na(MarkerName) & !is.na(BETA) & !is.na(SE) & !is.na(P) & SE > 0 & P > 0 & P <= 1]
  dt[, Z := BETA / SE]
  parse_marker(dt)
}

is_ambiguous_pair <- function(a1, a2) {
  pair <- paste0(a1, a2)
  pair %in% c("AT", "TA", "CG", "GC")
}

harmonize_pair <- function(sn_dt, bulk_dt) {
  sn <- sn_dt[, .(MarkerName, CHR, CHR_NUM, POS, REF, ALT,
                  A1_sn = Allele1, A2_sn = Allele2, BETA_sn = BETA, SE_sn = SE,
                  P_sn = P, Z_sn = Z, Direction_sn = Direction,
                  HetISq_sn = HetISq, HetPVal_sn = HetPVal)]
  bulk <- bulk_dt[, .(MarkerName,
                      A1_bulk = Allele1, A2_bulk = Allele2, BETA_bulk_raw = BETA,
                      SE_bulk = SE, P_bulk = P, Z_bulk_raw = Z,
                      Direction_bulk = Direction, HetISq_bulk = HetISq, HetPVal_bulk = HetPVal)]
  joined <- merge(sn, bulk, by = "MarkerName", all = FALSE)
  joined[, allele_match := A1_sn == A1_bulk & A2_sn == A2_bulk]
  joined[, allele_swap := A1_sn == A2_bulk & A2_sn == A1_bulk]
  joined[, ambiguous := is_ambiguous_pair(A1_sn, A2_sn)]
  joined[, harmonized_status := fifelse(ambiguous, "ambiguous",
                                 fifelse(allele_match, "match",
                                  fifelse(allele_swap, "swap", "mismatch")))]
  out <- joined[harmonized_status %in% c("match", "swap")]
  out[, `:=`(
    BETA_bulk = fifelse(harmonized_status == "swap", -BETA_bulk_raw, BETA_bulk_raw),
    Z_bulk = fifelse(harmonized_status == "swap", -Z_bulk_raw, Z_bulk_raw),
    sign_concordant = sign(BETA_sn) == sign(fifelse(harmonized_status == "swap", -BETA_bulk_raw, BETA_bulk_raw)),
    logP_sn = -log10(P_sn),
    logP_bulk = -log10(P_bulk)
  )]
  list(joined = joined, harmonized = out)
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

storey_pi1 <- function(p, lambda = 0.5) {
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  if (length(p) == 0) return(NA_real_)
  pi0 <- mean(p > lambda) / (1 - lambda)
  max(0, 1 - min(1, pi0))
}

lead_loci <- function(dt, threshold, window) {
  candidates <- dt[CHR_NUM %between% c(1L, 22L) & P < threshold]
  if (nrow(candidates) == 0) return(candidates[0])
  setorder(candidates, P)
  keep_idx <- logical(nrow(candidates))
  kept <- data.table(CHR_NUM = integer(), POS = integer())
  for (i in seq_len(nrow(candidates))) {
    chr <- candidates$CHR_NUM[i]
    pos <- candidates$POS[i]
    near <- nrow(kept[CHR_NUM == chr & abs(POS - pos) <= window]) > 0
    if (!near) {
      keep_idx[i] <- TRUE
      kept <- rbind(kept, data.table(CHR_NUM = chr, POS = pos))
    }
  }
  candidates[keep_idx]
}

summarize_harmonization <- function(joined, sn_ct, bulk_ct, sn_n, bulk_n) {
  data.table(
    sn_cell_type = sn_ct,
    bulk_cell_type = bulk_ct,
    sn_variants = sn_n,
    bulk_variants = bulk_n,
    shared_marker_variants = nrow(joined),
    retained_harmonized = sum(joined$harmonized_status %in% c("match", "swap")),
    allele_matches = sum(joined$harmonized_status == "match"),
    allele_swaps = sum(joined$harmonized_status == "swap"),
    ambiguous_dropped = sum(joined$harmonized_status == "ambiguous"),
    allele_mismatches_dropped = sum(joined$harmonized_status == "mismatch")
  )
}

top_concordance <- function(h, p_col, n_top) {
  if (nrow(h) == 0) return(NA_real_)
  idx <- order(h[[p_col]])[seq_len(min(n_top, nrow(h)))]
  mean(h$sign_concordant[idx], na.rm = TRUE)
}

metric_summary <- function(h, sn_ct, bulk_ct) {
  data.table(
    sn_cell_type = sn_ct,
    bulk_cell_type = bulk_ct,
    n_harmonized = nrow(h),
    z_cor_pearson = safe_cor(h$Z_sn, h$Z_bulk),
    z_cor_spearman = safe_cor(h$Z_sn, h$Z_bulk, "spearman"),
    beta_cor_pearson = safe_cor(h$BETA_sn, h$BETA_bulk),
    logp_cor_pearson = safe_cor(h$logP_sn, h$logP_bulk),
    sign_concordance_all = mean(h$sign_concordant, na.rm = TRUE),
    sign_concordance_sn_top100 = top_concordance(h, "P_sn", 100),
    sign_concordance_sn_top1000 = top_concordance(h, "P_sn", 1000),
    sign_concordance_bulk_top100 = top_concordance(h, "P_bulk", 100),
    sign_concordance_bulk_top1000 = top_concordance(h, "P_bulk", 1000)
  )
}

replication_rows <- function(h, sn_ct, bulk_ct, source_label, source_p, target_p, source_beta, target_beta) {
  rows <- list()
  thresholds <- c(opt$gw_thresh, opt$sugg_thresh, 1e-4)
  for (threshold in thresholds) {
    selected <- h[source_p < threshold]
    rows[[length(rows) + 1]] <- data.table(
      sn_cell_type = sn_ct,
      bulk_cell_type = bulk_ct,
      direction = source_label,
      selection = paste0("p<", threshold),
      n_selected = nrow(selected),
      median_target_logp = if (nrow(selected)) median(-log10(target_p[source_p < threshold]), na.rm = TRUE) else NA_real_,
      mean_target_logp = if (nrow(selected)) mean(-log10(target_p[source_p < threshold]), na.rm = TRUE) else NA_real_,
      target_nominal_fraction = if (nrow(selected)) mean(target_p[source_p < threshold] < 0.05, na.rm = TRUE) else NA_real_,
      target_fdr05_fraction = if (nrow(selected)) mean(p.adjust(target_p[source_p < threshold], "BH") < 0.05, na.rm = TRUE) else NA_real_,
      same_direction_fraction = if (nrow(selected)) mean(sign(source_beta[source_p < threshold]) == sign(target_beta[source_p < threshold]), na.rm = TRUE) else NA_real_,
      pi1_target = if (nrow(selected)) storey_pi1(target_p[source_p < threshold]) else NA_real_
    )
  }
  for (n_top in c(100, 1000)) {
    if (nrow(h) == 0) {
      selected_idx <- integer()
    } else {
      selected_idx <- order(source_p)[seq_len(min(n_top, nrow(h)))]
    }
    selected_target <- target_p[selected_idx]
    rows[[length(rows) + 1]] <- data.table(
      sn_cell_type = sn_ct,
      bulk_cell_type = bulk_ct,
      direction = source_label,
      selection = paste0("top", n_top),
      n_selected = length(selected_idx),
      median_target_logp = if (length(selected_idx)) median(-log10(selected_target), na.rm = TRUE) else NA_real_,
      mean_target_logp = if (length(selected_idx)) mean(-log10(selected_target), na.rm = TRUE) else NA_real_,
      target_nominal_fraction = if (length(selected_idx)) mean(selected_target < 0.05, na.rm = TRUE) else NA_real_,
      target_fdr05_fraction = if (length(selected_idx)) mean(p.adjust(selected_target, "BH") < 0.05, na.rm = TRUE) else NA_real_,
      same_direction_fraction = if (length(selected_idx)) mean(sign(source_beta[selected_idx]) == sign(target_beta[selected_idx]), na.rm = TRUE) else NA_real_,
      pi1_target = if (length(selected_idx)) storey_pi1(selected_target) else NA_real_
    )
  }
  rbindlist(rows, fill = TRUE)
}

lead_concordance_rows <- function(source_dt, target_dt, h, sn_ct, bulk_ct, direction, threshold) {
  leads <- lead_loci(source_dt, threshold, opt$locus_window)
  if (nrow(leads) == 0) return(data.table())
  target_by_marker <- target_dt[, .(MarkerName, P, BETA, SE, CHR_NUM, POS)]
  h_by_marker <- h[, .(MarkerName, BETA_sn, BETA_bulk, P_sn, P_bulk, sign_concordant)]
  rows <- lapply(seq_len(nrow(leads)), function(i) {
    lead <- leads[i]
    region <- target_by_marker[CHR_NUM == lead$CHR_NUM & abs(POS - lead$POS) <= opt$locus_window]
    same <- target_by_marker[MarkerName == lead$MarkerName]
    hz <- h_by_marker[MarkerName == lead$MarkerName]
    data.table(
      sn_cell_type = sn_ct,
      bulk_cell_type = bulk_ct,
      direction = direction,
      threshold = threshold,
      source_marker = lead$MarkerName,
      chr = lead$CHR,
      pos = lead$POS,
      source_p = lead$P,
      source_beta = lead$BETA,
      target_same_marker_p = if (nrow(same)) same$P[1] else NA_real_,
      target_same_marker_beta = if (nrow(same)) same$BETA[1] else NA_real_,
      same_marker_harmonized = nrow(hz) > 0,
      same_marker_direction_concordant = if (nrow(hz)) hz$sign_concordant[1] else NA,
      target_region_min_p = if (nrow(region)) min(region$P, na.rm = TRUE) else NA_real_,
      target_region_best_marker = if (nrow(region)) region$MarkerName[which.min(region$P)] else NA_character_,
      target_region_passes_threshold = if (nrow(region)) any(region$P < threshold, na.rm = TRUE) else FALSE
    )
  })
  rbindlist(rows, fill = TRUE)
}

prepare_rsid_lookup <- function(ld_dir) {
  message("Building HapMap3 rsID lookup from ", ld_dir)
  pieces <- lapply(1:22, function(chr) {
    f <- file.path(ld_dir, paste0(chr, ".l2.ldscore.gz"))
    if (!file.exists(f)) {
      warning("Missing LD score file: ", f)
      return(NULL)
    }
    fread(f, select = c("CHR", "SNP", "BP"), showProgress = FALSE)[
      , .(CHR_NUM = as.integer(CHR), POS = as.integer(BP), SNP)]
  })
  rbindlist(pieces, fill = TRUE)
}

variant_n_from_harmonized <- function(sn_ct) {
  files <- list.files(file.path(opt$sn_dir, "harmonized", sn_ct, "harmonized"),
                      pattern = "_harmonized\\.raw_p$", full.names = TRUE)
  if (length(files) == 0) return(data.table(MarkerName = character(), N = numeric()))
  pieces <- lapply(files, function(f) {
    fread(f, select = c("ID", "N"), showProgress = FALSE)[
      !is.na(ID) & !is.na(N), .(MarkerName = ID, N = as.numeric(N))]
  })
  rbindlist(pieces)[, .(N = sum(N, na.rm = TRUE)), by = MarkerName]
}

write_sn_ldsc_input <- function(sn_dt, sn_ct, rsid_lookup) {
  n_by_variant <- variant_n_from_harmonized(sn_ct)
  out <- merge(sn_dt, n_by_variant, by = "MarkerName", all = FALSE)
  out <- merge(out, rsid_lookup, by = c("CHR_NUM", "POS"), all = FALSE, allow.cartesian = TRUE)
  out <- out[!is.na(SNP) & SNP != "" & !duplicated(SNP)]
  out <- out[!is.na(N) & N > 0 & !is.na(P) & P > 0 & P <= 1 & !is.na(BETA)]
  out <- out[, .(SNP, A1 = Allele1, A2 = Allele2, P, BETA, N)]
  out_path <- file.path(opt$outdir, "ldsc", "sumstats", paste0("sn_", sn_ct, ".ldsc_input.tsv"))
  fwrite(out, out_path, sep = "\t")
  data.table(sn_cell_type = sn_ct, sn_ldsc_input = out_path, sn_ldsc_rows = nrow(out))
}

harm_summaries <- list()
metric_summaries <- list()
replication_summaries <- list()
lead_summaries <- list()
coloc_manifest <- list()
ldsc_manifest <- list()
plot_samples <- list()

rsid_lookup <- NULL
if (isTRUE(opt$prepare_ldsc_inputs)) {
  rsid_lookup <- prepare_rsid_lookup(opt$ld_dir)
}

for (i in seq_len(nrow(crosswalk))) {
  sn_ct <- crosswalk$sn_cell_type[i]
  bulk_ct <- crosswalk$bulk_cell_type[i]
  message("Processing ", sn_ct, " vs ", bulk_ct)

  sn_dt <- read_meta(crosswalk$sn_meta_tbl[i])
  bulk_dt <- read_meta(crosswalk$bulk_meta_tbl[i])

  hz <- harmonize_pair(sn_dt, bulk_dt)
  h <- hz$harmonized
  joined <- hz$joined

  harm_summaries[[sn_ct]] <- summarize_harmonization(joined, sn_ct, bulk_ct, nrow(sn_dt), nrow(bulk_dt))
  metric_summaries[[sn_ct]] <- metric_summary(h, sn_ct, bulk_ct)

  h_path <- file.path(opt$outdir, "harmonized_pairs", paste0(sn_ct, "_vs_", bulk_ct, ".tsv.gz"))
  fwrite(h, h_path, sep = "\t")

  if (nrow(h) > 0) {
    if (nrow(h) > opt$max_scatter_points) {
      set.seed(42 + i)
      sample_idx <- sample.int(nrow(h), opt$max_scatter_points)
      plot_samples[[sn_ct]] <- h[sample_idx, .(sn_cell_type = sn_ct, bulk_cell_type = bulk_ct, Z_sn, Z_bulk, logP_sn, logP_bulk)]
    } else {
      plot_samples[[sn_ct]] <- h[, .(sn_cell_type = sn_ct, bulk_cell_type = bulk_ct, Z_sn, Z_bulk, logP_sn, logP_bulk)]
    }
  }

  replication_summaries[[paste0(sn_ct, "_sn_to_bulk")]] <-
    replication_rows(h, sn_ct, bulk_ct, "sn_to_bulk", h$P_sn, h$P_bulk, h$BETA_sn, h$BETA_bulk)
  replication_summaries[[paste0(sn_ct, "_bulk_to_sn")]] <-
    replication_rows(h, sn_ct, bulk_ct, "bulk_to_sn", h$P_bulk, h$P_sn, h$BETA_bulk, h$BETA_sn)

  for (threshold in c(opt$gw_thresh, opt$sugg_thresh)) {
    lead_summaries[[paste(sn_ct, "sn", threshold)]] <-
      lead_concordance_rows(sn_dt, bulk_dt, h, sn_ct, bulk_ct, "sn_to_bulk", threshold)
    lead_summaries[[paste(sn_ct, "bulk", threshold)]] <-
      lead_concordance_rows(bulk_dt, sn_dt, h, sn_ct, bulk_ct, "bulk_to_sn", threshold)
  }

  coloc_leads <- unique(rbindlist(list(
    lead_loci(sn_dt, opt$sugg_thresh, opt$locus_window)[, .(MarkerName, CHR_NUM, CHR, POS, P, source = "sn")],
    lead_loci(bulk_dt, opt$sugg_thresh, opt$locus_window)[, .(MarkerName, CHR_NUM, CHR, POS, P, source = "bulk")]
  ), fill = TRUE))

  if (nrow(coloc_leads) > 0 && nrow(h) > 0) {
    for (j in seq_len(nrow(coloc_leads))) {
      lead <- coloc_leads[j]
      region <- h[CHR_NUM == lead$CHR_NUM & POS >= lead$POS - opt$coloc_window & POS <= lead$POS + opt$coloc_window]
      if (nrow(region) == 0) next
      start <- max(0L, lead$POS - opt$coloc_window)
      end <- lead$POS + opt$coloc_window
      locus_id <- paste(sn_ct, "vs", bulk_ct, paste0("chr", lead$CHR_NUM), start, end, sep = "_")
      locus_path <- file.path(opt$outdir, "coloc_ready_loci", paste0(locus_id, ".tsv.gz"))
      fwrite(region[, .(MarkerName, CHR, POS, REF, ALT, A1 = A1_sn, A2 = A2_sn,
                        BETA_sn, SE_sn, P_sn, BETA_bulk, SE_bulk, P_bulk,
                        Z_sn, Z_bulk, sign_concordant)], locus_path, sep = "\t")
      coloc_manifest[[length(coloc_manifest) + 1]] <- data.table(
        sn_cell_type = sn_ct,
        bulk_cell_type = bulk_ct,
        locus_id = locus_id,
        lead_source = lead$source,
        lead_marker = lead$MarkerName,
        chr = lead$CHR,
        lead_pos = lead$POS,
        lead_p = lead$P,
        window_start = start,
        window_end = end,
        n_harmonized_snps = nrow(region),
        file = locus_path
      )
    }
  }

  if (isTRUE(opt$prepare_ldsc_inputs)) {
    sn_ldsc <- write_sn_ldsc_input(sn_dt, sn_ct, rsid_lookup)
    bulk_ldsc_input <- file.path(opt$bulk_dir, "ldsc", "sumstats", paste0(bulk_ct, ".ldsc_input.tsv"))
    ldsc_manifest[[sn_ct]] <- cbind(
      data.table(sn_cell_type = sn_ct, bulk_cell_type = bulk_ct,
                 bulk_ldsc_input = bulk_ldsc_input,
                 bulk_ldsc_input_found = file.exists(bulk_ldsc_input)),
      sn_ldsc[, .(sn_ldsc_input, sn_ldsc_rows)]
    )
  }
}

harm_summary <- rbindlist(harm_summaries, fill = TRUE)
metric_summary_all <- rbindlist(metric_summaries, fill = TRUE)
replication_summary <- rbindlist(replication_summaries, fill = TRUE)
lead_summary <- rbindlist(lead_summaries, fill = TRUE)
coloc_manifest_dt <- rbindlist(coloc_manifest, fill = TRUE)
ldsc_manifest_dt <- rbindlist(ldsc_manifest, fill = TRUE)

fwrite(harm_summary, file.path(opt$outdir, "harmonization_summary.tsv"), sep = "\t")
fwrite(metric_summary_all, file.path(opt$outdir, "snp_similarity_summary.tsv"), sep = "\t")
fwrite(replication_summary, file.path(opt$outdir, "replication_enrichment_summary.tsv"), sep = "\t")
fwrite(lead_summary, file.path(opt$outdir, "lead_locus_concordance.tsv"), sep = "\t")
fwrite(coloc_manifest_dt, file.path(opt$outdir, "coloc_locus_manifest.tsv"), sep = "\t")
if (nrow(ldsc_manifest_dt) > 0) {
  fwrite(ldsc_manifest_dt, file.path(opt$outdir, "ldsc", "ldsc_input_manifest.tsv"), sep = "\t")
}

plot_theme <- theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())

p_z <- ggplot(metric_summary_all, aes(x = reorder(sn_cell_type, z_cor_pearson), y = z_cor_pearson)) +
  geom_col(fill = "#2171b5") +
  coord_flip() +
  labs(title = "Genome-wide Z-score correlation", x = NULL, y = "Pearson r") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$outdir, "plots", "z_correlation_by_cell_type.png"), p_z, width = 7, height = 5, dpi = 300)

p_sign <- ggplot(metric_summary_all, aes(x = reorder(sn_cell_type, sign_concordance_all), y = sign_concordance_all)) +
  geom_col(fill = "#238b45") +
  coord_flip() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = "Genome-wide effect direction concordance", x = NULL, y = "Fraction same sign") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$outdir, "plots", "direction_concordance_by_cell_type.png"), p_sign, width = 7, height = 5, dpi = 300)

if (length(plot_samples) > 0) {
  scatter_dt <- rbindlist(plot_samples, fill = TRUE)
  p_scatter <- ggplot(scatter_dt, aes(x = Z_sn, y = Z_bulk)) +
    geom_point(alpha = 0.08, size = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    facet_wrap(~sn_cell_type, scales = "free", ncol = 5) +
    labs(title = "sn vs bulk SNP Z-scores", x = "sn Z", y = "bulk Z") +
    theme_bw(base_size = 9)
  ggsave(file.path(opt$outdir, "plots", "z_score_scatter_facets.png"), p_scatter, width = 14, height = 9, dpi = 300)
}

if (nrow(replication_summary) > 0) {
  p_rep <- ggplot(replication_summary[selection %in% c("p<1e-05", "top1000")],
                  aes(x = sn_cell_type, y = target_nominal_fraction, fill = direction)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    facet_wrap(~selection, ncol = 1) +
    labs(title = "Replication-style nominal enrichment", x = NULL,
         y = "Fraction target P < 0.05", fill = NULL) +
    plot_theme
  ggsave(file.path(opt$outdir, "plots", "replication_nominal_enrichment.png"), p_rep, width = 11, height = 7, dpi = 300)
}

if (nrow(lead_summary) > 0) {
  locus_counts <- lead_summary[, .(
    n_leads = .N,
    n_same_marker = sum(!is.na(target_same_marker_p)),
    n_region_pass = sum(target_region_passes_threshold, na.rm = TRUE),
    same_marker_direction_concordance = mean(same_marker_direction_concordant, na.rm = TRUE)
  ), by = .(sn_cell_type, bulk_cell_type, direction, threshold)]
  fwrite(locus_counts, file.path(opt$outdir, "lead_locus_concordance_summary.tsv"), sep = "\t")
  p_locus <- ggplot(locus_counts, aes(x = sn_cell_type, y = n_region_pass, fill = direction)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    facet_wrap(~threshold, ncol = 1, scales = "free_y") +
    labs(title = "Lead loci with target signal in the same region", x = NULL,
         y = "Lead loci with target region P below threshold", fill = NULL) +
    plot_theme
  ggsave(file.path(opt$outdir, "plots", "lead_locus_region_concordance.png"), p_locus, width = 11, height = 7, dpi = 300)
}

message("Done. Outputs written to: ", opt$outdir)
