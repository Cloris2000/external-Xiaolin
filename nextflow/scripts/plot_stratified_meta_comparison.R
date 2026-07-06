#!/usr/bin/env Rscript
# Visualise pooled vs stratified meta comparison (AD / EUR / HBCC).
#
# Usage:
#   Rscript plot_stratified_meta_comparison.R \
#     --hit-overlap  results/meta_sensitivity/stratified_comparison_pac/stratified_hit_overlap.tsv \
#     --coloc-tsv    results/coloc/coloc_all_results.tsv \
#     --output-dir   results/meta_sensitivity/stratified_comparison_pac/figures

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args_raw <- commandArgs(trailingOnly = TRUE)
parse_args <- function(argv) {
  out <- list(
    hit_overlap = NULL, coloc_tsv = NULL, output_dir = "figures",
    max_loci_per_cell_type = 2L, max_total_loci = 36L, locus_window_bp = 500000L
  )
  i <- 1L
  while (i <= length(argv)) {
    flag <- argv[i]
    val  <- if (i + 1L <= length(argv)) argv[i + 1L] else NA_character_
    if      (flag == "--hit-overlap")  { out$hit_overlap  <- val; i <- i + 2L }
    else if (flag == "--coloc-tsv")    { out$coloc_tsv    <- val; i <- i + 2L }
    else if (flag == "--output-dir")   { out$output_dir   <- val; i <- i + 2L }
    else if (flag == "--max-loci-per-cell-type") { out$max_loci_per_cell_type <- as.integer(val); i <- i + 2L }
    else if (flag == "--max-total-loci")         { out$max_total_loci <- as.integer(val); i <- i + 2L }
    else if (flag == "--locus-window-bp")        { out$locus_window_bp <- as.integer(val); i <- i + 2L }
    else { i <- i + 1L }
  }
  out
}
args <- parse_args(args_raw)
if (is.null(args$hit_overlap)) stop("Required: --hit-overlap")
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
int1 <- function(x) as.integer(num(x) %||% 0)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

load_priority <- function(path) {
  pri <- list()
  if (is.null(path) || !file.exists(path)) return(pri)
  d <- read.delim(path, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(d))) {
    parts <- strsplit(d$locus_id[i], "_", fixed = FALSE)[[1]]
    if (length(parts) >= 2) {
      ct <- parts[1]
      loc <- paste(parts[-1], collapse = "_")
      pri[[ct]] <- unique(c(pri[[ct]], loc))
    }
  }
  pri
}

is_priority <- function(ct, marker, pri) {
  if (is.null(pri[[ct]])) return(FALSE)
  m <- sub(":.*", "", marker)
  m <- gsub(":", "_", m)
  any(sapply(pri[[ct]], function(p) grepl(p, m, fixed = TRUE)))
}

select_dumbbell_leads <- function(df, pri, max_per_ct = 2L, max_total = 36L, window_bp = 500000L) {
  interesting <- c(
    "all_three_strata", "AD_and_EUR_only", "AD_only", "EUR_only", "HBCC_only",
    "AD_and_HBCC_only", "EUR_and_HBCC_only"
  )

  df %>%
    mutate(
      chr = sub("^chr", "", sub(":.*", "", marker)),
      pos = suppressWarnings(as.numeric(sub(".*:([0-9]+).*", "\\1", marker))),
      locus_key = paste0(chr, ":", floor(pos / window_bp))
    ) %>%
    filter(!is.na(pos)) %>%
    filter(
      (pooled_is_gw == 1L & replication_pattern %in% interesting) |
        (is_priority == 1L & pooled_is_gw == 1L)
    ) %>%
    group_by(cell_type, locus_key) %>%
    slice_min(order_by = pooled_p, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(cell_type) %>%
    arrange(desc(is_priority), pooled_p) %>%
    filter(row_number() <= max_per_ct) %>%
    ungroup() %>%
    arrange(
      factor(replication_pattern, levels = c(interesting, "neither_stratum_sig")),
      desc(is_priority), pooled_p
    ) %>%
    slice_head(n = max_total)
}

hits <- read.delim(args$hit_overlap, stringsAsFactors = FALSE)
pri <- load_priority(args$coloc_tsv)

strata <- c("AD_neurological", "EUR_homogeneous", "HBCC_AFR_enriched")

classify <- function(df) {
  ad_sig  <- int1(df$AD_neurological_is_gw) == 1L | int1(df$AD_neurological_is_suggestive) == 1L
  eur_sig <- int1(df$EUR_homogeneous_is_gw) == 1L | int1(df$EUR_homogeneous_is_suggestive) == 1L
  hb_sig  <- int1(df$HBCC_AFR_enriched_is_gw) == 1L | int1(df$HBCC_AFR_enriched_is_suggestive) == 1L
  gw <- int1(df$pooled_is_gw) == 1L

  pat <- case_when(
    ad_sig & eur_sig & hb_sig ~ "all_three_strata",
    ad_sig & eur_sig & !hb_sig ~ "AD_and_EUR_only",
    ad_sig & !eur_sig & !hb_sig ~ "AD_only",
    !ad_sig & eur_sig & !hb_sig ~ "EUR_only",
    !ad_sig & !eur_sig & hb_sig ~ "HBCC_only",
    ad_sig & !eur_sig & hb_sig ~ "AD_and_HBCC_only",
    !ad_sig & eur_sig & hb_sig ~ "EUR_and_HBCC_only",
    TRUE ~ "neither_stratum_sig"
  )

  df %>%
    mutate(
      pooled_p = num(pooled_p),
      pooled_beta = num(pooled_beta),
      pooled_i2 = num(pooled_het_isq),
      pooled_is_gw = int1(pooled_is_gw),
      is_priority = as.integer(mapply(is_priority, cell_type, marker, MoreArgs = list(pri = pri))),
      replication_pattern = factor(pat),
      sig_level = ifelse(pooled_is_gw == 1L, "Genome-wide", "Suggestive")
    )
}

hits <- classify(hits)

pat_cols <- c(
  all_three_strata = "#1a9850",
  AD_and_EUR_only = "#4575b4",
  AD_only = "#74add1",
  EUR_only = "#313695",
  HBCC_only = "#f46d43",
  AD_and_HBCC_only = "#d73027",
  EUR_and_HBCC_only = "#762a83",
  neither_stratum_sig = "#bdbdbd"
)

save_plot <- function(p, path, w = 10, h = 7) {
  ggsave(path, plot = p, width = w, height = h, dpi = 200, bg = "white")
  message("Saved: ", path)
}

# 1. GW pattern summary
gw <- hits %>% filter(pooled_is_gw == 1L)
pat_counts <- gw %>% count(replication_pattern, name = "n")
p1 <- ggplot(pat_counts, aes(x = reorder(replication_pattern, n), y = n, fill = replication_pattern)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = pat_cols, guide = "none") +
  labs(
    title = "Genome-wide hits: replication across AD / EUR / HBCC strata",
    subtitle = "Pooled 15-cohort GW hits classified by which stratum metas are significant",
    x = NULL, y = "Number of GW hits"
  ) +
  theme_classic(base_size = 12)

save_plot(p1, file.path(args$output_dir, "gw_replication_pattern.png"), 10, 6)

# 2. Priority loci dumbbell (one lead SNP per locus; diverse cell types)
top <- select_dumbbell_leads(
  hits, pri,
  max_per_ct = args$max_loci_per_cell_type,
  max_total = args$max_total_loci,
  window_bp = args$locus_window_bp
)

if (nrow(top) > 0) {
  write.table(
    top %>% select(cell_type, marker, chr, pos, locus_key, pooled_p, pooled_beta, pooled_i2,
                   replication_pattern, is_priority),
    file.path(args$output_dir, "dumbbell_lead_snps.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  long <- top %>%
    mutate(
      hit_label = paste0(cell_type, " | ", sub("^chr([0-9]+):([0-9]+).*", "chr\\1:\\2", marker)),
      hit_label = factor(hit_label, levels = rev(unique(hit_label)))
    ) %>%
    pivot_longer(
      c(pooled_beta, AD_neurological_beta, EUR_homogeneous_beta, HBCC_AFR_enriched_beta),
      names_to = "stratum", values_to = "beta"
    ) %>%
    mutate(stratum = factor(stratum,
      levels = c("pooled_beta", "AD_neurological_beta", "EUR_homogeneous_beta", "HBCC_AFR_enriched_beta"),
      labels = c("Pooled", "AD meta", "EUR meta", "HBCC meta")
    ))

  p2 <- ggplot(long, aes(x = beta, y = hit_label, colour = stratum)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(aes(group = hit_label), colour = "grey75", linewidth = 0.6) +
    geom_point(size = 3) +
    scale_colour_manual(values = c(Pooled = "#333333", `AD meta` = "#e41a1c",
                                   `EUR meta` = "#2166AC", `HBCC meta` = "#D95F02"), name = "Stratum") +
    facet_wrap(~ replication_pattern, scales = "free_y", ncol = 1) +
    labs(
      title = "Lead SNPs per locus: AD / EUR / HBCC strata",
      subtitle = sprintf(
        "One lead SNP per %dkb locus; up to %d loci per cell type (GW + coloc-priority)",
        args$locus_window_bp %/% 1000L, args$max_loci_per_cell_type
      ),
      x = "Effect size (beta)", y = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

  save_plot(p2, file.path(args$output_dir, "priority_loci_dumbbell.png"),
            11, max(8, 0.22 * nrow(top)))
}

# 3. AD vs pooled concordance for GW
ad_gw <- gw %>%
  mutate(
    ad_sig = int1(AD_neurological_is_gw) == 1L,
    ad_beta = num(AD_neurological_beta),
    pattern = case_when(
      ad_sig ~ "GW in pooled + AD",
      int1(AD_neurological_is_suggestive) == 1L ~ "GW pooled, suggestive AD",
      TRUE ~ "GW pooled, absent AD"
    )
  )

ad_counts <- ad_gw %>% count(pattern, name = "n")
p3 <- ggplot(ad_counts, aes(x = pattern, y = n, fill = pattern)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.3, size = 4) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title = "AD-stratified meta replication of pooled GW hits",
    x = NULL, y = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(p3, file.path(args$output_dir, "ad_replication_summary.png"), 9, 5)

# 4. Export tables
write.table(
  gw %>% filter(replication_pattern %in% c("all_three_strata", "AD_and_EUR_only", "AD_only")) %>%
    select(cell_type, marker, pooled_p, pooled_beta, pooled_i2, replication_pattern,
           AD_neurological_p, AD_neurological_beta, EUR_homogeneous_p, HBCC_AFR_enriched_p),
  file.path(args$output_dir, "gw_hits_ad_replicated.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
)

message("Done.")
