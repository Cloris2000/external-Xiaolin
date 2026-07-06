#!/usr/bin/env Rscript
# Ancestry-matched diagnosis stratum comparison (AD_EUR vs psychiatric_EUR vs EUR ref).

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
    else if (flag == "--output-dir")  { out$output_dir   <- val; i <- i + 2L }
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
    if (length(parts) >= 2) pri[[parts[1]]] <- unique(c(pri[[parts[1]]], paste(parts[-1], collapse = "_")))
  }
  pri
}
is_priority <- function(ct, marker, pri) {
  if (is.null(pri[[ct]])) return(FALSE)
  m <- gsub(":", "_", sub(":.*", "", marker))
  any(sapply(pri[[ct]], function(p) grepl(p, m, fixed = TRUE)))
}

# One lead SNP per (cell type, locus window); diversify across cell types / chromosomes.
select_dumbbell_leads <- function(df, pri, max_per_ct = 2L, max_total = 36L, window_bp = 500000L) {
  interesting <- c(
    "AD_EUR_only", "psychiatric_EUR_only",
    "AD_and_psychiatric_EUR", "EUR_ref_only"
  )

  df %>%
    mutate(
      chr = sub("^chr", "", sub(":.*", "", marker)),
      pos = suppressWarnings(as.numeric(sub(".*:([0-9]+).*", "\\1", marker))),
      locus_key = paste0(chr, ":", floor(pos / window_bp))
    ) %>%
    filter(!is.na(pos)) %>%
    filter(
      (pooled_is_gw == 1L & diagnosis_pattern %in% interesting) |
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
      factor(diagnosis_pattern, levels = c(interesting, "neither_matched_stratum")),
      desc(is_priority), pooled_p
    ) %>%
    slice_head(n = max_total)
}

hits <- read.delim(args$hit_overlap, stringsAsFactors = FALSE)
pri <- load_priority(args$coloc_tsv)

classify <- function(df) {
  ad_sig  <- int1(df$AD_EUR_is_gw) == 1L | int1(df$AD_EUR_is_suggestive) == 1L
  psy_sig <- int1(df$psychiatric_EUR_is_gw) == 1L | int1(df$psychiatric_EUR_is_suggestive) == 1L
  eur_sig <- int1(df$EUR_homogeneous_is_gw) == 1L | int1(df$EUR_homogeneous_is_suggestive) == 1L

  pat <- case_when(
    ad_sig & psy_sig ~ "AD_and_psychiatric_EUR",
    ad_sig & !psy_sig ~ "AD_EUR_only",
    !ad_sig & psy_sig ~ "psychiatric_EUR_only",
    !ad_sig & !psy_sig & eur_sig ~ "EUR_ref_only",
    TRUE ~ "neither_matched_stratum"
  )

  df %>% mutate(
    pooled_p = num(pooled_p), pooled_beta = num(pooled_beta), pooled_i2 = num(pooled_het_isq),
    pooled_is_gw = int1(pooled_is_gw),
    is_priority = as.integer(mapply(is_priority, cell_type, marker, MoreArgs = list(pri = pri))),
    diagnosis_pattern = factor(pat)
  )
}
hits <- classify(hits)

pat_cols <- c(
  AD_and_psychiatric_EUR = "#1a9850",
  AD_EUR_only = "#d73027",
  psychiatric_EUR_only = "#4575b4",
  EUR_ref_only = "#fee08b",
  neither_matched_stratum = "#bdbdbd"
)

save_plot <- function(p, path, w = 10, h = 6) {
  ggsave(path, plot = p, width = w, height = h, dpi = 200, bg = "white")
  message("Saved: ", path)
}

gw <- hits %>% filter(pooled_is_gw == 1L)
pat_counts <- gw %>% count(diagnosis_pattern, name = "n")
p1 <- ggplot(pat_counts, aes(x = reorder(diagnosis_pattern, n), y = n, fill = diagnosis_pattern)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = pat_cols, guide = "none") +
  labs(
    title = "GW hits: ancestry-matched diagnosis replication",
    subtitle = "Pooled GW hits classified by AD_EUR (4 cohorts) vs psychiatric_EUR (CMC x3) significance",
    x = NULL, y = "Count"
  ) +
  theme_classic(base_size = 12)
save_plot(p1, file.path(args$output_dir, "diagnosis_eur_gw_pattern.png"), 10, 5)

top <- select_dumbbell_leads(
  hits,
  pri,
  max_per_ct = args$max_loci_per_cell_type,
  max_total = args$max_total_loci,
  window_bp = args$locus_window_bp
)

if (nrow(top) > 0) {
  write.table(
    top %>% select(cell_type, marker, chr, pos, locus_key, pooled_p, pooled_beta, pooled_i2,
                   diagnosis_pattern, is_priority, AD_EUR_p, psychiatric_EUR_p, EUR_homogeneous_p),
    file.path(args$output_dir, "dumbbell_lead_snps.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  long <- top %>%
    mutate(
      hit_label = paste0(cell_type, " | ", sub("^chr([0-9]+):([0-9]+).*", "chr\\1:\\2", marker)),
      hit_label = factor(hit_label, levels = rev(unique(hit_label)))
    ) %>%
    pivot_longer(
      c(pooled_beta, AD_EUR_beta, psychiatric_EUR_beta, EUR_homogeneous_beta),
      names_to = "stratum", values_to = "beta"
    ) %>%
    mutate(stratum = factor(stratum,
      levels = c("pooled_beta", "EUR_homogeneous_beta", "AD_EUR_beta", "psychiatric_EUR_beta"),
      labels = c("Pooled 15", "EUR 9-cohort", "AD EUR 4", "Psych EUR 3")
    ))

  p2 <- ggplot(long, aes(x = beta, y = hit_label, colour = stratum)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(aes(group = hit_label), colour = "grey75", linewidth = 0.6) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("Pooled 15" = "#333333", "EUR 9-cohort" = "#756bb1",
                                   "AD EUR 4" = "#e41a1c", "Psych EUR 3" = "#2166AC"), name = "Stratum") +
    facet_wrap(~ diagnosis_pattern, scales = "free_y", ncol = 1) +
    labs(
      title = "Lead SNPs per locus: ancestry-matched diagnosis strata",
      subtitle = sprintf(
        "One lead SNP per %dkb locus; up to %d loci per cell type (GW + coloc-priority)",
        args$locus_window_bp %/% 1000L, args$max_loci_per_cell_type
      ),
      x = "Effect size (beta)", y = NULL
    ) %>%
    theme_classic(base_size = 9) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
  save_plot(p2, file.path(args$output_dir, "diagnosis_eur_priority_dumbbell.png"),
            11, max(8, 0.22 * nrow(top)))
}

write.table(
  gw %>% filter(diagnosis_pattern %in% c("AD_and_psychiatric_EUR", "AD_EUR_only", "psychiatric_EUR_only")) %>%
    select(cell_type, marker, pooled_p, pooled_beta, pooled_i2, diagnosis_pattern,
           AD_EUR_p, AD_EUR_beta, psychiatric_EUR_p, psychiatric_EUR_beta, EUR_homogeneous_p),
  file.path(args$output_dir, "gw_diagnosis_eur_classified.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
message("Done.")
