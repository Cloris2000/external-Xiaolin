#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--cohort"), type = "character", default = NULL,
              help = "Cohort to benchmark: CMC_MSSM or ROSMAP"),
  make_option(c("--results_root"), type = "character", default = NULL,
              help = "Directory containing arm subdirectories (default: nextflow benchmark results root for cohort)"),
  make_option(c("--arms"), type = "character", default = "none,shared,full",
              help = "Comma-delimited arm names to compare (default: none,shared,full)"),
  make_option(c("--proportions_file"), type = "character", default = "cell_proportions_scaled.csv",
              help = "Per-arm MGP proportions filename (default: cell_proportions_scaled.csv)"),
  make_option(c("--metadata_file"), type = "character", default = "metadata_cleaned.csv",
              help = "Per-arm metadata filename (default: metadata_cleaned.csv)"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory (default: <results_root>/benchmark_summary)"),
  make_option(c("--ground_truth_cmc"), type = "character", default = NULL,
              help = "Path to CMC snRNA ground truth CSV (optional; default in script)"),
  make_option(c("--ground_truth_rosmap"), type = "character", default = NULL,
              help = "Path to ROSMAP snRNA ground truth CSV (optional; default in script)")
)

opt <- parse_args(OptionParser(option_list = option_list))

`%||%` <- function(a, b) if (is.null(a)) b else a

if (is.null(opt$cohort)) {
  stop("ERROR: --cohort is required")
}

if (!opt$cohort %in% c("CMC_MSSM", "ROSMAP")) {
  stop("ERROR: --cohort must be one of: CMC_MSSM, ROSMAP")
}

if (is.null(opt$results_root)) {
  opt$results_root <- file.path(
    "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/benchmark_tech_cov",
    opt$cohort
  )
}

if (is.null(opt$output_dir)) {
  opt$output_dir <- file.path(opt$results_root, "benchmark_summary")
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

arm_names <- trimws(unlist(strsplit(opt$arms, ",")))
arm_names <- arm_names[arm_names != ""]

normalize_cell_type <- function(x) {
  normalized <- tolower(x)
  normalized <- gsub("^exc[ _./-]*", "", normalized)
  normalized <- gsub("^inh[ _./-]*", "", normalized)
  normalized <- gsub("[ _./-]+", "", normalized)
  
  synonym_map <- c(
    "astro" = "astrocyte",
    "astrocyte" = "astrocyte",
    "endo" = "endothelial",
    "endothelial" = "endothelial",
    "micro" = "microglia",
    "microglia" = "microglia",
    "oligo" = "oligodendrocyte",
    "oligodendrocyte" = "oligodendrocyte",
    "it" = "it",
    "l4it" = "l4it",
    "l5et" = "l5et",
    "l56itcar3" = "l56itcar3",
    "l6itcar3" = "l56itcar3",
    "l56np" = "l56np",
    "l6ct" = "l6ct",
    "l6it" = "l6it",
    "l6b" = "l6b",
    "lamp5" = "lamp5",
    "lamp5lhx6" = "lamp5lhx6",
    "sncg" = "sncg",
    "sst" = "sst",
    "sstchodl" = "sstchodl",
    "chandelier" = "chandelier",
    "pax6" = "pax6",
    "pericyte" = "pericyte",
    "pvalb" = "pvalb",
    "vip" = "vip",
    "vlmc" = "vlmc",
    "opc" = "opc",
    "l23it" = "l23it"
  )
  
  ifelse(normalized %in% names(synonym_map), synonym_map[normalized], normalized)
}

find_id_column <- function(df) {
  candidate_cols <- c("specimenID", "synapseID", "sampleID", "id", "individualID", "Sample")
  matched_col <- candidate_cols[candidate_cols %in% colnames(df)][1]
  if (is.na(matched_col) || is.null(matched_col)) {
    matched_col <- colnames(df)[1]
  }
  matched_col
}

find_best_metadata_match <- function(metadata_df, sample_ids) {
  candidate_cols <- c("specimenID", "synapseID", "id", "individualID", "sampleID", "standardID")
  candidate_cols <- candidate_cols[candidate_cols %in% colnames(metadata_df)]
  
  if (length(candidate_cols) == 0) {
    stop("ERROR: Could not find a metadata sample ID column")
  }
  
  match_counts <- sapply(candidate_cols, function(col) {
    sum(metadata_df[[col]] %in% sample_ids, na.rm = TRUE)
  })
  
  candidate_cols[which.max(match_counts)]
}

load_ground_truth_cmc <- function(path = NULL) {
  path <- path %||% "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/PsychEncode_label_transferred_snCTP.csv"
  sn_ctp <- read_csv(path, show_col_types = FALSE)
  sn_ctp <- sn_ctp %>% filter(Institution == "MtSinai")
  sn_ctp$benchmark_id <- sn_ctp$subject_id
  cell_type_cols <- setdiff(colnames(sn_ctp), c("", "X1", "ID", "Institution", "subject_id", "benchmark_id"))
  sn_ctp %>%
    select(benchmark_id, all_of(cell_type_cols))
}

load_ground_truth_rosmap <- function(metadata_df, path = NULL) {
  # Ensure we have the required columns for ground truth matching
  required_cols <- c("specimenID", "projid")
  if (!all(required_cols %in% colnames(metadata_df))) {
    missing_cols <- setdiff(required_cols, colnames(metadata_df))
    stop("ERROR: ROSMAP metadata must include ", paste(missing_cols, collapse = ", "), " for ground-truth matching")
  }
  path <- path %||% "/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv"
  sn_raw <- read_csv(path, show_col_types = FALSE)
  
  matcher <- metadata_df %>%
    select(specimenID, projid) %>%
    distinct(specimenID, .keep_all = TRUE) %>%
    distinct(projid, .keep_all = TRUE)
  
  sn_ctp <- sn_raw %>%
    left_join(matcher, by = c("ID" = "projid")) %>%
    filter(!is.na(specimenID)) %>%
    rename(benchmark_id = specimenID)

  # Keep raw ROSMAP snRNA proportions; downstream matching uses overlap only.
  cell_type_cols <- setdiff(colnames(sn_ctp), c("", "X1", "ID", "benchmark_id"))
  sn_ctp %>%
    select(benchmark_id, all_of(cell_type_cols))
}

prepare_estimates <- function(arm_dir, proportions_file, metadata_file) {
  proportions_path <- file.path(arm_dir, proportions_file)
  metadata_path <- file.path(arm_dir, metadata_file)
  
  # Try to load combined_metrics.csv for ID matching (has specimenID and projid),
  # fall back to metadata_cleaned.csv if not available
  combined_metrics_path <- file.path(arm_dir, "combined_metrics.csv")
  metadata_for_id_matching <- NULL
  
  if (file.exists(combined_metrics_path)) {
    metadata_for_id_matching <- read.csv(combined_metrics_path, check.names = FALSE)
    cat("  Using combined_metrics.csv for ID matching\n")
  } else if (file.exists(metadata_path)) {
    metadata_for_id_matching <- read.csv(metadata_path, check.names = FALSE)
    cat("  Using", metadata_file, "for ID matching\n")
  } else {
    stop("ERROR: Missing both combined_metrics.csv and", metadata_file, "files in", arm_dir)
  }
  
  if (!file.exists(proportions_path)) {
    stop("ERROR: Missing proportions file: ", proportions_path)
  }
  
  estimates_df <- read.csv(proportions_path, check.names = FALSE)
  
  estimate_id_col <- find_id_column(estimates_df)
  metadata_match_col <- find_best_metadata_match(metadata_for_id_matching, estimates_df[[estimate_id_col]])
  
  mapping_df <- metadata_for_id_matching %>%
    mutate(match_id = .data[[metadata_match_col]]) %>%
    mutate(benchmark_id = if ("specimenID" %in% colnames(metadata_for_id_matching)) specimenID else match_id) %>%
    select(match_id, benchmark_id) %>%
    distinct(match_id, .keep_all = TRUE)
  
  estimates_df <- estimates_df %>%
    mutate(match_id = .data[[estimate_id_col]]) %>%
    left_join(mapping_df, by = "match_id") %>%
    filter(!is.na(benchmark_id))
  
  cell_type_cols <- setdiff(colnames(estimates_df), c(estimate_id_col, "match_id", "benchmark_id"))
  
  estimates_agg <- estimates_df %>%
    group_by(benchmark_id) %>%
    summarize(across(all_of(cell_type_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  list(
    estimates = estimates_agg,
    metadata = metadata_for_id_matching
  )
}

match_cell_types <- function(estimates_df, ground_truth_df) {
  estimate_cols <- setdiff(colnames(estimates_df), "benchmark_id")
  ground_truth_cols <- setdiff(colnames(ground_truth_df), "benchmark_id")
  
  estimate_norm <- setNames(vapply(estimate_cols, normalize_cell_type, character(1)), estimate_cols)
  ground_truth_norm <- setNames(vapply(ground_truth_cols, normalize_cell_type, character(1)), ground_truth_cols)
  shared_norm <- intersect(unname(estimate_norm), unname(ground_truth_norm))
  
  mapping_df <- lapply(shared_norm, function(norm_name) {
    data.frame(
      canonical_cell_type = norm_name,
      estimate_cell_type = names(estimate_norm)[match(norm_name, estimate_norm)],
      ground_truth_cell_type = names(ground_truth_norm)[match(norm_name, ground_truth_norm)],
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  mapping_df
}

compute_correlations <- function(estimates_df, ground_truth_df, cell_type_mapping, arm_name) {
  merged_df <- inner_join(estimates_df, ground_truth_df, by = "benchmark_id", suffix = c(".est", ".gt"))
  
  if (nrow(merged_df) == 0) {
    stop("ERROR: No samples matched between estimates and ground truth for arm ", arm_name)
  }
  
  lapply(seq_len(nrow(cell_type_mapping)), function(i) {
    estimate_col <- cell_type_mapping$estimate_cell_type[i]
    ground_truth_col <- cell_type_mapping$ground_truth_cell_type[i]
    
    # Find the actual column names in merged_df (may have .est/.gt suffixes)
    est_col_actual <- if (estimate_col %in% colnames(merged_df)) {
      estimate_col
    } else if (paste0(estimate_col, ".est") %in% colnames(merged_df)) {
      paste0(estimate_col, ".est")
    } else {
      NA_character_
    }
    
    gt_col_actual <- if (ground_truth_col %in% colnames(merged_df)) {
      ground_truth_col
    } else if (paste0(ground_truth_col, ".gt") %in% colnames(merged_df)) {
      paste0(ground_truth_col, ".gt")
    } else {
      NA_character_
    }
    
    if (is.na(est_col_actual) || is.na(gt_col_actual)) {
      cat("WARNING: Column mismatch for arm", arm_name, "- estimate_col:", estimate_col, "-> ", est_col_actual, "| ground_truth_col:", ground_truth_col, "-> ", gt_col_actual, "\n")
      return(NULL)
    }
    
    data.frame(
      arm = arm_name,
      canonical_cell_type = cell_type_mapping$canonical_cell_type[i],
      estimate_cell_type = estimate_col,
      ground_truth_cell_type = ground_truth_col,
      n_samples = sum(complete.cases(merged_df[[est_col_actual]], merged_df[[gt_col_actual]])),
      correlation = cor(merged_df[[est_col_actual]], merged_df[[gt_col_actual]], use = "pairwise.complete.obs", method = "pearson"),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

first_arm_inputs <- prepare_estimates(
  file.path(opt$results_root, arm_names[1]),
  opt$proportions_file,
  opt$metadata_file
)

ground_truth_df <- if (opt$cohort == "CMC_MSSM") {
  load_ground_truth_cmc(opt$ground_truth_cmc)
} else {
  load_ground_truth_rosmap(first_arm_inputs$metadata, opt$ground_truth_rosmap)
}

arm_inputs_list <- lapply(arm_names, function(arm_name) {
  arm_dir <- file.path(opt$results_root, arm_name)
  prepare_estimates(arm_dir, opt$proportions_file, opt$metadata_file)
})
names(arm_inputs_list) <- arm_names

arm_cell_type_mappings <- lapply(arm_names, function(arm_name) {
  mapping_df <- match_cell_types(arm_inputs_list[[arm_name]]$estimates, ground_truth_df)
  
  if (nrow(mapping_df) == 0) {
    stop("ERROR: No shared cell types found for arm ", arm_name)
  }
  
  mapping_df
})
names(arm_cell_type_mappings) <- arm_names

common_benchmark_ids <- Reduce(intersect, lapply(arm_names, function(arm_name) {
  intersect(arm_inputs_list[[arm_name]]$estimates$benchmark_id, ground_truth_df$benchmark_id)
}))

if (length(common_benchmark_ids) == 0) {
  stop("ERROR: No shared benchmark samples found across all requested arms")
}

common_canonical_cell_types <- Reduce(intersect, lapply(arm_cell_type_mappings, function(mapping_df) {
  mapping_df$canonical_cell_type
}))

if (length(common_canonical_cell_types) == 0) {
  stop("ERROR: No shared cell types are common to all requested arms")
}

common_canonical_cell_types <- sort(common_canonical_cell_types)

cell_type_mapping_summary <- lapply(arm_names, function(arm_name) {
  arm_cell_type_mappings[[arm_name]] %>%
    filter(canonical_cell_type %in% common_canonical_cell_types) %>%
    mutate(arm = arm_name) %>%
    select(arm, canonical_cell_type, estimate_cell_type, ground_truth_cell_type)
}) %>% bind_rows()

all_results <- lapply(arm_names, function(arm_name) {
  arm_inputs <- arm_inputs_list[[arm_name]]
  cell_type_mapping <- arm_cell_type_mappings[[arm_name]] %>%
    filter(canonical_cell_type %in% common_canonical_cell_types)
  
  estimates_common <- arm_inputs$estimates %>%
    filter(benchmark_id %in% common_benchmark_ids)
  ground_truth_common <- ground_truth_df %>%
    filter(benchmark_id %in% common_benchmark_ids)
  
  compute_correlations(estimates_common, ground_truth_common, cell_type_mapping, arm_name)
}) %>% bind_rows()

arm_summary <- all_results %>%
  group_by(arm) %>%
  summarize(
    n_cell_types = n(),
    n_common_samples = length(common_benchmark_ids),
    median_correlation = median(correlation, na.rm = TRUE),
    iqr_correlation = IQR(correlation, na.rm = TRUE),
    n_positive_correlations = sum(correlation > 0, na.rm = TRUE),
    n_missing_correlations = sum(is.na(correlation)),
    .groups = "drop"
  )

baseline_arm <- if ("none" %in% arm_names) "none" else arm_names[1]

delta_df <- all_results %>%
  select(arm, canonical_cell_type, correlation) %>%
  pivot_wider(names_from = arm, values_from = correlation)

comparison_arms <- setdiff(arm_names, baseline_arm)
if (length(comparison_arms) > 0 && baseline_arm %in% colnames(delta_df)) {
  for (comparison_arm in comparison_arms) {
    if (comparison_arm %in% colnames(delta_df)) {
      delta_df[[paste0("delta_vs_", baseline_arm, "_", comparison_arm)]] <- delta_df[[comparison_arm]] - delta_df[[baseline_arm]]
    }
  }
}

improvement_summary <- if (length(comparison_arms) > 0 && baseline_arm %in% colnames(delta_df)) {
  lapply(comparison_arms, function(comparison_arm) {
    delta_col <- paste0("delta_vs_", baseline_arm, "_", comparison_arm)
    if (!delta_col %in% colnames(delta_df)) {
      return(NULL)
    }
    data.frame(
      baseline_arm = baseline_arm,
      comparison_arm = comparison_arm,
      improved_cell_types = sum(delta_df[[delta_col]] > 0, na.rm = TRUE),
      worse_cell_types = sum(delta_df[[delta_col]] < 0, na.rm = TRUE),
      unchanged_cell_types = sum(delta_df[[delta_col]] == 0, na.rm = TRUE),
      evaluable_cell_types = sum(!is.na(delta_df[[delta_col]])),
      median_delta = median(delta_df[[delta_col]], na.rm = TRUE),
      win_rate = sum(delta_df[[delta_col]] > 0, na.rm = TRUE) / sum(!is.na(delta_df[[delta_col]])),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
} else {
  data.frame()
}

display_levels <- all_results %>%
  group_by(canonical_cell_type) %>%
  summarize(median_correlation = median(correlation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_correlation)) %>%
  pull(canonical_cell_type)

bar_plot <- all_results %>%
  mutate(canonical_cell_type = factor(canonical_cell_type, levels = rev(display_levels))) %>%
  ggplot(aes(x = canonical_cell_type, y = correlation, fill = arm)) +
  geom_col(position = "dodge", na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  coord_flip() +
  labs(
    title = paste0(opt$cohort, ": MGP accuracy across pre-MGP tech-covariate arms"),
    x = "Cell type",
    y = "Pearson correlation",
    fill = "Arm"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(
  file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_accuracy_barplot.png")),
  bar_plot,
  width = 12,
  height = 10,
  dpi = 300
)

if (length(comparison_arms) > 0 && nrow(improvement_summary) > 0) {
  delta_long <- delta_df %>%
    select(canonical_cell_type, starts_with("delta_vs_")) %>%
    pivot_longer(
      cols = starts_with("delta_vs_"),
      names_to = "comparison",
      values_to = "delta"
    )
  
  delta_plot <- delta_long %>%
    mutate(canonical_cell_type = factor(canonical_cell_type, levels = rev(display_levels))) %>%
    ggplot(aes(x = canonical_cell_type, y = delta, fill = comparison)) +
    geom_col(position = "dodge", na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    coord_flip() +
    labs(
      title = paste0(opt$cohort, ": correlation delta versus ", baseline_arm),
      x = "Cell type",
      y = "Delta Pearson correlation",
      fill = "Comparison"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(
    file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_accuracy_delta.png")),
    delta_plot,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# Heatmap view is more robust for comparing methods than relying on mean summaries.
heatmap_plot <- all_results %>%
  mutate(canonical_cell_type = factor(canonical_cell_type, levels = rev(display_levels))) %>%
  ggplot(aes(x = arm, y = canonical_cell_type, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey80") +
  labs(
    title = paste0(opt$cohort, ": correlation heatmap by arm and cell type"),
    x = "Arm",
    y = "Cell type",
    fill = "Pearson r"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(
  file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_accuracy_heatmap.png")),
  heatmap_plot,
  width = 9,
  height = 10,
  dpi = 300
)

write.csv(all_results, file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_correlations_by_cell_type.csv")), row.names = FALSE)
write.csv(arm_summary, file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_arm_summary.csv")), row.names = FALSE)
write.csv(delta_df, file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_delta_summary.csv")), row.names = FALSE)
write.csv(improvement_summary, file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_improvement_summary.csv")), row.names = FALSE)
write.csv(cell_type_mapping_summary, file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_cell_type_mapping.csv")), row.names = FALSE)
write.csv(data.frame(benchmark_id = common_benchmark_ids), file.path(opt$output_dir, paste0(tolower(opt$cohort), "_tech_cov_common_samples.csv")), row.names = FALSE)

cat("\nPer-arm summary:\n")
print(arm_summary)
if (nrow(improvement_summary) > 0) {
  cat("\nImprovement summary:\n")
  print(improvement_summary)
}
cat("\nCommon samples across arms:", length(common_benchmark_ids), "\n")
cat("Common cell types across arms:", length(common_canonical_cell_types), "\n")
cat("\nOutputs written to:", opt$output_dir, "\n")
