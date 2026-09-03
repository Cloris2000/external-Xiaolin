# =============================================================================
# Figure 2. Bulk-derived cell-type proportion phenotypes and snRNA-seq validation
#
# Panels:
#   A — Validation scatter plots for VIP, L5.6.IT.Car3, Microglia  (top row)
#   B — All-cell-type Pearson r bar plot, mean ± SD across cohorts  (bottom row)
#
# Validation cohorts: ROSMAP, HBCC, MSBB — paired to Hodge label-transferred
# snRNA-seq proportions (unified 19-subclass taxonomy matching bulk MGP).
# Layout: A on top (accuracy) / B on bottom (scatter grid)
#
# Output:
#   figure2_celltype_validation_focus.png / .pdf / .svg
#   combined_bulk_snrna_paired.tsv
#   figure2_celltype_accuracy.tsv
# =============================================================================

# =============================================================================
# SECTION 1 — File paths and parameters (edit here)
# =============================================================================

DATA_DIR    <- "/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow"
RESULTS_DIR <- file.path(DATA_DIR, "results")
OUT_DIR     <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure"

COHORT_LIST <- c(
  "ROSMAP", "ROSMAP_array",
  "MSBB", "Mayo",
  "CMC_MSSM", "CMC_PENN", "CMC_PITT",
  "GTEx_v10", "NABEC", "GVEX",
  "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
  "AMP_AD_Mayo", "AMP_AD_Rush"
)
PROP_FILENAME <- "cell_proportions.csv"

# Focus cell types for scatter panel (clean names)
FOCUS_CELLTYPES <- c("vip", "l5_6_it_car3", "microglia")

# Drop from both validation panels (clean names)
EXCLUDE_CELLTYPES <- c("l4_it")

# --------------------------------------------------------------------------
# Optional: pre-built combined paired file (all cohorts).
# If NULL, auto-generate from ROSMAP + HBCC + MSBB below.
# --------------------------------------------------------------------------
VALID_PAIRED_FILE <- NULL

# Hodge label-transferred snRNA-seq cell proportions (QC'd, wide format)
HODGE_SN_DIR <- file.path(DATA_DIR, "data/snRNAseq_hodge_label_cell_prop")

# --- ROSMAP ---
# Join: bulk syn* → provenance specimenID → meta individualID (R#######) → sn
ROSMAP_BULK_FILE  <- file.path(RESULTS_DIR, "ROSMAP", "cell_proportions.csv")
ROSMAP_SN_FILE    <- file.path(HODGE_SN_DIR, "rosmap_green_cell_proportions_qc.csv")
ROSMAP_META_FILE  <- paste0(
  "/project/rrg-shreejoy/zhoux156/HBCC_ROSMAP_MSBB_bulkRNAseq/rosmap/metadata/",
  "RNAseq_Harmonization_ROSMAP_combined_metadata.csv")
# Provenance tables (syn* barcode → specimenID) — transferred with the nextflow
# working directory; expected under DATA_DIR/data/rosmap_provenance/ or the
# Synapse download location below. Update if your copy lives elsewhere.
ROSMAP_PROV_DIR   <- file.path(DATA_DIR,
  "data/rosmap_provenance/Rosmap_Gene_Quantification")
ROSMAP_PROV_FILES <- c(
  file.path(ROSMAP_PROV_DIR, "Rosmap_Batch1_Stranded/ROSMAP_batch1_provenance.csv"),
  file.path(ROSMAP_PROV_DIR, "Rosmap_Batch2_Stranded/ROSMAP_batch2_provenance.csv"),
  file.path(ROSMAP_PROV_DIR, "Rosmap_Batch3_Stranded/ROSMAP_batch3_provenance.csv"),
  file.path(ROSMAP_PROV_DIR, "Rosmap_Batch4_Stranded/ROSMAP_batch4_provenance.csv")
)

# --- HBCC ---
# Join: sn AMPAD_HBCC → CMC_HBCC (crosswalk) → RNA map → bulk
HBCC_BULK_FILE      <- file.path(RESULTS_DIR, "NIMH_HBCC_1M", "cell_proportions.csv")
HBCC_SN_FILE        <- file.path(HODGE_SN_DIR, "psychad_hbcc_cell_proportions_qc.csv")
HBCC_MAP_FILE       <- file.path(DATA_DIR, "data_input/nimh_hbcc/HBCC_rna_wgs_id_mapping.csv")
HBCC_CROSSWALK_FILE <- file.path(DATA_DIR, "data_input/sn_hbcc/ampad_to_cmc_crosswalk.csv")

# --- MSBB ---
# Join: bulk → RNA meta individualID → bridge SubID_export_synapse → sn AMPAD_MSSM
MSBB_BULK_FILE   <- file.path(RESULTS_DIR, "MSBB", "cell_proportions.csv")
MSBB_SN_FILE     <- file.path(HODGE_SN_DIR, "psychad_mssm_cell_proportions_qc.csv")
MSBB_RNA_META    <- paste0(
  "/project/rrg-shreejoy/zhoux156/HBCC_ROSMAP_MSBB_bulkRNAseq/msbb/metadata/",
  "RNAseq_Harmonization_MSBB_combined_metadata.csv")
MSBB_BRIDGE_FILE <- file.path(RESULTS_DIR, "MSBB_sn/msbb_sn_wgs_bridge_224.tsv")

# --- MATHYS (ROSMAP Mathys Hodge label-transferred) ---
# Join: same ROSMAP bulk + provenance as the ROSMAP section; only the snRNA file differs.
# Mathys snRNA was built from mathys_hodge_subclass_fine_hodge19_labels.csv by
# build_mathys_proportions.py; individual IDs are R* (same ROSMAP namespace).
MATHYS_SN_FILE <- file.path(DATA_DIR,
  "data_input/sn_rosmap_mathys_hodge_wgs/cell_proportions.csv")
# ROSMAP_BULK_FILE, ROSMAP_PROV_FILES, ROSMAP_META_FILE reused from the ROSMAP block above.

# --- RUZICKA (CMC MSSM / Ruzicka et al. 2024) ---
# Join: direct individual_id match — Ruzicka snRNA uses CMC_MSSM_NNN IDs that
# coincide with the individualID column in the CMC_MSSM bulk pipeline output.
RUZ_BULK_FILE <- file.path(RESULTS_DIR, "CMC_MSSM", "cell_proportions.csv")
RUZ_SN_FILE   <- file.path(DATA_DIR,
  "data_input/sn_ruz_mssm/cell_proportions.csv")

# Cohort colours — Okabe-Ito colorblind-safe palette (no red or green)
COHORT_COLORS <- c(
  ROSMAP  = "#0072B2",   # blue
  HBCC    = "#E69F00",   # orange
  MSBB    = "#CC79A7",   # mauve/pink
  Mathys  = "#56B4E9",   # sky blue
  Ruzicka = "#009E73"    # bluish green
)

# Figure dimensions
FIG_WIDTH  <- 15
FIG_HEIGHT <- 11
FIG_DPI    <- 300

# Reliability threshold (dashed guide) for estimation-accuracy panel
REL_THRESHOLD <- 0.3

# =============================================================================
# SECTION 2 — Package loading
# =============================================================================

required_pkgs <- c("tidyverse", "ggplot2", "cowplot", "scales",
                   "RColorBrewer", "readr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  Installing: ", pkg); install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

HAS_JANITOR <- requireNamespace("janitor", quietly = TRUE)
if (HAS_JANITOR) {
  suppressPackageStartupMessages(library(janitor))
} else {
  message("Note: 'janitor' not found; using base-R fallback.")
}

clean_names_fn <- function(df) {
  if (HAS_JANITOR) return(janitor::clean_names(df))
  nms <- sub("^_+|_+$", "",
             tolower(gsub("[^a-zA-Z0-9]+", "_",
                          gsub("([a-z])([A-Z])", "\\1_\\2", names(df)))))
  names(df) <- nms; df
}
make_clean_fn <- function(x) {
  if (HAS_JANITOR) return(janitor::make_clean_names(x))
  sub("^_+|_+$", "",
      tolower(gsub("[^a-zA-Z0-9]+", "_",
                   gsub("([a-z])([A-Z])", "\\1_\\2", x))))
}

cat("Packages loaded.\n")

# =============================================================================
# SECTION 3 — Helper functions
# =============================================================================

read_tabular <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported extension: ", path)
  }
}

make_placeholder <- function(msg, title = "") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 3.8, hjust = 0.5, vjust = 0.5,
             color = "grey45", fontface = "italic") +
    labs(title = title) +
    theme_void() +
    theme(panel.border    = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
          plot.title      = element_text(size = 11, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

pearson_stats <- function(x, y) {
  ok  <- complete.cases(x, y)
  ct  <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))
}

# =============================================================================
# SECTION 4 — Load and clean bulk proportions
# =============================================================================

cat("\n--- Loading bulk proportions ---\n")

prop_list <- lapply(COHORT_LIST, function(cohort) {
  path <- file.path(RESULTS_DIR, cohort, PROP_FILENAME)
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  df  <- read_tabular(path)
  df  <- clean_names_fn(df)
  sid <- grep("^specimen_?id$|^sample_?id$", names(df),
              value = TRUE, ignore.case = TRUE)[1]
  if (is.na(sid)) { warning("No ID column in ", cohort); return(NULL) }
  df  <- dplyr::rename(df, sample_id = !!sid)
  df$cohort <- cohort; df
})
prop_list <- Filter(Negate(is.null), prop_list)
prop_wide <- dplyr::bind_rows(prop_list)
cat("  Donors:", nrow(prop_wide), "| Cohorts:", length(prop_list), "\n")

ct_cols   <- setdiff(names(prop_wide), c("sample_id", "cohort"))
prop_long <- prop_wide %>%
  tidyr::pivot_longer(cols = all_of(ct_cols),
                      names_to = "cell_type", values_to = "proportion") %>%
  dplyr::filter(!is.na(sample_id), !is.na(cohort),
                !is.na(cell_type), !is.na(proportion))

max_prop <- max(prop_long$proportion, na.rm = TRUE)
has_neg  <- any(prop_long$proportion < 0, na.rm = TRUE)
if (max_prop <= 1 && !has_neg) {
  prop_long$proportion_pct <- prop_long$proportion * 100
  BULK_XLAB <- "Bulk-derived cell-type proportion (%)"
} else {
  prop_long$proportion_pct <- prop_long$proportion
  BULK_XLAB <- "Bulk-derived cell-type proportion estimate (MGP)"
}
prop_long$cohort <- factor(prop_long$cohort, levels = COHORT_LIST)

# =============================================================================
# SECTION 5 — Cell-type ordering, labels, colors, base theme
# =============================================================================

CELLTYPE_ORDER <- c(
  "oligodendrocyte", "opc", "astrocyte", "microglia",
  "endothelial", "pericyte", "vlmc",
  "it", "l4_it", "l5_6_np", "l5_et", "l6_ct", "l5_6_it_car3", "l6b",
  "pvalb", "sst", "vip", "lamp5", "pax6"
)
present_ct <- unique(prop_long$cell_type)
ordered_ct <- c(CELLTYPE_ORDER[CELLTYPE_ORDER %in% present_ct],
                setdiff(present_ct, CELLTYPE_ORDER))
prop_long$cell_type <- factor(prop_long$cell_type, levels = ordered_ct)

# Display labels — Title-case for full-word cell classes; keep acronyms
# (OPC, VIP, SST, PVALB, LAMP5, PAX6, VLMC) and layer nomenclature uppercase.
CT_DISPLAY <- c(
  oligodendrocyte = "Oligodendrocyte",
  opc             = "OPC",
  astrocyte       = "Astrocyte",
  microglia       = "Microglia",
  endothelial     = "Endothelial",
  pericyte        = "Pericyte",
  vlmc            = "VLMC",
  it              = "IT",
  l4_it           = "L4 IT",
  l5_6_np         = "L5/6 NP",
  l5_et           = "L5 ET",
  l6_ct           = "L6 CT",
  l5_6_it_car3    = "L5/6 IT Car3",
  l6b             = "L6b",
  pvalb           = "PVALB",
  sst             = "SST",
  vip             = "VIP",
  lamp5           = "LAMP5",
  pax6            = "PAX6"
)
# Fallback: Title-case any cell type not explicitly mapped above
missing_ct <- setdiff(ordered_ct, names(CT_DISPLAY))
if (length(missing_ct) > 0) {
  fb <- vapply(gsub("_", " ", missing_ct), function(s)
    paste0(toupper(substring(s, 1, 1)), substring(s, 2)), character(1))
  CT_DISPLAY <- c(CT_DISPLAY, setNames(fb, missing_ct))
}
# Focus display labels
FOCUS_DISPLAY <- CT_DISPLAY[FOCUS_CELLTYPES]

# Broad cell-class mapping (for grouping/ordering the panels)
CT_CLASS <- c(
  oligodendrocyte = "Non-neuronal", opc = "Non-neuronal", astrocyte = "Non-neuronal",
  microglia = "Non-neuronal", endothelial = "Non-neuronal", pericyte = "Non-neuronal",
  vlmc = "Non-neuronal",
  it = "Excitatory", l4_it = "Excitatory", l5_6_np = "Excitatory", l5_et = "Excitatory",
  l6_ct = "Excitatory", l5_6_it_car3 = "Excitatory", l6b = "Excitatory",
  pvalb = "Inhibitory", sst = "Inhibitory", vip = "Inhibitory", lamp5 = "Inhibitory",
  pax6 = "Inhibitory"
)
CLASS_LEVELS <- c("Excitatory", "Inhibitory", "Non-neuronal")
CLASS_COLORS <- c("Excitatory"   = "#D55E00",   # vermillion
                  "Inhibitory"   = "#009E73",   # bluish green
                  "Non-neuronal" = "#7570B3")   # muted purple

# Colors
pal_g <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[3:8])(7)
pal_e <- colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges")[3:8])(7)
pal_i <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[3:8])(5)
ct_colors <- setNames(c(pal_g, pal_e, pal_i)[seq_len(length(ordered_ct))], ordered_ct)
FOCUS_COLORS <- ct_colors[FOCUS_CELLTYPES]

# Base ggplot theme
base_th <- theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text       = element_text(size = 10),
    axis.title      = element_text(size = 12),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 12, face = "bold")
  )

cat("Cell types:", length(ordered_ct), "\n")

# =============================================================================
# SECTION 6 — Build paired bulk–snRNA validation data (ROSMAP + HBCC + MSBB)
#   snRNA source: Hodge label-transferred proportions (unified 19 subclasses).
# =============================================================================

cat("\n--- Building paired validation data (Hodge-labelled snRNA) ---\n")

all_paired <- list()

# Helper: build paired long-format from a wide merged data frame
make_paired_long <- function(merged_wide, sn_to_bulk, cohort_label) {
  rows <- lapply(names(sn_to_bulk), function(sn_ct) {
    bulk_ct <- sn_to_bulk[[sn_ct]]
    if (!bulk_ct %in% names(merged_wide) || !sn_ct %in% names(merged_wide)) return(NULL)
    df <- data.frame(
      sample_id         = merged_wide$join_id,
      validation_cohort = cohort_label,
      cell_type         = bulk_ct,
      bulk_proportion   = merged_wide[[bulk_ct]],
      snrna_proportion  = merged_wide[[sn_ct]],
      stringsAsFactors  = FALSE
    )
    df[!is.na(df$bulk_proportion) & !is.na(df$snrna_proportion), ]
  })
  dplyr::bind_rows(Filter(Negate(is.null), rows))
}

# Hodge columns match bulk names after clean_names. Prefix sn CT columns with
# sn__ so they do not collide with bulk columns during the join, and return
# the sn_col → bulk_col map used by make_paired_long.
prefix_sn_cts <- function(sn_df, bulk_df, id_cols_sn) {
  sn_cts   <- setdiff(names(sn_df), id_cols_sn)
  bulk_cts <- names(bulk_df)
  shared   <- intersect(sn_cts, bulk_cts)
  if (length(shared) == 0L)
    stop("No shared cell-type columns between snRNA and bulk after clean_names.")
  sn_df <- dplyr::rename_with(
    sn_df, ~ paste0("sn__", .x), .cols = dplyr::all_of(shared)
  )
  sn_to_bulk <- setNames(shared, paste0("sn__", shared))
  list(sn_df = sn_df, sn_to_bulk = sn_to_bulk, n_cts = length(shared))
}

# --- User-provided file (overrides auto-generation) ---
if (!is.null(VALID_PAIRED_FILE) && file.exists(VALID_PAIRED_FILE)) {
  df <- read_tabular(VALID_PAIRED_FILE)
  df <- clean_names_fn(df)
  if (all(c("sample_id","cell_type","bulk_proportion","snrna_proportion") %in% names(df))) {
    if (!"validation_cohort" %in% names(df)) df$validation_cohort <- "Validation"
    df$cell_type <- make_clean_fn(df$cell_type)
    all_paired[["user"]] <- df
    cat("  Loaded user file:", VALID_PAIRED_FILE, "\n")
  }
}

# --- ROSMAP ---
cat("  Auto-generating ROSMAP...\n")
tryCatch({
  prov_exist <- file.exists(ROSMAP_PROV_FILES)
  if (!any(prov_exist)) stop("No ROSMAP provenance files found.")
  prov <- dplyr::bind_rows(lapply(ROSMAP_PROV_FILES[prov_exist],
    function(f) readr::read_csv(f, show_col_types = FALSE)))
  prov_link <- prov %>%
    dplyr::select(synapse_id = id, specimen_id = specimenID) %>%
    dplyr::distinct(synapse_id, .keep_all = TRUE)
  # Hodge rosmap_green uses meta individualID (R#######), not projid
  meta_r <- readr::read_csv(ROSMAP_META_FILE, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(individualID), nchar(as.character(individualID)) > 0) %>%
    dplyr::select(specimen_id = specimenID, individualID) %>%
    dplyr::mutate(individualID = as.character(individualID)) %>%
    dplyr::distinct(specimen_id, .keep_all = TRUE)
  synid_ind <- dplyr::inner_join(prov_link, meta_r, by = "specimen_id") %>%
    dplyr::select(synapse_id, individual_id = individualID) %>%
    dplyr::mutate(individual_id = as.character(individual_id)) %>%
    dplyr::distinct(synapse_id, .keep_all = TRUE)

  sn_r <- clean_names_fn(readr::read_csv(ROSMAP_SN_FILE, show_col_types = FALSE))
  ind_col_r <- grep("^individual_?id$", names(sn_r), value = TRUE)[1]
  if (is.na(ind_col_r)) ind_col_r <- names(sn_r)[1]
  sn_r <- dplyr::rename(sn_r, individual_id = !!ind_col_r) %>%
    dplyr::mutate(individual_id = as.character(individual_id))

  bulk_r <- readr::read_csv(ROSMAP_BULK_FILE, show_col_types = FALSE)
  bulk_c <- clean_names_fn(bulk_r)
  sid_r  <- grep("^specimen_?id$|^sample_?id$", names(bulk_c), value = TRUE)[1]
  bulk_c <- dplyr::rename(bulk_c, specimen_id = !!sid_r)

  pref <- prefix_sn_cts(sn_r, bulk_c, id_cols_sn = "individual_id")
  sn_r <- pref$sn_df
  cat("    Shared cell types:", pref$n_cts, "\n")

  merged_r <- bulk_c %>%
    dplyr::inner_join(
      dplyr::rename(synid_ind, specimen_id = synapse_id),
      by = "specimen_id"
    ) %>%
    dplyr::inner_join(sn_r, by = "individual_id") %>%
    dplyr::rename(join_id = specimen_id)
  cat("    ROSMAP merged:", nrow(merged_r), "donors\n")
  all_paired[["ROSMAP"]] <- make_paired_long(merged_r, pref$sn_to_bulk, "ROSMAP")
}, error = function(e) message("  ROSMAP failed: ", conditionMessage(e)))

# --- HBCC ---
cat("  Auto-generating HBCC...\n")
tryCatch({
  bulk_h   <- clean_names_fn(readr::read_csv(HBCC_BULK_FILE, show_col_types = FALSE))
  map_h    <- readr::read_csv(HBCC_MAP_FILE, show_col_types = FALSE)
  xwalk_h  <- readr::read_csv(HBCC_CROSSWALK_FILE, show_col_types = FALSE)
  sn_h_raw <- readr::read_csv(HBCC_SN_FILE, show_col_types = FALSE)

  # Remap AMPAD_HBCC → CMC_HBCC before joining to the RNA specimen map
  sn_h_raw <- sn_h_raw %>%
    dplyr::inner_join(xwalk_h, by = c("individualID" = "ampad_id")) %>%
    dplyr::mutate(individualID = cmc_id) %>%
    dplyr::select(-cmc_id)
  cat("    After AMPAD→CMC crosswalk:", nrow(sn_h_raw), "donors\n")

  sn_h <- clean_names_fn(sn_h_raw)
  ind_col_h <- grep("^individual_?id$", names(sn_h), value = TRUE)[1]
  if (is.na(ind_col_h)) ind_col_h <- names(sn_h)[1]
  sn_h <- dplyr::rename(sn_h, individual_id = !!ind_col_h)

  rna_to_ind <- map_h %>%
    dplyr::select(specimen_id = RNA_specimenID, individual_id = individualID) %>%
    dplyr::distinct(specimen_id, .keep_all = TRUE)
  sid_h  <- grep("^specimen_?id$|^sample_?id$", names(bulk_h), value = TRUE)[1]
  bulk_h <- dplyr::rename(bulk_h, specimen_id = !!sid_h)

  pref <- prefix_sn_cts(sn_h, bulk_h, id_cols_sn = "individual_id")
  sn_h <- pref$sn_df
  cat("    Shared cell types:", pref$n_cts, "\n")

  merged_h <- bulk_h %>%
    dplyr::inner_join(rna_to_ind, by = "specimen_id") %>%
    dplyr::inner_join(sn_h, by = "individual_id") %>%
    dplyr::rename(join_id = specimen_id)
  cat("    HBCC merged:", nrow(merged_h), "donors\n")
  all_paired[["HBCC"]] <- make_paired_long(merged_h, pref$sn_to_bulk, "HBCC")
}, error = function(e) message("  HBCC failed: ", conditionMessage(e)))

# --- MSBB ---
cat("  Auto-generating MSBB...\n")
tryCatch({
  bulk_m   <- clean_names_fn(readr::read_csv(MSBB_BULK_FILE, show_col_types = FALSE))
  meta_m   <- readr::read_csv(MSBB_RNA_META, show_col_types = FALSE)
  bridge_m <- readr::read_tsv(MSBB_BRIDGE_FILE, show_col_types = FALSE)
  sn_m     <- clean_names_fn(readr::read_csv(MSBB_SN_FILE, show_col_types = FALSE))

  ind_col_m <- grep("^individual_?id$", names(sn_m), value = TRUE)[1]
  if (is.na(ind_col_m)) ind_col_m <- names(sn_m)[1]
  sn_m <- dplyr::rename(sn_m, individual_id = !!ind_col_m)

  rna_to_ind_m <- meta_m %>%
    dplyr::select(specimen_id = specimenID, individual_id = individualID) %>%
    dplyr::filter(!is.na(individual_id), nchar(as.character(individual_id)) > 0) %>%
    dplyr::mutate(individual_id = as.character(individual_id)) %>%
    dplyr::distinct(specimen_id, .keep_all = TRUE)
  # Restrict to WGS-bridged donors; sn IDs are already AMPAD_MSSM (= SubID)
  bridge_link <- bridge_m %>%
    dplyr::select(individual_id = SubID_export_synapse) %>%
    dplyr::mutate(individual_id = as.character(individual_id)) %>%
    dplyr::distinct(individual_id)

  sid_m  <- grep("^specimen_?id$|^sample_?id$", names(bulk_m), value = TRUE)[1]
  bulk_m <- dplyr::rename(bulk_m, specimen_id = !!sid_m)

  pref <- prefix_sn_cts(sn_m, bulk_m, id_cols_sn = "individual_id")
  sn_m <- pref$sn_df
  cat("    Shared cell types:", pref$n_cts, "\n")

  merged_m <- bulk_m %>%
    dplyr::inner_join(rna_to_ind_m, by = "specimen_id") %>%
    dplyr::inner_join(bridge_link,  by = "individual_id") %>%
    dplyr::inner_join(sn_m,         by = "individual_id") %>%
    dplyr::rename(join_id = specimen_id)
  cat("    MSBB merged:", nrow(merged_m), "donors\n")
  all_paired[["MSBB"]] <- make_paired_long(merged_m, pref$sn_to_bulk, "MSBB")
}, error = function(e) message("  MSBB failed: ", conditionMessage(e)))

# --- MATHYS ---
# Re-uses the ROSMAP bulk + provenance tables (same ROSMAP donors, different snRNA).
cat("  Auto-generating Mathys...\n")
tryCatch({
  if (!file.exists(MATHYS_SN_FILE)) stop("Mathys snRNA file not found: ", MATHYS_SN_FILE)
  prov_exist_m2 <- file.exists(ROSMAP_PROV_FILES)
  if (!any(prov_exist_m2)) stop("No ROSMAP provenance files found.")
  prov_m2 <- dplyr::bind_rows(lapply(ROSMAP_PROV_FILES[prov_exist_m2],
    function(f) readr::read_csv(f, show_col_types = FALSE)))
  prov_link_m2 <- prov_m2 %>%
    dplyr::select(synapse_id = id, specimen_id = specimenID) %>%
    dplyr::distinct(synapse_id, .keep_all = TRUE)
  meta_m2 <- readr::read_csv(ROSMAP_META_FILE, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(individualID), nchar(as.character(individualID)) > 0) %>%
    dplyr::select(specimen_id = specimenID, individualID) %>%
    dplyr::mutate(individualID = as.character(individualID)) %>%
    dplyr::distinct(specimen_id, .keep_all = TRUE)
  synid_ind_m2 <- dplyr::inner_join(prov_link_m2, meta_m2, by = "specimen_id") %>%
    dplyr::select(synapse_id, individual_id = individualID) %>%
    dplyr::mutate(individual_id = as.character(individual_id)) %>%
    dplyr::distinct(synapse_id, .keep_all = TRUE)

  sn_mathys <- clean_names_fn(readr::read_csv(MATHYS_SN_FILE, show_col_types = FALSE))
  ind_col_m2 <- grep("^individual_?id$", names(sn_mathys), value = TRUE)[1]
  if (is.na(ind_col_m2)) ind_col_m2 <- names(sn_mathys)[1]
  sn_mathys <- dplyr::rename(sn_mathys, individual_id = !!ind_col_m2) %>%
    dplyr::mutate(individual_id = as.character(individual_id))

  bulk_m2 <- readr::read_csv(ROSMAP_BULK_FILE, show_col_types = FALSE)
  bulk_c2  <- clean_names_fn(bulk_m2)
  sid_m2   <- grep("^specimen_?id$|^sample_?id$", names(bulk_c2), value = TRUE)[1]
  bulk_c2  <- dplyr::rename(bulk_c2, specimen_id = !!sid_m2)

  pref_m2 <- prefix_sn_cts(sn_mathys, bulk_c2, id_cols_sn = "individual_id")
  sn_mathys <- pref_m2$sn_df
  cat("    Shared cell types:", pref_m2$n_cts, "\n")

  merged_m2 <- bulk_c2 %>%
    dplyr::inner_join(
      dplyr::rename(synid_ind_m2, specimen_id = synapse_id),
      by = "specimen_id"
    ) %>%
    dplyr::inner_join(sn_mathys, by = "individual_id") %>%
    dplyr::rename(join_id = specimen_id)
  cat("    Mathys merged:", nrow(merged_m2), "donors\n")
  all_paired[["Mathys"]] <- make_paired_long(merged_m2, pref_m2$sn_to_bulk, "Mathys")
}, error = function(e) message("  Mathys failed: ", conditionMessage(e)))

# --- RUZICKA ---
# Direct join: Ruzicka snRNA uses CMC_MSSM_NNN individual IDs that match the
# individualID column in the CMC_MSSM bulk cell_proportions.csv.
cat("  Auto-generating Ruzicka...\n")
tryCatch({
  if (!file.exists(RUZ_BULK_FILE)) stop("CMC_MSSM bulk file not found: ", RUZ_BULK_FILE)
  if (!file.exists(RUZ_SN_FILE))   stop("Ruzicka snRNA file not found: ", RUZ_SN_FILE)

  bulk_ruz <- clean_names_fn(readr::read_csv(RUZ_BULK_FILE, show_col_types = FALSE))
  sn_ruz   <- clean_names_fn(readr::read_csv(RUZ_SN_FILE,   show_col_types = FALSE))

  # CMC_MSSM bulk: ID column may be named individual_id, specimen_id, or sample_id
  sid_ruz <- grep("^specimen_?id$|^sample_?id$|^individual_?id$",
                  names(bulk_ruz), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(sid_ruz)) sid_ruz <- names(bulk_ruz)[1]
  bulk_ruz <- dplyr::rename(bulk_ruz, individual_id = !!sid_ruz) %>%
    dplyr::mutate(individual_id = as.character(individual_id))

  ind_col_ruz <- grep("^individual_?id$", names(sn_ruz), value = TRUE)[1]
  if (is.na(ind_col_ruz)) ind_col_ruz <- names(sn_ruz)[1]
  sn_ruz <- dplyr::rename(sn_ruz, individual_id = !!ind_col_ruz) %>%
    dplyr::mutate(individual_id = as.character(individual_id))

  pref_ruz <- prefix_sn_cts(sn_ruz, bulk_ruz, id_cols_sn = "individual_id")
  sn_ruz   <- pref_ruz$sn_df
  cat("    Shared cell types:", pref_ruz$n_cts, "\n")

  merged_ruz <- bulk_ruz %>%
    dplyr::inner_join(sn_ruz, by = "individual_id") %>%
    dplyr::rename(join_id = individual_id)
  cat("    Ruzicka merged:", nrow(merged_ruz), "donors\n")
  all_paired[["Ruzicka"]] <- make_paired_long(merged_ruz, pref_ruz$sn_to_bulk, "Ruzicka")
}, error = function(e) message("  Ruzicka failed: ", conditionMessage(e)))

# --- Combine all cohorts ---
paired_val <- dplyr::bind_rows(all_paired)

# Drop degenerate snRNA samples where a single cell type accounts for ~100% of
# nuclei (biologically implausible; indicates a failed/near-empty library).
# One such ROSMAP endothelial donor (snRNA prop = 100%) otherwise dominates the
# endothelial fit. This is applied at the source so the accuracy table and both
# panels stay consistent.
.n_before <- nrow(paired_val)
paired_val <- dplyr::filter(paired_val, snrna_proportion < 0.99)
if (nrow(paired_val) < .n_before)
  cat("  Removed", .n_before - nrow(paired_val),
      "degenerate snRNA row(s) with proportion >= 99%\n")

# Drop excluded cell types from both panels / accuracy table
if (length(EXCLUDE_CELLTYPES) > 0L) {
  .n_ex <- sum(paired_val$cell_type %in% EXCLUDE_CELLTYPES)
  paired_val <- dplyr::filter(paired_val, !cell_type %in% EXCLUDE_CELLTYPES)
  if (.n_ex > 0L)
    cat("  Excluded cell types:", paste(EXCLUDE_CELLTYPES, collapse = ", "),
        "(", .n_ex, "rows)\n")
}

HAS_PAIRED <- nrow(paired_val) > 0

if (HAS_PAIRED) {
  for (coh in names(all_paired)) {
    n_d <- length(unique(all_paired[[coh]]$sample_id))
    n_c <- length(unique(all_paired[[coh]]$cell_type))
    cat("  ", coh, ":", n_d, "donors,", n_c, "cell types\n")
  }
  readr::write_tsv(paired_val, file.path(OUT_DIR, "combined_bulk_snrna_paired.tsv"))
  cat("  Saved: combined_bulk_snrna_paired.tsv\n")
}
cat("  Total paired rows:", nrow(paired_val), "\n")

# =============================================================================
# SECTION 7 — Compute accuracy metrics per validation_cohort × cell_type
#   Two scale-invariant metrics are used: Pearson r and Spearman rho.
#   (RMSE is intentionally NOT used: bulk values are MGP scores on a relative
#    scale, not 0-1 proportions, so an absolute-error metric is not meaningful.)
# =============================================================================

if (HAS_PAIRED) {
  acc_df <- paired_val %>%
    dplyr::group_by(validation_cohort, cell_type) %>%
    dplyr::summarise(
      pearson_r    = suppressWarnings(
        cor(bulk_proportion, snrna_proportion, method = "pearson",
            use = "pairwise.complete.obs")),
      spearman_rho = suppressWarnings(
        cor(bulk_proportion, snrna_proportion, method = "spearman",
            use = "pairwise.complete.obs")),
      n_donors     = dplyr::n(),
      .groups      = "drop"
    ) %>%
    dplyr::filter(!is.na(pearson_r))

  # Mean metrics across cohorts for each cell type
  acc_mean_df <- acc_df %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      mean_r       = mean(pearson_r,    na.rm = TRUE),
      sd_r         = sd(pearson_r,      na.rm = TRUE),
      mean_rho     = mean(spearman_rho, na.rm = TRUE),
      sd_rho       = sd(spearman_rho,   na.rm = TRUE),
      n_cohorts    = dplyr::n(),
      total_donors = sum(n_donors, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    dplyr::mutate(
      ct_clean  = make_clean_fn(cell_type),
      is_focus  = ct_clean %in% FOCUS_CELLTYPES
    ) %>%
    dplyr::arrange(mean_r)

  readr::write_tsv(acc_mean_df, file.path(OUT_DIR, "figure2_celltype_accuracy.tsv"))
  cat("  Accuracy summary rows:", nrow(acc_mean_df), "\n")
}

# =============================================================================
# SECTION 9 — Panel A: Estimation accuracy across ALL cell types (r and rho)
#   Comprehensive, even-handed benchmark (no cell type singled out). Diamond =
#   mean across cohorts; small dots = per-cohort values; whiskers = ±SD.
# =============================================================================

cat("--- Building Panel A (accuracy metrics) ---\n")

panel_A <- tryCatch({
  if (!HAS_PAIRED) stop("No paired validation data available.")

  # Order within class by mean Pearson r (desc); classes as labelled blocks
  ord_tbl <- acc_mean_df %>%
    dplyr::mutate(class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
    dplyr::arrange(class_f, dplyr::desc(mean_r))
  ct_order_a <- ord_tbl$ct_clean

  mean_a <- acc_mean_df %>%
    dplyr::mutate(ct_ord  = factor(ct_clean, levels = ct_order_a),
                  class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS))

  coh_a <- acc_df %>%
    dplyr::mutate(ct_clean = make_clean_fn(cell_type),
                  ct_ord   = factor(ct_clean, levels = ct_order_a),
                  class_f  = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
    dplyr::filter(!is.na(ct_ord))

  has_se_a <- any(!is.na(mean_a$sd_r) & mean_a$sd_r > 0)

  ggplot(mean_a, aes(x = ct_ord, y = mean_r)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey55") +
    { if (has_se_a)
      geom_errorbar(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r),
                    width = 0.25, linewidth = 0.45, color = "grey60") } +
    geom_point(data = coh_a,
               aes(x = ct_ord, y = pearson_r, color = validation_cohort),
               size = 1.7, alpha = 0.85, shape = 16, inherit.aes = FALSE) +
    geom_point(size = 3.0, shape = 18, color = "#08519c") +
    facet_grid(cols = vars(class_f), scales = "free_x", space = "free_x") +
    scale_color_manual(values = COHORT_COLORS, name = "Cohort",
                       breaks = names(COHORT_COLORS),
                       guide = guide_legend(override.aes = list(size = 3))) +
    scale_x_discrete(labels = function(x) {
      lbl <- CT_DISPLAY[x]; ifelse(is.na(lbl), x, lbl)
    }) +
    labs(x = NULL, y = "Pearson r") +
    base_th +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
          axis.text.x        = element_text(size = 10, angle = 45, hjust = 1),
          strip.text         = element_text(size = 11, face = "bold"),
          strip.background   = element_rect(fill = "grey94", color = NA),
          panel.spacing.x    = unit(0.4, "lines"),
          legend.position    = "none",
          panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4))

}, error = function(e) {
  message("Panel A: ", conditionMessage(e))
  make_placeholder("Validation data not available.", "")
})

# =============================================================================
# SECTION 10 — Panel B: Per-cell-type scatter grid (small multiples, raw units)
#   Actual values shown: snRNA-seq proportion (%) on x, bulk MGP score (AU) on y.
#   Free per-facet scales so each cell type spans its own native range.
#   Per-facet r label is the cross-cohort mean r (identical to Panel A).
#   Facets grouped by class, then mean r (desc) within class; coloured by class.
# =============================================================================

cat("--- Building Panel B (per-cell-type scatter grid) ---\n")

panel_B <- tryCatch({
  if (!HAS_PAIRED) stop("No paired validation data available.")

  zdat <- paired_val %>%
    dplyr::mutate(ct_clean = make_clean_fn(cell_type),
                  sn_pct   = snrna_proportion * 100,   # fraction -> %
                  bulk_au  = bulk_proportion,          # MGP score (AU)
                  ct_lab   = CT_DISPLAY[ct_clean],
                  ct_lab   = ifelse(is.na(ct_lab), ct_clean, ct_lab)) %>%
    dplyr::filter(is.finite(sn_pct), is.finite(bulk_au))

  # Facet order: grouped by class, then mean Pearson r (desc) within class
  ct_order_b <- acc_mean_df %>%
    dplyr::mutate(class_f = factor(CT_CLASS[ct_clean], levels = CLASS_LEVELS)) %>%
    dplyr::arrange(class_f, dplyr::desc(mean_r)) %>%
    dplyr::pull(ct_clean)
  lab_levels <- unname(CT_DISPLAY[ct_order_b]); lab_levels <- ifelse(is.na(lab_levels), ct_order_b, lab_levels)
  zdat$ct_lab   <- factor(zdat$ct_lab, levels = lab_levels)
  zdat$class_f  <- factor(CT_CLASS[zdat$ct_clean], levels = CLASS_LEVELS)

  # Per-facet label uses the SAME cross-cohort mean r reported in Panel A
  r_facet <- acc_mean_df %>%
    dplyr::mutate(ct_lab = ifelse(is.na(CT_DISPLAY[ct_clean]), ct_clean, CT_DISPLAY[ct_clean]),
                  ct_lab = factor(ct_lab, levels = lab_levels),
                  lbl    = paste0("r = ", sprintf("%.2f", mean_r))) %>%
    dplyr::select(ct_lab, lbl) %>%
    dplyr::filter(!is.na(ct_lab))

  # Place label at top-left of each free panel
  pos_df <- zdat %>%
    dplyr::group_by(ct_lab) %>%
    dplyr::summarise(xpos = min(sn_pct, na.rm = TRUE),
                     ypos = max(bulk_au, na.rm = TRUE), .groups = "drop")
  r_facet <- dplyr::left_join(r_facet, pos_df, by = "ct_lab")

  ggplot(zdat, aes(x = sn_pct, y = bulk_au, color = validation_cohort)) +
    geom_point(size = 0.5, alpha = 0.28) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 0.8) +
    geom_text(data = r_facet,
              aes(x = xpos, y = ypos, label = lbl),
              inherit.aes = FALSE, hjust = 0, vjust = 1,
              size = 3.2, fontface = "bold", color = "grey15") +
    facet_wrap(~ ct_lab, nrow = 2, scales = "free") +
    scale_color_manual(values = COHORT_COLORS, name = "Cohort",
                       breaks = names(COHORT_COLORS),
                       guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "snRNA-seq-derived proportion (%)",
         y = "Bulk-derived proportion (AU)") +
    base_th +
    theme(strip.text       = element_text(size = 9.5, face = "bold"),
          strip.background = element_rect(fill = "grey94", color = NA),
          panel.spacing    = unit(0.55, "lines"),
          axis.text        = element_text(size = 7.5),
          legend.position  = "bottom",
          panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.4))

}, error = function(e) {
  message("Panel B: ", conditionMessage(e))
  make_placeholder("Validation data not available.",
                   "B. Bulk vs snRNA-seq proportion, per cell type")
})

# =============================================================================
# SECTION 11 — Combine and save
# =============================================================================

cat("\n--- Combining panels ---\n")

# Panel A (accuracy metrics, all cell types) on top; Panel B (combined
# standardized scatter) below. Each keeps its own legend (cohort / cell type).
fig2 <- cowplot::plot_grid(
  panel_A, panel_B,
  nrow          = 2,
  rel_heights   = c(0.72, 1.0),
  labels        = c("A", "B"),
  label_size    = 15,
  label_fontface = "bold"
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

out_png <- file.path(OUT_DIR, "figure2_celltype_validation_focus.png")
out_pdf <- file.path(OUT_DIR, "figure2_celltype_validation_focus.pdf")
out_svg <- file.path(OUT_DIR, "figure2_celltype_validation_focus.svg")

ggsave(out_png, fig2, width = FIG_WIDTH, height = FIG_HEIGHT,
       dpi = FIG_DPI, bg = "white")
cat("Saved PNG:", out_png, "\n")

ggsave(out_pdf, fig2, width = FIG_WIDTH, height = FIG_HEIGHT, bg = "white")
cat("Saved PDF:", out_pdf, "\n")

# SVG vector output
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(out_svg, fig2, width = FIG_WIDTH, height = FIG_HEIGHT,
         device = svglite::svglite)
  cat("Saved SVG:", out_svg, "\n")
} else {
  grDevices::svg(out_svg, width = FIG_WIDTH, height = FIG_HEIGHT)
  print(fig2)
  grDevices::dev.off()
  cat("Saved SVG:", out_svg, "\n")
}

cat("Done.\n")
