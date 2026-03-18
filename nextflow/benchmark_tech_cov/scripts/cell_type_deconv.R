#!/usr/bin/env Rscript
# Cell Type Proportion Estimation
# Flexible script supporting multiple deconvolution tools and reference taxonomies
# Supported tools: MGP, dtangle, CIBERSORT, MuSiC, Bisque, DWLS
# Supported taxonomies: sonny_markers, aging_paper, custom

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyverse)
  library(readr)
  library(data.table)
})

# Helper function to find column names with flexible matching
find_column <- function(df, possible_names, required = FALSE) {
  for (name in possible_names) {
    if (name %in% colnames(df)) {
      return(name)
    }
  }
  if (required) {
    stop(paste("Required column not found. Tried:", paste(possible_names, collapse = ", ")))
  }
  return(NULL)
}

# Parse command line arguments
option_list <- list(
  make_option(c("--corrected_data"), type="character", default=NULL,
              help="Path to corrected z-score data RData file"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Path to cleaned metadata CSV file"),
  make_option(c("--deconv_tool"), type="character", default="MGP",
              help="Deconvolution tool: MGP, dtangle, CIBERSORT, MuSiC, Bisque, DWLS (default: MGP)"),
  make_option(c("--reference_taxonomy"), type="character", default="sonny_markers",
              help="Reference taxonomy: sonny_markers, aging_paper, custom (default: sonny_markers)"),
  make_option(c("--marker_file"), type="character", default=NULL,
              help="Path or URL to marker file (required for custom taxonomy)"),
  make_option(c("--hgnc_mapping_file"), type="character", default=NULL,
              help="Path to HGNC mapping file (optional, for gene ID conversion)"),
  make_option(c("--marker_format"), type="character", default="auto",
              help="Marker file format: auto, sonny, aging, custom_csv, custom_tsv (default: auto)"),
  make_option(c("--cell_type_col"), type="character", default=NULL,
              help="Column name for cell type in marker file (auto-detect if not specified)"),
  make_option(c("--gene_col"), type="character", default=NULL,
              help="Column name for gene identifier in marker file (auto-detect if not specified)"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--proportions_output"), type="character", default="cell_proportions.csv",
              help="Output filename for cell proportions"),
  make_option(c("--proportions_scaled_output"), type="character", default="cell_proportions_scaled.csv",
              help="Output filename for scaled cell proportions")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("=== Cell Type Deconvolution Pipeline ===\n")
cat("Tool:", opt$deconv_tool, "\n")
cat("Reference taxonomy:", opt$reference_taxonomy, "\n")

# Load corrected data
cat("Loading corrected expression data...\n")
zscore_data_removedBatchEff_cov <- readRDS(opt$corrected_data)
cat("  Data dimensions:", dim(zscore_data_removedBatchEff_cov), "\n")
cat("  Data class:", class(zscore_data_removedBatchEff_cov), "\n")
cat("  Number of samples (cols):", ncol(zscore_data_removedBatchEff_cov), "\n")
cat("  Number of genes (rows):", nrow(zscore_data_removedBatchEff_cov), "\n")

# Helper function to validate gene IDs (Ensembl IDs)
validate_gene_ids <- function(gene_ids, context = "", min_ensembl_pct = 0.8) {
  if (is.null(gene_ids) || length(gene_ids) == 0) {
    stop("ERROR [", context, "]: Gene IDs are NULL or empty")
  }
  
  # Check if they're numeric (bad)
  if (all(grepl("^[0-9]+$", head(gene_ids, min(100, length(gene_ids)))))) {
    stop("ERROR [", context, "]: Gene IDs are numeric (1, 2, 3...) instead of Ensembl IDs. ",
         "First 10 IDs: ", paste(head(gene_ids, 10), collapse=", "))
  }
  
  # Check if they look like Ensembl IDs
  ensembl_count <- sum(grepl("^ENSG", gene_ids, ignore.case = TRUE))
  ensembl_pct <- ensembl_count / length(gene_ids)
  
  if (ensembl_pct < min_ensembl_pct) {
    warning("WARNING [", context, "]: Only ", round(ensembl_pct * 100, 1), 
            "% of gene IDs look like Ensembl IDs (expected at least ", round(min_ensembl_pct * 100), "%). ",
            "First 10 IDs: ", paste(head(gene_ids, 10), collapse=", "))
    if (ensembl_pct < 0.5) {
      stop("ERROR [", context, "]: Less than 50% of gene IDs are Ensembl IDs. ",
           "This suggests a serious problem with gene ID extraction.")
    }
  } else {
    cat("  ✓ Validated: ", round(ensembl_pct * 100, 1), "% of gene IDs are Ensembl IDs\n", sep="")
  }
  
  return(TRUE)
}

# Check rownames - critical for marker matching!
cat("  Checking rownames...\n")
cat("    is.null(rownames()):", is.null(rownames(zscore_data_removedBatchEff_cov)), "\n")
cat("    length(rownames()):", length(rownames(zscore_data_removedBatchEff_cov)), "\n")
cat("    dimnames()[[1]] is null:", is.null(dimnames(zscore_data_removedBatchEff_cov)[[1]]), "\n")
if (!is.null(dimnames(zscore_data_removedBatchEff_cov)[[1]])) {
  cat("    length(dimnames()[[1]]):", length(dimnames(zscore_data_removedBatchEff_cov)[[1]]), "\n")
  cat("    First 5 dimnames[[1]]:", paste(head(dimnames(zscore_data_removedBatchEff_cov)[[1]], 5), collapse=", "), "\n")
}

# If rownames are NULL or empty, try to get them from dimnames
if (is.null(rownames(zscore_data_removedBatchEff_cov)) || length(rownames(zscore_data_removedBatchEff_cov)) == 0) {
  cat("    WARNING: rownames() is NULL or empty!\n")
  # Try to get rownames from dimnames
  if (!is.null(dimnames(zscore_data_removedBatchEff_cov)[[1]]) && length(dimnames(zscore_data_removedBatchEff_cov)[[1]]) > 0) {
    cat("    Using dimnames[[1]] as rownames\n")
    rownames(zscore_data_removedBatchEff_cov) <- dimnames(zscore_data_removedBatchEff_cov)[[1]]
  } else if (is.list(zscore_data_removedBatchEff_cov) && "genes" %in% names(zscore_data_removedBatchEff_cov)) {
    cat("    Data is a list, checking for 'genes' element\n")
    rownames(zscore_data_removedBatchEff_cov) <- zscore_data_removedBatchEff_cov$genes
  } else if (is.data.frame(zscore_data_removedBatchEff_cov) && "gene" %in% colnames(zscore_data_removedBatchEff_cov)) {
    cat("    Data is data.frame, using 'gene' column as rownames\n")
    rownames(zscore_data_removedBatchEff_cov) <- zscore_data_removedBatchEff_cov$gene
    zscore_data_removedBatchEff_cov <- zscore_data_removedBatchEff_cov[, !colnames(zscore_data_removedBatchEff_cov) %in% "gene"]
  } else {
    cat("    ERROR: Cannot find rownames in any form!\n")
    stop("ERROR: Expression data has no rownames (gene IDs). Cannot match markers. Check how corrected_data is saved in remove_tech_covar.R")
  }
}

# Verify rownames are now set
if (length(rownames(zscore_data_removedBatchEff_cov)) == 0 || is.null(rownames(zscore_data_removedBatchEff_cov))) {
  cat("    ERROR: Still no rownames after attempts to fix!\n")
  stop("ERROR: Expression data has no rownames. Cannot proceed with marker matching.")
} else {
  cat("  Rownames successfully set. Length:", length(rownames(zscore_data_removedBatchEff_cov)), "\n")
  cat("  First 5 rownames:", paste(head(rownames(zscore_data_removedBatchEff_cov), 5), collapse=", "), "\n")
  
  # IMPORTANT: Strip version suffixes from Ensembl IDs (e.g. ENSG00000122012.10 -> ENSG00000122012)
  # This ensures marker matching works correctly with unversioned marker files
  original_rownames <- rownames(zscore_data_removedBatchEff_cov)
  rownames(zscore_data_removedBatchEff_cov) <- sub("\\.[0-9]+$", "", original_rownames)
  cat("  After removing version suffixes - First 5 rownames:", paste(head(rownames(zscore_data_removedBatchEff_cov), 5), collapse=", "), "\n")
  
  # CHECKPOINT: Validate gene IDs after loading corrected_data
  validate_gene_ids(rownames(zscore_data_removedBatchEff_cov), "After loading corrected_data")
}

# Load metadata if provided
metadata <- NULL
if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
  metadata <- read.csv(opt$metadata)
}

# ============================================================================
# Load and prepare marker genes based on reference taxonomy
# ============================================================================
cat("Loading marker gene files...\n")

load_markers <- function(reference_taxonomy, marker_file, marker_format, hgnc_mapping_file, cell_type_col = NULL, gene_col = NULL) {
  markers <- NULL
  detected_cell_type_col <- cell_type_col
  detected_gene_col <- gene_col
  
  if (reference_taxonomy == "sonny_markers") {
    # Sonny markers format
    if (is.null(marker_file)) {
      marker_url <- 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'
      cat("Downloading Sonny markers from:", marker_url, "\n")
      markers <- read_csv(marker_url, show_col_types = FALSE)
    } else {
      if (startsWith(marker_file, "http")) {
        markers <- read_csv(marker_file, show_col_types = FALSE)
      } else {
        markers <- read_csv(marker_file, show_col_types = FALSE)
      }
    }
    
    # Standardize column names
    colnames(markers) <- colnames(markers) %>% make.names() %>% tolower()
    
    # Map with HGNC if available
    if (!is.null(hgnc_mapping_file) && file.exists(hgnc_mapping_file)) {
      cat("  Loading HGNC mapping file:", hgnc_mapping_file, "\n")
      hgnc_mapping <- read_tsv(hgnc_mapping_file, show_col_types = FALSE)
      cat("  HGNC mapping file loaded. Columns:", paste(colnames(hgnc_mapping), collapse=", "), "\n")
      cat("  HGNC mapping rows:", nrow(hgnc_mapping), "\n")
      
      # Find entrez ID column
      entrez_col <- find_column(markers, c("entrez.gene.id", "entrez_id", "entrezid", "entrez"))
      if (!is.null(entrez_col)) {
        cat("  Found entrez ID column:", entrez_col, "\n")
        cat("  Markers before HGNC merge:", nrow(markers), "rows\n")
        cat("  Markers with non-NA entrez IDs:", sum(!is.na(markers[[entrez_col]])), "\n")
        
        markers <- left_join(
          markers %>% dplyr::rename(entrez_id = !!sym(entrez_col)),
          hgnc_mapping %>% distinct(entrez_id, .keep_all = TRUE) %>%
            dplyr::select(entrez_id, ensembl_gene_id) %>%
            dplyr::rename(ensembl_id = ensembl_gene_id),
          by = "entrez_id"
        )
        cat("  Markers after HGNC merge:", nrow(markers), "rows\n")
        cat("  Markers with ensembl_id:", sum(!is.na(markers$ensembl_id)), "\n")
        cat("  Sample ensembl_ids:", paste(head(markers$ensembl_id[!is.na(markers$ensembl_id)], 5), collapse=", "), "\n")
      } else {
        cat("  WARNING: Could not find entrez ID column in markers. HGNC mapping skipped.\n")
      }
    } else {
      cat("  WARNING: HGNC mapping file not provided or not found. Ensembl ID mapping will not be performed.\n")
    }
    
  } else if (reference_taxonomy == "aging_paper") {
    # Aging paper markers format (example - adjust based on actual format)
    if (is.null(marker_file)) {
      stop("Marker file is required for aging_paper taxonomy")
    }
    if (startsWith(marker_file, "http")) {
      markers <- read_csv(marker_file, show_col_types = FALSE)
    } else {
      markers <- read_csv(marker_file, show_col_types = FALSE)
    }
    colnames(markers) <- colnames(markers) %>% make.names() %>% tolower()
    
  } else if (reference_taxonomy == "custom") {
    # Custom marker file
    if (is.null(marker_file)) {
      stop("Marker file is required for custom taxonomy")
    }
    
    # Determine file format
    if (marker_format == "auto") {
      if (endsWith(marker_file, ".tsv") || endsWith(marker_file, ".txt")) {
        markers <- read_tsv(marker_file, show_col_types = FALSE)
      } else {
        markers <- read_csv(marker_file, show_col_types = FALSE)
      }
    } else if (marker_format == "custom_csv") {
      markers <- read_csv(marker_file, show_col_types = FALSE)
    } else if (marker_format == "custom_tsv") {
      markers <- read_tsv(marker_file, show_col_types = FALSE)
    } else {
      markers <- read_csv(marker_file, show_col_types = FALSE)
    }
    
    # Auto-detect cell type and gene columns if not specified
    if (is.null(detected_cell_type_col)) {
      detected_cell_type_col <- find_column(markers, c("cell_type", "CellType", "cellType", "subclass", "Subclass", 
                                                         "cluster", "Cluster", "type", "Type"), required = TRUE)
    }
    
    if (is.null(detected_gene_col)) {
      detected_gene_col <- find_column(markers, c("gene", "Gene", "gene_symbol", "GeneSymbol", "geneSymbol",
                                                   "ensembl_id", "EnsemblID", "ensembl_gene_id", "gene_id"), required = TRUE)
    }
  } else {
    stop(paste("Unknown reference taxonomy:", reference_taxonomy))
  }
  
  return(list(markers = markers, cell_type_col = detected_cell_type_col, gene_col = detected_gene_col))
}

marker_data <- load_markers(opt$reference_taxonomy, opt$marker_file, opt$marker_format, opt$hgnc_mapping_file, 
                            opt$cell_type_col, opt$gene_col)
markers <- marker_data$markers
detected_cell_type_col <- marker_data$cell_type_col
detected_gene_col <- marker_data$gene_col

# ============================================================================
# Prepare marker list for deconvolution
# ============================================================================
cat("Preparing marker gene list...\n")

prepare_marker_list <- function(markers, reference_taxonomy, expr_rownames, cell_type_col = NULL, gene_col = NULL) {
  cat("  In prepare_marker_list function:\n")
  cat("    Expression rownames length:", length(expr_rownames), "\n")
  cat("    Expression rownames sample (first 5):", paste(head(expr_rownames, 5), collapse=", "), "\n")
  
  if (reference_taxonomy == "sonny_markers") {
    # Filter for markers used in MGP
    used_col <- find_column(markers, c("used.in.mgp", "used_in_mgp", "used"))
    if (!is.null(used_col)) {
      cat("    Filtering for markers with used_col:", used_col, "\n")
      markers_before_filter <- nrow(markers)
      markers <- markers %>% filter(!!sym(used_col) == "TRUE" | !!sym(used_col) == TRUE)
      cat("    Markers before filter:", markers_before_filter, ", after filter:", nrow(markers), "\n")
    } else {
      cat("    WARNING: Could not find 'used.in.mgp' column. Using all markers.\n")
    }
    
    subclass_col <- find_column(markers, c("subclass", "Subclass", "cell_type", "CellType"))
    if (is.null(subclass_col)) {
      stop("Cannot find cell type column in Sonny markers")
    }
    cat("    Using subclass column:", subclass_col, "\n")
    
    new_cell_types <- markers %>% filter(!is.na(!!sym(subclass_col))) %>% pull(!!sym(subclass_col)) %>% unique()
    cat("    Found", length(new_cell_types), "cell types\n")
    
    # Create marker list
    marker_list <- lapply(new_cell_types, function(cell_type) {
      # Try Ensembl ID first
      if ("ensembl_id" %in% colnames(markers)) {
        # Get all markers for this cell type (before filtering by expr_rownames)
        all_markers_for_ct <- markers %>%
          filter(!!sym(subclass_col) == cell_type, !is.na(ensembl_id)) %>%
          pull(ensembl_id)
        
        # Filter to only those in expression matrix
        marker_genes <- all_markers_for_ct[all_markers_for_ct %in% expr_rownames]
        
        # Debug output for first cell type
        if (cell_type == new_cell_types[1]) {
          cat("    Example for", cell_type, ":\n")
          cat("      All markers (with ensembl_id):", length(all_markers_for_ct), "\n")
          cat("      Markers in expression matrix:", length(marker_genes), "\n")
          if (length(all_markers_for_ct) > 0 && length(marker_genes) == 0) {
            cat("      First marker Ensembl ID:", head(all_markers_for_ct, 1), "\n")
            cat("      Is it in expr_rownames?", head(all_markers_for_ct, 1) %in% expr_rownames, "\n")
            cat("      Sample expr_rownames:", paste(head(expr_rownames, 3), collapse=", "), "\n")
          }
        }
      } else {
        # Fallback to gene names
        gene_col <- find_column(markers, c("gene", "Gene", "gene_symbol", "GeneSymbol"))
        if (!is.null(gene_col)) {
          marker_genes <- markers %>%
            filter(!!sym(subclass_col) == cell_type,
                   !!sym(gene_col) %in% expr_rownames) %>%
            pull(!!sym(gene_col))
        } else {
          stop("Cannot find gene identifier column in markers")
        }
      }
      return(marker_genes)
    })
    names(marker_list) <- new_cell_types
    
  } else if (reference_taxonomy %in% c("aging_paper", "custom")) {
    # Generic format: cell_type_col and gene_col
    if (is.null(cell_type_col)) {
      cell_type_col <- find_column(markers, c("cell_type", "CellType", "cellType", "subclass", "Subclass", 
                                               "cluster", "Cluster", "type", "Type"), required = TRUE)
    }
    if (is.null(gene_col)) {
      gene_col <- find_column(markers, c("gene", "Gene", "gene_symbol", "GeneSymbol", "geneSymbol",
                                         "ensembl_id", "EnsemblID", "ensembl_gene_id", "gene_id"), required = TRUE)
    }
    
    new_cell_types <- markers %>% filter(!is.na(!!sym(cell_type_col))) %>% pull(!!sym(cell_type_col)) %>% unique()
    
    marker_list <- lapply(new_cell_types, function(cell_type) {
      marker_genes <- markers %>%
        filter(!!sym(cell_type_col) == cell_type,
               !!sym(gene_col) %in% expr_rownames) %>%
        pull(!!sym(gene_col))
      return(marker_genes)
    })
    names(marker_list) <- new_cell_types
  }
  
  return(marker_list)
}

marker_list <- prepare_marker_list(markers, opt$reference_taxonomy, 
                                    rownames(zscore_data_removedBatchEff_cov),
                                    detected_cell_type_col, detected_gene_col)

# CHECKPOINT: Validate marker list has valid Ensembl IDs
all_marker_genes_check <- unique(unlist(marker_list))
if (length(all_marker_genes_check) > 0) {
  validate_gene_ids(all_marker_genes_check, "In marker_list after creation")
} else {
  stop("ERROR: marker_list is empty. No markers found for any cell type.")
}

cell_types <- names(marker_list)
cat("Cell types to estimate:", paste(cell_types, collapse = ", "), "\n")
cat("Number of cell types:", length(cell_types), "\n")

# ============================================================================
# Run deconvolution based on selected tool
# ============================================================================
cat("Running", opt$deconv_tool, "deconvolution...\n")

if (opt$deconv_tool == "MGP") {
  # MGP deconvolution
  if (!require("markerGeneProfile", quietly = TRUE)) {
    stop("markerGeneProfile package is required for MGP deconvolution. Please install it.")
  }
  
  # Prepare expression data
  cat("Preparing expression data for MGP...\n")
  matrix <- as.data.frame(zscore_data_removedBatchEff_cov)
  cat("  Matrix dimensions after conversion:", dim(matrix), "\n")
  
  gene_id_col <- ifelse("ensembl_id" %in% colnames(markers), "ensembl_id", 
                       find_column(markers, c("gene", "Gene", "gene_symbol")))
  
  if (gene_id_col == "ensembl_id" && all(rownames(matrix) %in% markers$ensembl_id)) {
    matrix <- matrix %>% rownames_to_column(var = "gene_symbol")
  } else {
    matrix <- matrix %>% rownames_to_column(var = "gene_symbol")
  }
  
  genes_only <- matrix %>% subset(gene_symbol != "")
  if (length(which(duplicated(genes_only$gene_symbol))) != 0) {
    genes_only <- genes_only[-which(duplicated(genes_only$gene_symbol)), ]
  }
  cat("  Genes_only dimensions:", dim(genes_only), "\n")
  cat("  Number of samples in genes_only:", ncol(genes_only) - 1, "\n")  # -1 for gene_symbol column
  
  # Debug: Check marker gene matching
  cat("  Checking marker gene matching...\n")
  cat("  Expression matrix rownames (first 5):", paste(head(rownames(zscore_data_removedBatchEff_cov), 5), collapse=", "), "\n")
  cat("  Expression matrix rownames format check:\n")
  expr_rownames_sample <- head(rownames(zscore_data_removedBatchEff_cov), 10)
  cat("    Sample rownames:", paste(expr_rownames_sample, collapse=", "), "\n")
  cat("    Rownames look like Ensembl IDs:", all(grepl("^ENSG", expr_rownames_sample)), "\n")
  
  all_marker_genes <- unique(unlist(marker_list))
  cat("    Total unique marker genes in marker_list:", length(all_marker_genes), "\n")
  if (length(all_marker_genes) > 0) {
    cat("    First few marker genes:", paste(head(all_marker_genes, 5), collapse=", "), "\n")
    cat("    Marker genes look like Ensembl IDs:", all(grepl("^ENSG", head(all_marker_genes, min(10, length(all_marker_genes))))), "\n")
  }
  
  matched_genes <- intersect(all_marker_genes, genes_only$gene_symbol)
  cat("    Matched marker genes in expression data:", length(matched_genes), "\n")
  
  # CHECKPOINT: Validate marker gene matching
  if (length(matched_genes) == 0) {
    cat("    ERROR: No marker genes matched! Check gene ID format.\n")
    cat("    This will cause MGP to fail.\n")
    # Check if it's a format mismatch
    if (length(all_marker_genes) > 0 && length(expr_rownames_sample) > 0) {
      marker_sample <- head(all_marker_genes, 1)
      expr_sample <- head(rownames(zscore_data_removedBatchEff_cov), 1)
      cat("    Marker gene example:", marker_sample, "\n")
      cat("    Expression gene example:", expr_sample, "\n")
      if (!grepl("^ENSG", marker_sample) && grepl("^ENSG", expr_sample)) {
        cat("    ISSUE: Markers are not Ensembl IDs but expression matrix uses Ensembl IDs!\n")
      } else if (grepl("^ENSG", marker_sample) && !grepl("^ENSG", expr_sample)) {
        cat("    ISSUE: Markers are Ensembl IDs but expression matrix does not use Ensembl IDs!\n")
      }
    }
    stop("ERROR: No marker genes matched. Cannot proceed with MGP estimation. ",
         "Marker genes: ", length(all_marker_genes), ", Expression genes: ", length(genes_only$gene_symbol),
         ". Check that both use Ensembl IDs.")
  } else {
    cat("    First few matched genes:", paste(head(matched_genes, 5), collapse=", "), "\n")
    # Validate matched genes are Ensembl IDs
    if (length(matched_genes) > 0) {
      ensembl_pct <- sum(grepl("^ENSG", matched_genes, ignore.case = TRUE)) / length(matched_genes)
      if (ensembl_pct < 0.8) {
        warning("WARNING: Only ", round(ensembl_pct * 100, 1), "% of matched marker genes are Ensembl IDs")
      } else {
        cat("    ✓ Validated: ", round(ensembl_pct * 100, 1), "% of matched markers are Ensembl IDs\n", sep="")
      }
    }
  }
  
  # Check marker list structure
  cat("  Marker list structure:\n")
  for (i in seq_along(marker_list)) {
    ct_name <- names(marker_list)[i]
    ct_genes <- marker_list[[i]]
    matched_ct_genes <- intersect(ct_genes, genes_only$gene_symbol)
    cat("    ", ct_name, ": ", length(ct_genes), " markers, ", length(matched_ct_genes), " matched\n", sep="")
    if (length(matched_ct_genes) == 0 && length(ct_genes) > 0) {
      cat("      WARNING: No markers matched for this cell type!\n")
      cat("      First marker:", head(ct_genes, 1), "\n")
    }
  }
  
  # Run MGP
  cat("Running MGP deconvolution...\n")
  tryCatch({
    estimations <- markerGeneProfile::mgpEstimate(
      exprData = genes_only,
      genes = marker_list,
      geneColName = 'gene_symbol',
      outlierSampleRemove = FALSE,
      geneTransform = NULL,
      groups = NULL,
      seekConsensus = FALSE,
      removeMinority = FALSE
    )
    cat("  MGP estimation completed successfully\n")
  }, error = function(e) {
    cat("  ERROR in MGP estimation:", e$message, "\n")
    stop("MGP estimation failed: ", e$message)
  })
  
  # MGP returns a list - extract the estimates matrix
  # Based on reference script (rosmap_MGP_estimation.r line 87):
  # rosmap_estimations$estimates %>% as.data.frame()
  # This creates a data.frame with samples as rows and cell types as columns
  cat("  MGP output structure:\n")
  cat("    Top-level names:", names(estimations), "\n")
  cat("    Estimates class:", class(estimations$estimates), "\n")
  
  # Check if simpleScaledEstimation exists (this might be the actual output)
  if ("simpleScaledEstimation" %in% names(estimations)) {
    cat("    Found simpleScaledEstimation\n")
    cat("    simpleScaledEstimation class:", class(estimations$simpleScaledEstimation), "\n")
    if (is.matrix(estimations$simpleScaledEstimation) || is.data.frame(estimations$simpleScaledEstimation)) {
      cat("    simpleScaledEstimation dimensions:", paste(dim(estimations$simpleScaledEstimation), collapse=" x "), "\n")
    }
  }
  
  # Debug: Inspect the structure of estimates
  if (is.list(estimations$estimates)) {
    cat("    Estimates is a list with", length(estimations$estimates), "elements\n")
    if (length(estimations$estimates) > 0) {
      first_el <- estimations$estimates[[1]]
      cat("    First element class:", class(first_el), "\n")
      cat("    First element type:", typeof(first_el), "\n")
      if (is.vector(first_el) || is.matrix(first_el) || is.array(first_el)) {
        cat("    First element dimensions/length:", if(is.vector(first_el)) length(first_el) else paste(dim(first_el), collapse=" x "), "\n")
        if (length(first_el) > 0 && length(first_el) <= 5 && !is.logical(first_el)) {
          cat("    First element values:", paste(head(first_el, 5), collapse=", "), "\n")
        }
      }
      # Check if it's a nested list
      if (is.list(first_el) && length(first_el) > 0) {
        cat("    First element is a nested list with", length(first_el), "elements\n")
        cat("    First nested element class:", class(first_el[[1]]), "\n")
        if (is.matrix(first_el[[1]]) || is.array(first_el[[1]])) {
          cat("    First nested element dimensions:", paste(dim(first_el[[1]]), collapse=" x "), "\n")
        }
      }
    }
  }
  
  # Try to extract estimates - check multiple possible structures
  cell_proportions_df <- NULL
  
  # Method 1: Check if simpleScaledEstimation exists and is usable
  if ("simpleScaledEstimation" %in% names(estimations)) {
    if (is.matrix(estimations$simpleScaledEstimation)) {
      cat("  Using simpleScaledEstimation (matrix)\n")
      cell_proportions_df <- as.data.frame(estimations$simpleScaledEstimation)
      # Check orientation: if more columns than rows, samples are columns
      if (ncol(cell_proportions_df) > nrow(cell_proportions_df)) {
        cat("  Transposing simpleScaledEstimation (samples were columns)\n")
        cell_proportions_df <- as.data.frame(t(estimations$simpleScaledEstimation))
      }
    } else if (is.data.frame(estimations$simpleScaledEstimation)) {
      cat("  Using simpleScaledEstimation (data.frame)\n")
      cell_proportions_df <- estimations$simpleScaledEstimation
      if (ncol(cell_proportions_df) > nrow(cell_proportions_df)) {
        cat("  Transposing simpleScaledEstimation (samples were columns)\n")
        cell_proportions_df <- as.data.frame(t(as.matrix(estimations$simpleScaledEstimation)))
      }
    }
  }
  
  # Method 2: Try to extract from estimates list if it contains nested structures
  if (is.null(cell_proportions_df) && is.list(estimations$estimates)) {
    cat("  Attempting to extract from estimates list...\n")
    # Check if list elements are nested lists with matrices
    first_el <- estimations$estimates[[1]]
    if (is.list(first_el) && length(first_el) > 0) {
      # Nested list structure - try to extract matrices
      sample_names <- NULL
      cell_type_data <- list()
      
      for (i in seq_along(estimations$estimates)) {
        ct_list <- estimations$estimates[[i]]
        ct_name <- names(estimations$estimates)[i]
        if (is.null(ct_name)) ct_name <- paste0("CellType_", i)
        
        # Look for a matrix in the nested list
        if (is.list(ct_list)) {
          for (j in seq_along(ct_list)) {
            if (is.matrix(ct_list[[j]]) && nrow(ct_list[[j]]) > 0) {
              # Extract first row (estimates)
              if (is.null(sample_names)) {
                sample_names <- colnames(ct_list[[j]])
                if (is.null(sample_names)) {
                  sample_names <- paste0("Sample_", seq_len(ncol(ct_list[[j]])))
                }
              }
              cell_type_data[[ct_name]] <- ct_list[[j]][1, , drop = TRUE]
              break
            }
          }
        } else if (is.matrix(ct_list) && nrow(ct_list) > 0) {
          if (is.null(sample_names)) {
            sample_names <- colnames(ct_list)
            if (is.null(sample_names)) {
              sample_names <- paste0("Sample_", seq_len(ncol(ct_list)))
            }
          }
          cell_type_data[[ct_name]] <- ct_list[1, , drop = TRUE]
        } else if (is.vector(ct_list) && length(ct_list) > 1 && !is.logical(ct_list)) {
          if (is.null(sample_names)) {
            sample_names <- names(ct_list)
            if (is.null(sample_names)) {
              sample_names <- paste0("Sample_", seq_along(ct_list))
            }
          }
          cell_type_data[[ct_name]] <- ct_list
        }
      }
      
      if (length(cell_type_data) > 0) {
        cat("  Extracted", length(cell_type_data), "cell types with", length(cell_type_data[[1]]), "samples each\n")
        cell_proportions_df <- as.data.frame(cell_type_data)
        if (!is.null(sample_names) && length(sample_names) == nrow(cell_proportions_df)) {
          rownames(cell_proportions_df) <- sample_names
        }
      }
    }
  }
  
  # Method 3: Try direct conversion (as in reference script)
  if (is.null(cell_proportions_df)) {
    cat("  Attempting direct conversion as in reference script...\n")
    cell_proportions_df <- as.data.frame(estimations$estimates)
    cat("  After as.data.frame() - dimensions:", paste(dim(cell_proportions_df), collapse=" x "), "\n")
    
    # If we got 1 row, the structure is wrong
    if (nrow(cell_proportions_df) == 1 && ncol(cell_proportions_df) == length(estimations$estimates)) {
      cat("  ERROR: Direct conversion failed - only 1 row detected\n")
      cat("  This suggests MGP estimation may have failed or structure is unexpected\n")
      stop("Cannot extract cell proportions from MGP output. The estimates list contains logical NAs instead of numeric vectors. Check MGP estimation input and marker gene matching.")
    }
  }
  
  cat("  Final data.frame dimensions:", paste(dim(cell_proportions_df), collapse=" x "), "\n")
  cat("  Number of samples (rows):", nrow(cell_proportions_df), "\n")
  cat("  Number of cell types (cols):", ncol(cell_proportions_df), "\n")
  cat("  Cell type names:", paste(head(colnames(cell_proportions_df), 5), collapse=", "), "...\n")
  
  # Convert to matrix (already in correct orientation: samples as rows, cell types as columns)
  cell_proportions <- as.matrix(cell_proportions_df)
  
  cat("  Final cell_proportions dimensions:", dim(cell_proportions), "\n")
  cat("  Final cell_proportions class:", class(cell_proportions), "\n")
  cat("  Samples (rows):", nrow(cell_proportions), "\n")
  cat("  Cell types (cols):", ncol(cell_proportions), "\n")
  
} else if (opt$deconv_tool == "dtangle") {
  # dtangle deconvolution
  if (!require("dtangle", quietly = TRUE)) {
    stop("dtangle package is required for dtangle deconvolution. Please install it.")
  }
  
  # dtangle requires reference expression matrix - this is a simplified version
  # In practice, you'd need a reference single-cell dataset
  stop("dtangle requires a reference expression matrix. Please provide reference data or use MGP.")
  
} else if (opt$deconv_tool == "CIBERSORT") {
  # CIBERSORT deconvolution
  stop("CIBERSORT requires external CIBERSORT software. Please run CIBERSORT separately or use MGP.")
  
} else if (opt$deconv_tool == "MuSiC") {
  # MuSiC deconvolution
  if (!require("MuSiC", quietly = TRUE)) {
    stop("MuSiC package is required for MuSiC deconvolution. Please install it.")
  }
  
  stop("MuSiC requires reference single-cell data. Please provide reference data or use MGP.")
  
} else if (opt$deconv_tool == "Bisque") {
  # Bisque deconvolution
  if (!require("BisqueRNA", quietly = TRUE)) {
    stop("BisqueRNA package is required for Bisque deconvolution. Please install it.")
  }
  
  stop("Bisque requires reference single-cell data. Please provide reference data or use MGP.")
  
} else if (opt$deconv_tool == "DWLS") {
  # DWLS deconvolution
  stop("DWLS requires reference single-cell data and additional setup. Please use MGP or provide reference data.")
  
} else {
  stop(paste("Unknown deconvolution tool:", opt$deconv_tool))
}

# ============================================================================
# Process and save results
# ============================================================================
cat("Processing deconvolution results...\n")
cat("  Input cell_proportions dimensions:", dim(cell_proportions), "\n")
cat("  Input cell_proportions class:", class(cell_proportions), "\n")

# Convert to data frame - ensure samples are rows and cell types are columns
if (is.matrix(cell_proportions)) {
  # Transpose if needed: MGP typically returns samples as columns, cell types as rows
  # But we want samples as rows for the output
  if (ncol(cell_proportions) < nrow(cell_proportions)) {
    # Likely samples are rows already
    cell_proportions_df <- as.data.frame(cell_proportions)
  } else {
    # Likely samples are columns, need to transpose
    cat("  Transposing cell_proportions (samples were columns)\n")
    cell_proportions_df <- as.data.frame(t(cell_proportions))
  }
} else {
  cell_proportions_df <- as.data.frame(cell_proportions)
}

cat("  After conversion - cell_proportions_df dimensions:", dim(cell_proportions_df), "\n")
cat("  After conversion - cell_proportions_df columns:", paste(head(colnames(cell_proportions_df), 5), collapse=", "), "\n")

# Determine sample ID column
specimen_col <- find_column(cell_proportions_df, c("specimenID", "SpecimenID", "synapseID", "SynapseID", 
                                                    "sampleID", "SampleID", "id", "ID"))
if (is.null(specimen_col)) {
  cell_proportions_df <- cell_proportions_df %>% tibble::rownames_to_column(var = "specimenID")
  specimen_col <- "specimenID"
}

# Save unscaled proportions (save to current directory for Nextflow, also copy to output_dir)
write.csv(cell_proportions_df,
          opt$proportions_output,
          row.names = FALSE)

# Scale proportions
cell_proportions_scaled <- cell_proportions_df
numeric_cols <- setdiff(colnames(cell_proportions_scaled), specimen_col)
cell_proportions_scaled[, numeric_cols] <- scale(cell_proportions_scaled[, numeric_cols, drop = FALSE])

# Save scaled proportions
write.csv(cell_proportions_scaled,
          opt$proportions_scaled_output,
          row.names = FALSE)

# Also copy to output_dir if different from current directory
if (opt$output_dir != "." && opt$output_dir != getwd()) {
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(opt$proportions_output, file.path(opt$output_dir, opt$proportions_output), overwrite = TRUE)
  file.copy(opt$proportions_scaled_output, file.path(opt$output_dir, opt$proportions_scaled_output), overwrite = TRUE)
}

# Generate summary
summary_text <- paste0(
  "Cell Type Deconvolution Summary\n",
  "==============================\n",
  "Deconvolution tool: ", opt$deconv_tool, "\n",
  "Reference taxonomy: ", opt$reference_taxonomy, "\n",
  "Number of samples: ", nrow(cell_proportions_df), "\n",
  "Number of cell types: ", length(cell_types), "\n",
  "Cell types: ", paste(cell_types, collapse = ", "), "\n",
  "Marker file: ", ifelse(is.null(opt$marker_file), "default", opt$marker_file), "\n"
)

cat("\n", summary_text, "\n")
writeLines(summary_text, file.path(opt$output_dir, "deconv_summary.txt"))

cat("Cell type deconvolution completed!\n")
cat("Output files saved to:", opt$output_dir, "\n")

