#!/usr/bin/env Rscript
# Cell Type Proportion Estimation using MGP (Marker Gene Profile)
# Estimates cell type proportions from corrected expression data

suppressPackageStartupMessages({
  library(optparse)
  library(GEOquery)
  library(dtangle)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(limma)
  library(ggplot2)
  library(reshape2)
  library(edgeR)
  library(matrixStats)
  library(dplyr)
  library(tidyverse)
  library(broom)
  library(knitr)
  library(ggpubr)
  library(biomaRt)
  library(ggrepel)
  library(patchwork)
  library(ggsignif)
  library(modelr)
  library(cowplot)
  library(gridExtra)
  library(RColorBrewer)
  library(markerGeneProfile)
  library(data.table)
  library(readr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--corrected_data"), type="character", default=NULL,
              help="Path to corrected z-score data RData file"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Path to cleaned metadata CSV file"),
  make_option(c("--deconv_tool"), type="character", default="MGP",
              help="Deconvolution tool to use (default: MGP)"),
  make_option(c("--reference_taxonomy"), type="character", default="sonny_markers",
              help="Reference taxonomy to use (default: sonny_markers)"),
  make_option(c("--marker_file"), type="character", default=NULL,
              help="Path or URL to marker file"),
  make_option(c("--hgnc_mapping_file"), type="character", default=NULL,
              help="Path to HGNC mapping file"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--proportions_output"), type="character", default="cell_proportions.csv",
              help="Output filename for cell proportions"),
  make_option(c("--proportions_scaled_output"), type="character", default="cell_proportions_scaled.csv",
              help="Output filename for scaled cell proportions")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Loading corrected data and metadata...\n")
# Load corrected data
zscore_data_removedBatchEff_cov <- readRDS(opt$corrected_data)

# Load metadata if provided
metadata <- NULL
if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
  metadata <- read.csv(opt$metadata)
}

cat("Loading marker gene files...\n")
# Load marker genes based on reference taxonomy
if (opt$reference_taxonomy == "sonny_markers") {
  if (is.null(opt$marker_file)) {
    marker_url <- 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'
    cat("Downloading markers from:", marker_url, "\n")
    sonny_markers <- read_csv(marker_url, show_col_types = FALSE)
  } else {
    if (startsWith(opt$marker_file, "http")) {
      sonny_markers <- read_csv(opt$marker_file, show_col_types = FALSE)
    } else {
      sonny_markers <- read_csv(opt$marker_file, show_col_types = FALSE)
    }
  }
  
  colnames(sonny_markers) <- colnames(sonny_markers) %>% make.names() %>% tolower()
  
  # Load HGNC mapping
  if (!is.null(opt$hgnc_mapping_file) && file.exists(opt$hgnc_mapping_file)) {
    hgnc_mapping <- read_tsv(opt$hgnc_mapping_file, show_col_types = FALSE)
    
    # Merge markers with HGNC mapping
    sonny_hgnc_merged_markers <- left_join(
      sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id),
      hgnc_mapping %>% distinct(entrez_id, .keep_all = TRUE) %>%
        dplyr::select(entrez_id, ensembl_gene_id) %>%
        dplyr::rename(ensembl_id = ensembl_gene_id),
      by = "entrez_id"
    ) %>%
      dplyr::select(gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>%
      group_by(subclass) %>%
      arrange(subclass, -average.log.fold.change) %>%
      ungroup()
  } else {
    cat("Warning: HGNC mapping file not found. Using markers without Ensembl ID mapping.\n")
    sonny_hgnc_merged_markers <- sonny_markers
  }
} else {
  stop("Only 'sonny_markers' reference taxonomy is currently supported.")
}

cat("Preparing marker gene list...\n")
# Get markers for MGP
new_markers <- sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types <- new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique()

# Create marker list for each cell type
new_marker_list <- lapply(new_cell_types, function(cell_type) {
  if ("ensembl_id" %in% colnames(new_markers)) {
    cell_type_marker_list <- new_markers %>%
      filter(subclass == cell_type, 
             ensembl_id %in% rownames(zscore_data_removedBatchEff_cov)) %>%
      pull(ensembl_id)
  } else {
    # Fallback to gene names if no Ensembl IDs
    cell_type_marker_list <- new_markers %>%
      filter(subclass == cell_type,
             gene %in% rownames(zscore_data_removedBatchEff_cov)) %>%
      pull(gene)
  }
  return(cell_type_marker_list)
})
names(new_marker_list) <- new_cell_types

cat("Cell types to estimate:", paste(new_cell_types, collapse = ", "), "\n")

# Prepare expression data
cat("Preparing expression data for MGP estimation...\n")
matrix <- as.data.frame(zscore_data_removedBatchEff_cov)

# Determine gene identifier column name
if ("ensembl_id" %in% colnames(new_markers)) {
  gene_col_name <- ifelse(all(rownames(matrix) %in% new_markers$ensembl_id), 
                          "gene_symbol", "gene_symbol")
  matrix <- matrix %>% rownames_to_column(var = "gene_symbol")
} else {
  matrix <- matrix %>% rownames_to_column(var = "gene_symbol")
}

# Remove rows with empty gene symbols
genes_only <- matrix %>% subset(gene_symbol != "")

# Remove duplicates
if (length(which(duplicated(genes_only$gene_symbol))) != 0) {
  genes_only <- genes_only[-which(duplicated(genes_only$gene_symbol)), ]
}

cat("Running MGP estimation...\n")
# Run MGP estimation
rosmap_estimations <- mgpEstimate(
  exprData = genes_only,
  genes = new_marker_list,
  geneColName = 'gene_symbol',
  outlierSampleRemove = FALSE,
  geneTransform = NULL,
  groups = NULL,
  seekConsensus = FALSE,
  removeMinority = FALSE
)

cat("Processing MGP results...\n")
# Coerce estimations to data frame and scale
rosmap_estimations_scaled <- rosmap_estimations$estimates %>%
  as.data.frame() %>%
  scale() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "specimenID")

# Save scaled proportions
write.csv(rosmap_estimations_scaled, 
          file.path(opt$output_dir, opt$proportions_scaled_output), 
          row.names = FALSE)

# Also save unscaled proportions
rosmap_estimations_unscaled <- rosmap_estimations$estimates %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "specimenID")

write.csv(rosmap_estimations_unscaled,
          file.path(opt$output_dir, opt$proportions_output),
          row.names = FALSE)

# Generate summary
summary_text <- paste0(
  "MGP Estimation Summary\n",
  "=====================\n",
  "Number of samples: ", nrow(rosmap_estimations_scaled), "\n",
  "Number of cell types: ", length(new_cell_types), "\n",
  "Cell types: ", paste(new_cell_types, collapse = ", "), "\n",
  "Reference taxonomy: ", opt$reference_taxonomy, "\n",
  "Deconvolution tool: ", opt$deconv_tool, "\n"
)

cat(summary_text)
writeLines(summary_text, file.path(opt$output_dir, "deconv_summary.txt"))

cat("MGP estimation completed!\n")
cat("Output files saved to:", opt$output_dir, "\n")



