#!/usr/bin/env Rscript
# One-time setup: Build GVEX RNA_ID <-> WGS_ID mapping from SYNAPSE_METADATA_MANIFEST
# and BrainGVEX imputed VCF sample IDs. Output used for phenotype/genotype sample matching.
# Run from project root or set paths below.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Paths (adjust if not run from nextflow project root)
project_dir <- Sys.getenv("PROJECT_DIR", "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
manifest_path <- "/external/rprshnas01/external_data/psychencode/PsychENCODE/BrainGVEX/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv"
vcf_path <- "/external/rprshnas01/external_data/psychencode/PsychENCODE/genotypes_BrainGVEX/DNA/BrainGVEX.psy.chr1.dose.vcf.gz"
bcftools_path <- "/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools"
out_dir <- file.path(project_dir, "data_input", "gvex")
out_file <- file.path(out_dir, "gvex_rna_wgs_id_mapping.tsv")

# RNA IDs from manifest (geneExpression only)
manifest <- read_tsv(manifest_path, show_col_types = FALSE)
gvex_rna_ids <- manifest %>%
  filter(dataType == "geneExpression") %>%
  distinct(individualID) %>%
  pull(individualID)

# WGS sample IDs from VCF (via bcftools)
stopifnot("bcftools not found" = file.exists(bcftools_path))
vcf_ids <- system2(bcftools_path, args = c("query", "-l", vcf_path), stdout = TRUE)

# Overlap = subjects with both RNA and WGS
overlap_ids <- intersect(gvex_rna_ids, vcf_ids)

mapping <- tibble(
  RNA_ID = overlap_ids,
  WGS_ID = overlap_ids
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(mapping, out_file)

message("Wrote ", nrow(mapping), " rows to ", out_file)
