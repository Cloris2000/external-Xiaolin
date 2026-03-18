#!/usr/bin/env Rscript
# One-time prep for NIMH HBCC:
# 1. Subset CMC count matrix to HBCC samples only (columns CMC_HBCC_RNA_PFC_*).
# 2. Build HBCC metadata from SYNAPSE_TABLE_QUERY: specimenID <-> individualID.
# Output: HBCC_count_matrix.csv, HBCC_metadata.csv in data_input/nimh_hbcc/
#
# Mapping chain for GWAS:
#   Count matrix columns = specimenID (e.g. CMC_HBCC_RNA_PFC_2942)
#   SYNAPSE_TABLE_QUERY: specimenID -> individualID (e.g. CMC_HBCC_010)
#   CMC_Human_SNP_metadata: Individual_ID -> Genotyping_Sample_ID (e.g. 4040296003_A)
#   Normalized VCF/PSAM: sample ID = 0_Genotyping_Sample_ID (pipeline adds 0_ for CMC)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

project_dir <- Sys.getenv("PROJECT_DIR", "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
out_dir <- file.path(project_dir, "data_input", "nimh_hbcc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Paths
cmc_count_path <- "/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/CMC_count_matrix.csv"
synapse_table_path <- "/external/rprshnas01/external_data/psychencode/PsychENCODE/CMC/Metadata/SYNAPSE_TABLE_QUERY_123020650.csv"

# ---- 1. Subset count matrix to HBCC only ----
cat("Reading CMC count matrix...\n")
cmc <- read_csv(cmc_count_path, show_col_types = FALSE)
# First column is gene ID (empty name or "gene")
gene_col <- names(cmc)[1]
hbcc_cols <- names(cmc)[grepl("^CMC_HBCC_RNA_PFC_", names(cmc))]
cat("  HBCC columns found:", length(hbcc_cols), "\n")

hbcc_count <- cmc %>% select(all_of(c(gene_col, hbcc_cols)))
# Remove Ensembl version suffix from gene IDs (e.g. ENSG00000000003.10 -> ENSG00000000003)
hbcc_count[[gene_col]] <- sub("\\.[0-9]+$", "", hbcc_count[[gene_col]])
cat("  Stripped version numbers from gene IDs (e.g. ENSG...10 -> ENSG...)\n")
write_csv(hbcc_count, file.path(out_dir, "HBCC_count_matrix.csv"))
cat("  Wrote", nrow(hbcc_count), "genes x", ncol(hbcc_count) - 1, "samples ->", file.path(out_dir, "HBCC_count_matrix.csv"), "\n")

# ---- 2. Build HBCC metadata (specimenID, individualID) from Synapse table ----
cat("Reading Synapse table...\n")
syn <- read_csv(synapse_table_path, show_col_types = FALSE)
# Filter: rnaSeq, geneExpression, and specimenID like CMC_HBCC_RNA_PFC_*
hbcc_rna <- syn %>%
  filter(assay == "rnaSeq", dataType == "geneExpression") %>%
  filter(grepl("^CMC_HBCC_RNA_PFC_", specimenID))
cat("  HBCC RNA rows:", nrow(hbcc_rna), "\n")

# Keep one row per specimen; include tech-cov-like columns (similar to ROSMAP):
# platform, pH, libraryPrep, readLength, runType, BrodmannArea (from Synapse)
meta_cols <- c("individualID", "specimenID", "tissue", "PMI", "RIN", "organ", "hemisphere")
# Add columns that exist and are useful as tech covs (like ROSMAP)
for (c in c("platform", "pH", "libraryPrep", "readLength", "runType", "BrodmannArea")) {
  if (c %in% names(hbcc_rna)) meta_cols <- c(meta_cols, c)
}
meta <- hbcc_rna %>%
  select(any_of(meta_cols)) %>%
  distinct(specimenID, .keep_all = TRUE) %>%
  filter(specimenID %in% hbcc_cols)
cat("  Metadata rows (in count matrix):", nrow(meta), "\n")

# Join Capstone clinical for sex/age (required by pheno_cov_prep)
capstone_clinical_path <- "/external/rprshnas01/external_data/psychencode/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"
if (file.exists(capstone_clinical_path)) {
  clinical <- read_csv(capstone_clinical_path, show_col_types = FALSE)
  clinical_sub <- clinical %>% select(individualID, reportedGender, ageDeath, ethnicity, primaryDiagnosis)
  meta <- meta %>% left_join(clinical_sub, by = "individualID")
  cat("  Joined Capstone clinical; non-NA ageDeath:", sum(!is.na(meta$ageDeath)), ", reportedGender:", sum(!is.na(meta$reportedGender)), "\n")
}

write_csv(meta, file.path(out_dir, "HBCC_metadata.csv"))
cat("  Wrote", file.path(out_dir, "HBCC_metadata.csv"), "\n")

# ---- 3. Optional: write RNA <-> WGS ID mapping for reference ----
biospec_path <- "/external/rprshnas01/netdata_kcni/stlab/CMC_genotypes/SNPs/Release3/Metadata/CMC_Human_SNP_metadata.csv"
if (file.exists(biospec_path)) {
  biospec <- read_csv(biospec_path, show_col_types = FALSE)
  biospec_hbcc <- biospec %>% filter(Study == "CMC_HBCC") %>% select(Individual_ID, Genotyping_Sample_ID)
  # Join: individualID (metadata) -> Individual_ID -> Genotyping_Sample_ID
  mapping <- meta %>% select(RNA_specimenID = specimenID, individualID) %>%
    left_join(biospec_hbcc, by = c("individualID" = "Individual_ID")) %>%
    rename(WGS_ID = Genotyping_Sample_ID)
  write_csv(mapping, file.path(out_dir, "HBCC_rna_wgs_id_mapping.csv"))
  cat("  Wrote RNA<->WGS mapping (", sum(!is.na(mapping$WGS_ID)), " with genotype)\n", sep = "")
}

cat("Done. Use col_sample_id_for_matching = 'specimenID' and biospecimen_file = CMC_Human_SNP_metadata.csv for HBCC.\n")
