#!/usr/bin/env Rscript
# =============================================================================
# topmed_samples_not_in_pipeline.R
#
# Identifies ROSMAP TOPmed imputed samples that are NOT present in the
# current WGS-based Nextflow pipeline analysis, then checks which of those
# have matching bulk RNA-seq data — both across all tissues and for DLPFC
# specifically — using the actual count matrix sample list as the source of
# truth for what RNA-seq samples exist.
#
# ID bridging logic:
#   TOPmed n1686 (Illumina): MAP/ROS##### IDs  -> clinical Study+projid -> individualID
#   TOPmed n381 (Affymetrix): 11AD##### IDs    -> biospecimen snpArray  -> individualID
#   Pipeline WGS: SM-##### IDs                 -> biospecimen WGS       -> individualID
#   Count matrix: syn##### column names         -> combined_metrics synapseID -> individualID
#
# Outputs (written to OUTPUT_DIR):
#   topmed_not_in_pipeline.txt              — all TOPmed-only individuals
#   topmed_not_in_pipeline_rnaseq_any.txt   — subset with RNA-seq (any tissue)
#   topmed_not_in_pipeline_rnaseq_dlpfc.txt — subset with DLPFC RNA-seq
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# =============================================================================
# PATHS — edit if files move
# =============================================================================

TOPMED_FAM        <- "/external/rprshnas01/external_data/rosmap/genotype/TOPmed_imputed/vcf/merged/merged_overlap.fam"
PIPELINE_PSAM     <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/ROSMAP.QC.1.psam"
BIOSPECIMEN_FILE  <- "/nethome/kcni/xzhou/GWAS_tut/AMP-AD/ROSMAP_biospecimen_metadata.csv"
CLINICAL_FILE     <- "/external/rprshnas01/external_data/rosmap/metadata/ROSMAP_clinical.csv"
COMBINED_METRICS  <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_combined_metrics.csv"
COUNT_MATRIX      <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_batch_all.csv"

OUTPUT_DIR        <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP"

# =============================================================================
# 1. Load ID mapping tables
# =============================================================================

message("Loading biospecimen metadata...")
bio <- fread(BIOSPECIMEN_FILE, data.table = FALSE)

# WGS: specimenID -> individualID  (for SM-##### pipeline IDs)
wgs_map <- bio %>%
  filter(assay == "wholeGenomeSeq") %>%
  select(specimenID, individualID) %>%
  distinct()

# snpArray: specimenID -> individualID  (for 11AD##### TOPmed IDs)
snp_map <- bio %>%
  filter(assay == "snpArray") %>%
  select(specimenID, individualID) %>%
  distinct()

message("Loading clinical metadata...")
# colClasses ensures projid is read as character so leading zeros are preserved
# (e.g. projid "00285563" stays "00285563", not 285563). This is critical
# because TOPmed n1686 IDs use zero-padded 8-digit projids (MAP00285563).
clinical <- fread(CLINICAL_FILE, data.table = FALSE,
                  colClasses = list(character = "projid"))

# Study+projid -> individualID  (for MAP/ROS##### TOPmed and pipeline IDs)
clinical_map <- clinical %>%
  filter(!is.na(projid), !is.na(Study), !is.na(individualID)) %>%
  mutate(topmed_id = paste0(Study, projid)) %>%
  select(topmed_id, individualID) %>%
  distinct()

# =============================================================================
# 2. Load pipeline QC samples and resolve to individualID
# =============================================================================

message("Loading pipeline QC samples (PSAM)...")
psam <- fread(PIPELINE_PSAM, data.table = FALSE)
colnames(psam)[1] <- "FID"
pipeline_spec_ids <- psam$IID

# SM-##### IDs: WGS biospecimen;  MAP/ROS##### IDs: clinical map
pipeline_indivs <- bind_rows(
  data.frame(specimenID = pipeline_spec_ids) %>%
    inner_join(wgs_map, by = "specimenID") %>%
    select(specimenID, individualID),
  data.frame(topmed_id = pipeline_spec_ids) %>%
    inner_join(clinical_map, by = "topmed_id") %>%
    rename(specimenID = topmed_id) %>%
    select(specimenID, individualID)
) %>% distinct()

pipeline_indiv_set <- unique(pipeline_indivs$individualID)

message(sprintf("Pipeline QC samples: %d  |  resolved to %d unique individuals",
                length(pipeline_spec_ids), length(pipeline_indiv_set)))

# =============================================================================
# 3. Load TOPmed merged samples and resolve to individualID
# =============================================================================

message("Loading TOPmed merged FAM file...")
fam <- fread(TOPMED_FAM, header = FALSE, data.table = FALSE,
             col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
topmed_spec_ids <- fam$IID

# MAP/ROS##### -> clinical_map;  11AD##### -> snp_map
topmed_indivs <- bind_rows(
  data.frame(topmed_id = topmed_spec_ids) %>%
    inner_join(clinical_map, by = "topmed_id") %>%
    rename(specimenID = topmed_id) %>%
    select(specimenID, individualID),
  data.frame(specimenID = topmed_spec_ids) %>%
    inner_join(snp_map, by = "specimenID") %>%
    select(specimenID, individualID)
) %>%
  distinct() %>%
  group_by(specimenID) %>%
  slice(1) %>%
  ungroup()

no_link <- setdiff(topmed_spec_ids, topmed_indivs$specimenID)
message(sprintf("TOPmed merged samples: %d  |  resolved to %d individuals  |  no ID link: %d",
                length(topmed_spec_ids), n_distinct(topmed_indivs$individualID), length(no_link)))
if (length(no_link) > 0)
  message("  No-link samples: ", paste(no_link, collapse = ", "))

# =============================================================================
# 4. Identify TOPmed samples NOT in pipeline
# =============================================================================

topmed_not_in_pipeline <- topmed_indivs %>%
  filter(!individualID %in% pipeline_indiv_set)

n_illumina <- sum(grepl("^(MAP|ROS)", topmed_not_in_pipeline$specimenID))
n_affy     <- sum(grepl("^11AD",      topmed_not_in_pipeline$specimenID))

message(sprintf("\n=== TOPmed samples NOT in pipeline ==="))
message(sprintf("  Total: %d individuals (%d Illumina MAP/ROS, %d Affymetrix 11AD)",
                nrow(topmed_not_in_pipeline), n_illumina, n_affy))

# =============================================================================
# 5. Build RNA-seq sample table from count matrix + combined_metrics
#
#    Count matrix columns are syn##### IDs.
#    combined_metrics maps synapseID -> individualID + tissue info.
#    This gives us the definitive list of samples with actual expression data.
# =============================================================================

message("\nReading count matrix column names (RNA-seq samples)...")
# Read only the header row — fast even for large files
count_header <- scan(COUNT_MATRIX, what = character(), nlines = 1,
                     sep = ",", quiet = TRUE)
# First two entries are "" and "gene_id"; the rest are syn##### sample IDs
count_syn_ids <- count_header[-(1:2)]
count_syn_ids <- gsub('"', '', count_syn_ids)  # strip any surrounding quotes
message(sprintf("  Count matrix samples (syn IDs): %d", length(count_syn_ids)))

message("Loading combined metrics for RNA-seq metadata...")
metrics <- fread(COMBINED_METRICS, data.table = FALSE)

# Keep only rows whose synapseID appears in the count matrix
rnaseq_meta <- metrics %>%
  filter(assay == "rnaSeq", synapseID %in% count_syn_ids, !is.na(individualID)) %>%
  select(synapseID, specimenID, individualID, tissue, exclude) %>%
  distinct()

message(sprintf("  RNA-seq entries matched to count matrix: %d  (%d unique individuals)",
                nrow(rnaseq_meta), n_distinct(rnaseq_meta$individualID)))

# Tissue breakdown across all samples (for reference)
message("  Tissue breakdown in count matrix RNA-seq:")
rnaseq_meta %>%
  count(tissue, name = "n_samples") %>%
  arrange(desc(n_samples)) %>%
  { for (i in seq_len(nrow(.))) message(sprintf("    %4d  %s", .[[i,"n_samples"]], .[[i,"tissue"]])) }

# Any tissue
rnaseq_any <- rnaseq_meta %>%
  rename(rnaseq_synapseID  = synapseID,
         rnaseq_specimenID = specimenID,
         rnaseq_tissue     = tissue,
         rnaseq_exclude    = exclude)

# DLPFC only
rnaseq_dlpfc <- rnaseq_meta %>%
  filter(tolower(tissue) == "dorsolateral prefrontal cortex") %>%
  rename(rnaseq_synapseID  = synapseID,
         rnaseq_specimenID = specimenID,
         rnaseq_tissue     = tissue,
         rnaseq_exclude    = exclude)

# =============================================================================
# 6. Match TOPmed-not-in-pipeline individuals to RNA-seq
# =============================================================================

not_in_pipeline_set <- unique(topmed_not_in_pipeline$individualID)

topmed_with_rnaseq_any <- topmed_not_in_pipeline %>%
  inner_join(rnaseq_any, by = "individualID") %>%
  rename(topmed_specimenID = specimenID)

topmed_with_rnaseq_dlpfc <- topmed_not_in_pipeline %>%
  inner_join(rnaseq_dlpfc, by = "individualID") %>%
  rename(topmed_specimenID = specimenID)

n_any_indivs   <- n_distinct(topmed_with_rnaseq_any$individualID)
n_dlpfc_indivs <- n_distinct(topmed_with_rnaseq_dlpfc$individualID)

message(sprintf("\n=== RNA-seq matching results ==="))
message(sprintf("  With bulk RNA-seq in count matrix (any tissue):  %d unique individuals", n_any_indivs))
message(sprintf("  With bulk RNA-seq in count matrix (DLPFC only):  %d unique individuals", n_dlpfc_indivs))
message(sprintf("  Without any RNA-seq match:                       %d unique individuals",
                length(not_in_pipeline_set) - n_any_indivs))

tissue_summary <- topmed_with_rnaseq_any %>%
  count(rnaseq_tissue, name = "n_rna_samples") %>%
  arrange(desc(n_rna_samples))
message("\n  RNA-seq tissue breakdown for matched individuals (samples, not individuals):")
for (i in seq_len(nrow(tissue_summary)))
  message(sprintf("    %4d  %s", tissue_summary[[i, "n_rna_samples"]], tissue_summary[[i, "rnaseq_tissue"]]))

# =============================================================================
# 7. Write outputs
# =============================================================================

message("\nWriting output files...")

out1 <- topmed_not_in_pipeline %>%
  rename(topmed_specimenID = specimenID) %>%
  arrange(topmed_specimenID)

out2 <- topmed_with_rnaseq_any %>%
  arrange(topmed_specimenID, rnaseq_specimenID)

out3 <- topmed_with_rnaseq_dlpfc %>%
  arrange(topmed_specimenID, rnaseq_specimenID)

fwrite(out1, file.path(OUTPUT_DIR, "topmed_not_in_pipeline.txt"),
       sep = "\t", quote = FALSE, na = "NA")
fwrite(out2, file.path(OUTPUT_DIR, "topmed_not_in_pipeline_rnaseq_any.txt"),
       sep = "\t", quote = FALSE, na = "NA")
fwrite(out3, file.path(OUTPUT_DIR, "topmed_not_in_pipeline_rnaseq_dlpfc.txt"),
       sep = "\t", quote = FALSE, na = "NA")

message(sprintf("\nDone. Output written to: %s", OUTPUT_DIR))
message(sprintf("  topmed_not_in_pipeline.txt              — %d rows (%d unique individuals)",
                nrow(out1), n_distinct(out1$individualID)))
message(sprintf("  topmed_not_in_pipeline_rnaseq_any.txt   — %d rows (%d unique individuals)",
                nrow(out2), n_distinct(out2$individualID)))
message(sprintf("  topmed_not_in_pipeline_rnaseq_dlpfc.txt — %d rows (%d unique individuals)",
                nrow(out3), n_distinct(out3$individualID)))
