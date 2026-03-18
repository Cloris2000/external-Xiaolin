#!/usr/bin/env Rscript
# =============================================================================
# prep_rosmap_array_samples.R
#
# Maps the 170 ROSMAP TOPmed DLPFC samples to their Study+projid FIDs so
# that fid_method = 'study_projid' works in the pipeline without changes to
# pheno_cov_prep.R.
#
# Logic:
#   MAP/ROS##### (n1686 Illumina) -> directly matches Study+projid in clinical
#   11AD#####    (n381 Affymetrix) -> biospecimen snpArray -> individualID
#                                  -> clinical Study+projid
#
# Output: results/ROSMAP_array/samples_to_keep.txt
#   Two-column tab-separated: FID  IID  (both = Study+projid, e.g. MAP12345678)
#   Also includes the original topmed_specimenID as a comment column for tracing.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

DLPFC_FILE       <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/topmed_not_in_pipeline_rnaseq_dlpfc.txt"
CLINICAL_FILE    <- "/external/rprshnas01/external_data/rosmap/metadata/ROSMAP_clinical.csv"
BIOSPECIMEN_FILE <- "/nethome/kcni/xzhou/GWAS_tut/AMP-AD/ROSMAP_biospecimen_metadata.csv"
OUTPUT_DIR       <- "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP_array"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Load clinical: Study+projid -> individualID ----
message("Loading clinical data...")
clinical <- fread(CLINICAL_FILE, data.table = FALSE,
                  colClasses = list(character = "projid"))
clinical_map <- clinical %>%
  filter(!is.na(projid), !is.na(Study), !is.na(individualID)) %>%
  mutate(fid = paste0(Study, projid)) %>%
  select(individualID, fid) %>%
  distinct()

# Also build reverse: fid -> individualID (for MAP/ROS IDs that ARE the FID)
fid_to_indiv <- setNames(clinical_map$individualID, clinical_map$fid)
indiv_to_fid <- setNames(clinical_map$fid, clinical_map$individualID)

# ---- Load biospecimen: snpArray specimenID -> individualID ----
message("Loading biospecimen data...")
bio <- fread(BIOSPECIMEN_FILE, data.table = FALSE)
snp_spec_to_indiv <- bio %>%
  filter(assay == "snpArray") %>%
  select(specimenID, individualID) %>%
  distinct() %>%
  { setNames(.$individualID, .$specimenID) }

# ---- Load 170 DLPFC samples ----
message("Loading DLPFC sample list...")
dlpfc <- fread(DLPFC_FILE, data.table = FALSE, sep = "\t")
message(sprintf("  Input: %d rows, %d unique topmed_specimenIDs",
                nrow(dlpfc), n_distinct(dlpfc$topmed_specimenID)))

# ---- Map each sample to its Study+projid FID ----
result <- dlpfc %>%
  select(topmed_specimenID, individualID) %>%
  distinct() %>%
  rowwise() %>%
  mutate(
    fid = {
      spec   <- topmed_specimenID
      indiv  <- individualID
      # MAP/ROS##### IDs: directly look up as fid in clinical_map
      if (grepl("^(MAP|ROS)", spec) && spec %in% names(fid_to_indiv)) {
        spec  # the topmed_specimenID IS the FID (Study+projid)
      } else if (grepl("^(MAP|ROS)", spec)) {
        # Try via individualID route
        indiv_to_fid[indiv] %||% NA_character_
      } else {
        # 11AD##### IDs: go via biospecimen -> individualID -> FID
        bio_indiv <- snp_spec_to_indiv[spec]
        if (!is.na(bio_indiv) && bio_indiv %in% names(indiv_to_fid)) {
          indiv_to_fid[bio_indiv]
        } else {
          # Fall back to individualID from dlpfc file
          indiv_to_fid[indiv] %||% NA_character_
        }
      }
    }
  ) %>%
  ungroup()

# ---- Report mapping results ----
n_mapped   <- sum(!is.na(result$fid))
n_unmapped <- sum(is.na(result$fid))
message(sprintf("\nMapping results:"))
message(sprintf("  Total unique genotype IDs: %d", nrow(result)))
message(sprintf("  Successfully mapped to Study+projid FID: %d", n_mapped))
message(sprintf("  Could NOT be mapped: %d", n_unmapped))

if (n_unmapped > 0) {
  unmapped <- result %>% filter(is.na(fid))
  message("  Unmapped samples:")
  for (i in seq_len(nrow(unmapped))) {
    message(sprintf("    %s (individualID: %s)",
                    unmapped$topmed_specimenID[i], unmapped$individualID[i]))
  }
}

# Show breakdown by platform
n_illumina <- sum(grepl("^(MAP|ROS)", result$topmed_specimenID))
n_affy     <- sum(grepl("^11AD",      result$topmed_specimenID))
message(sprintf("\n  Illumina (MAP/ROS) IDs: %d", n_illumina))
message(sprintf("  Affymetrix (11AD) IDs:  %d", n_affy))
message("\nSample FID examples:")
message(sprintf("  Illumina: %s -> %s",
                result$topmed_specimenID[grepl("^(MAP|ROS)", result$topmed_specimenID)][1],
                result$fid[grepl("^(MAP|ROS)", result$topmed_specimenID)][1]))
message(sprintf("  Affymetrix: %s -> %s",
                result$topmed_specimenID[grepl("^11AD", result$topmed_specimenID)][1],
                result$fid[grepl("^11AD", result$topmed_specimenID)][1]))

# ---- Write samples_to_keep.txt ----
# Format: FID  IID  (both = Study+projid)
# plink2 --keep expects FID IID; pipeline also uses this to match phenotypes
mapped <- result %>%
  filter(!is.na(fid)) %>%
  transmute(
    `#FID` = fid,
    IID    = fid
  ) %>%
  distinct() %>%
  arrange(`#FID`)

out_file <- file.path(OUTPUT_DIR, "samples_to_keep.txt")
fwrite(mapped, out_file, sep = "\t", quote = FALSE)
message(sprintf("\nWrote %d samples to: %s", nrow(mapped), out_file))

# Also write a mapping reference file (topmed_specID -> FID) for tracing
mapping_ref <- result %>%
  filter(!is.na(fid)) %>%
  select(topmed_specimenID, individualID, fid) %>%
  arrange(topmed_specimenID)
fwrite(mapping_ref,
       file.path(OUTPUT_DIR, "sample_id_mapping.txt"),
       sep = "\t", quote = FALSE)
message(sprintf("Wrote ID mapping reference to: %s/sample_id_mapping.txt", OUTPUT_DIR))
