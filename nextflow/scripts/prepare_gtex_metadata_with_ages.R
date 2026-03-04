#!/usr/bin/env Rscript
# Prepare GTEx metadata with merged specific ages
# This script replicates the age merging from your original script

library(readr)
library(dplyr)

# Load main metadata
FC_sample_metadata_cleaned <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/data/GTEx_FC_sample_metadata_cleaned.csv")

# Load GTEx age file with specific ages
GTEx_meta_age <- read_tsv("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/GTEx_WGS/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt",
                         comment = "#",
                         skip = 1,
                         col_names = TRUE,
                         show_col_types = FALSE)

# Extract subject ID from SAMPID (first two parts: GTEX-XXXXX)
FC_sample_metadata_cleaned$subject_id <- sapply(strsplit(FC_sample_metadata_cleaned$SAMPID, "-"), 
                                                 function(x) paste(x[1], x[2], sep = "-"))

# Merge with age data
FC_sample_metadata_merged <- left_join(FC_sample_metadata_cleaned, GTEx_meta_age, 
                                       by = c("subject_id" = "SUBJID"))

# Save merged metadata
write.csv(FC_sample_metadata_merged, 
          "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/data/GTEx_FC_sample_metadata_with_ages.csv",
          row.names = FALSE)

cat("✅ Created merged metadata file with", nrow(FC_sample_metadata_merged), "samples\n")
cat("   Columns include: SEX.y, AGE.y (specific numeric ages)\n")
cat("   File: GTEx_FC_sample_metadata_with_ages.csv\n")
