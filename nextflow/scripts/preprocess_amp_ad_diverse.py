#!/usr/bin/env python3
"""
Preprocessing script for AMP-AD Diverse Cohort pipeline integration.

Generates per-sub-cohort metadata and cleaned count matrices for use with
combined_pipeline_v2.nf. Must be run once before launching any
nextflow.config.combined.amp_ad_* pipeline config.

Outputs (written to data_input/amp_ad_diverse/):
  Columbia_DLPFC_metadata.csv           -- metadata keyed by individualID
  MSSM_DLPFC_metadata.csv               -- metadata keyed by specimenID (RNA)
  Rush_DLPFC_metadata.csv               -- metadata keyed by specimenID (RNA = Div_*)
  Mayo_DLPFC_metadata.csv               -- metadata keyed by specimenID (RNA = {id}_DLPFC)
  MSSM_Gene_Count_Matrix_clean.csv      -- count matrix with .final suffix stripped
  Rush_gene_counts_matrix_clean.csv     -- count matrix with _S{N} suffix stripped
  Mayo_Emory_Gene_Count_Matrix.csv      -- Mayo/Emory shared count matrix (columns = {id}_DLPFC)

ID mapping chain per sub-cohort:
  Columbia : count matrix col = individualID (NYBB_100)
             -> biospecimen WGS: individualID -> specimenID (NYBB_100WGS) = VCF ID
  MSSM     : count matrix col (after strip .final) = RNA specimenID (122253)
             -> biospecimen rnaSeq: specimenID -> individualID (33690)
             -> biospecimen WGS:    individualID -> specimenID (122253-D)  = VCF ID
  Rush     : count matrix col (after strip _S*) = RNA specimenID (Div_100)
             -> biospecimen rnaSeq: specimenID -> individualID (R5738585)
             -> WGS specimenID = individualID (R5738585) = VCF ID  (direct)
  Mayo     : count matrix col = RNA specimenID ({id}_DLPFC, e.g. 1005_DLPFC)
             -> biospecimen WGS: individualID ({id}) -> specimenID ({id}_DLPFC_WGS) = VCF ID

Usage:
  python3 scripts/preprocess_amp_ad_diverse.py
  python3 scripts/preprocess_amp_ad_diverse.py --out_dir /path/to/output
"""

import argparse
import os
import re
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
META_DIR = "/external/rprshnas01/netdata_kcni/stlab/AMP_AD_Diverse/Metadata"
BULK_DIR = "/external/rprshnas01/netdata_kcni/stlab/AMP_AD_Diverse/Bulk_RNA_Seq_and_QC"
DEFAULT_OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data_input", "amp_ad_diverse"
)

DLPFC_LABEL = "dorsolateral prefrontal cortex"


def load_shared_metadata(meta_dir: str):
    individual = pd.read_csv(os.path.join(meta_dir, "AMP-AD_DiverseCohorts_individual_metadata.csv"))
    biospecimen = pd.read_csv(os.path.join(meta_dir, "AMP-AD_DiverseCohorts_biospecimen_metadata.csv"))
    rna_assay = pd.read_csv(os.path.join(meta_dir, "AMP-AD_DiverseCohorts_assay_RNAseq_metadata.csv"))
    print(f"  individual metadata:  {len(individual)} rows")
    print(f"  biospecimen metadata: {len(biospecimen)} rows")
    print(f"  RNA assay metadata:   {len(rna_assay)} rows")
    return individual, biospecimen, rna_assay


def build_combined_metadata(
    individual: pd.DataFrame,
    biospecimen: pd.DataFrame,
    rna_assay: pd.DataFrame,
    cohort_name: str,
    tissue: str = DLPFC_LABEL,
) -> pd.DataFrame:
    """
    Build a combined metadata DataFrame for one cohort+tissue.

    Steps:
      1. Filter biospecimen to rnaSeq + tissue to get (specimenID, individualID) map.
      2. Join RNA assay metadata onto that map (specimenID is the key).
      3. Join individual metadata (clinical covariates) by individualID.
      4. Filter to the requested cohort.
    """
    rna_bio = biospecimen[
        (biospecimen["assay"] == "rnaSeq") &
        (biospecimen["tissue"].str.lower() == tissue.lower())
    ][["individualID", "specimenID"]].drop_duplicates()

    # Attach RNA assay batch/QC info
    rna_merged = rna_assay.merge(rna_bio, on="specimenID", how="inner")

    # Attach individual (clinical) metadata
    combined = rna_merged.merge(individual, on="individualID", how="inner")

    # Filter to the requested cohort
    if "dataContributionGroup" in combined.columns:
        combined = combined[combined["dataContributionGroup"] == cohort_name]

    print(f"  {cohort_name} {tissue[:5]} combined metadata: {len(combined)} rows")
    return combined.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Columbia
# ---------------------------------------------------------------------------
def process_columbia(individual, biospecimen, rna_assay, out_dir):
    print("\n=== Columbia ===")
    meta = build_combined_metadata(individual, biospecimen, rna_assay, "Columbia")
    # Count matrix columns = individualID (NYBB_100); metadata is keyed by individualID.
    # The pipeline uses col_sample_id_for_matching = "individualID".
    out_path = os.path.join(out_dir, "Columbia_DLPFC_metadata.csv")
    meta.to_csv(out_path, index=False)
    print(f"  Written: {out_path}")


# ---------------------------------------------------------------------------
# MSSM
# ---------------------------------------------------------------------------
def process_mssm(individual, biospecimen, rna_assay, out_dir):
    print("\n=== MSSM ===")

    # Build combined metadata keyed by specimenID (RNA)
    meta = build_combined_metadata(individual, biospecimen, rna_assay, "MSSM")
    out_meta = os.path.join(out_dir, "MSSM_DLPFC_metadata.csv")
    meta.to_csv(out_meta, index=False)
    print(f"  Written: {out_meta}")

    # Clean count matrix: strip '.final' suffix from all column names except the first (GeneID)
    cm_path = os.path.join(BULK_DIR, "MSSM", "MSSM_Gene_Count_Matrix.csv")
    cm = pd.read_csv(cm_path, index_col=0)
    original_cols = list(cm.columns)
    new_cols = [re.sub(r"\.final$", "", c) for c in original_cols]
    n_renamed = sum(1 for a, b in zip(original_cols, new_cols) if a != b)
    cm.columns = new_cols
    out_cm = os.path.join(out_dir, "MSSM_Gene_Count_Matrix_clean.csv")
    cm.to_csv(out_cm)
    print(f"  MSSM count matrix: renamed {n_renamed}/{len(original_cols)} columns (stripped .final)")
    print(f"  Written: {out_cm}")


# ---------------------------------------------------------------------------
# Rush
# ---------------------------------------------------------------------------
def process_rush(individual, biospecimen, rna_assay, out_dir):
    print("\n=== Rush ===")

    # Build combined metadata keyed by specimenID (RNA = Div_*)
    meta = build_combined_metadata(individual, biospecimen, rna_assay, "Rush")
    out_meta = os.path.join(out_dir, "Rush_DLPFC_metadata.csv")
    meta.to_csv(out_meta, index=False)
    print(f"  Written: {out_meta}")

    # Clean count matrix: strip '_S{N}' suffix from sample column names
    cm_path = os.path.join(BULK_DIR, "rush_cohort_output", "gene_counts_matrix.tsv")
    cm = pd.read_csv(cm_path, sep="\t", index_col=0)
    original_cols = list(cm.columns)
    new_cols = [re.sub(r"_S\d+$", "", c) for c in original_cols]
    n_renamed = sum(1 for a, b in zip(original_cols, new_cols) if a != b)
    # Verify uniqueness after stripping (no collisions expected)
    if len(new_cols) != len(set(new_cols)):
        dupes = [c for c in new_cols if new_cols.count(c) > 1]
        print(f"  WARNING: {len(set(dupes))} duplicate column names after stripping _S suffix: {set(dupes)}")
    cm.columns = new_cols
    out_cm = os.path.join(out_dir, "Rush_gene_counts_matrix_clean.csv")
    cm.to_csv(out_cm)
    print(f"  Rush count matrix: renamed {n_renamed}/{len(original_cols)} columns (stripped _S{{N}})")
    print(f"  Written: {out_cm}")


# ---------------------------------------------------------------------------
# Mayo
# ---------------------------------------------------------------------------
def process_mayo(individual, biospecimen, rna_assay, out_dir):
    """
    Mayo DLPFC preprocessing.

    Count matrix columns: {id}_DLPFC  (e.g. 1005_DLPFC, 114_DLPFC)
    These match the RNA specimenID directly.
    WGS specimenID in VCF: {id}_DLPFC_WGS (e.g. 1005_DLPFC_WGS)

    The pipeline uses:
      col_sample_id_for_matching = "specimenID"  -> matches count matrix col = RNA specimenID
      fid_method = "biospec_specimen" + biospec_assay_filter = "wholeGenomeSeq"
        -> FID = WGS specimenID ({id}_DLPFC_WGS) which is the VCF sample ID.

    The count matrix also contains Emory samples (29476_DLPFC, etc.) which are NOT
    in the Mayo metadata and will be naturally filtered out by the pipeline's
    tissue+metadata join.
    """
    print("\n=== Mayo ===")

    meta = build_combined_metadata(individual, biospecimen, rna_assay, "Mayo")
    out_meta = os.path.join(out_dir, "Mayo_DLPFC_metadata.csv")
    meta.to_csv(out_meta, index=False)
    print(f"  Written: {out_meta}")

    # Count matrix shared with Emory — no column renaming needed
    # (columns are already in {id}_DLPFC format matching specimenID)
    cm_path = os.path.join(BULK_DIR, "Mayo_emory", "Gene_Count_Matrix.csv")
    cm = pd.read_csv(cm_path, index_col=0)
    out_cm = os.path.join(out_dir, "Mayo_Emory_Gene_Count_Matrix.csv")
    cm.to_csv(out_cm)
    print(f"  Mayo/Emory count matrix: {cm.shape[1]} samples, {cm.shape[0]} genes")
    print(f"  Written: {out_cm}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--out_dir", default=DEFAULT_OUT,
        help=f"Output directory (default: {DEFAULT_OUT})"
    )
    parser.add_argument(
        "--meta_dir", default=META_DIR,
        help=f"AMP-AD Diverse metadata directory (default: {META_DIR})"
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"Output directory: {args.out_dir}")

    print("\nLoading shared metadata...")
    individual, biospecimen, rna_assay = load_shared_metadata(args.meta_dir)

    process_columbia(individual, biospecimen, rna_assay, args.out_dir)
    process_mssm(individual, biospecimen, rna_assay, args.out_dir)
    process_rush(individual, biospecimen, rna_assay, args.out_dir)
    process_mayo(individual, biospecimen, rna_assay, args.out_dir)

    print("\nDone. Pre-processed files:")
    for f in sorted(os.listdir(args.out_dir)):
        fpath = os.path.join(args.out_dir, f)
        size_kb = os.path.getsize(fpath) // 1024
        print(f"  {f}  ({size_kb} KB)")


if __name__ == "__main__":
    main()
