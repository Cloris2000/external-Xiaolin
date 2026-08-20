#!/usr/bin/env python3
"""
Build ROSMAP_Mathys single-nucleus cell-type proportions (19-class Hodge) from
the per-cell label-transfer table:

  data/Mathys/mathys_hodge_subclass_fine_hodge19_labels.csv

Rules
-----
- Use only cells with passed_QC == True.
- Use the final 19-class column `Hodge_subclass_fine_hodge19` as the annotation.
- Map its label spellings to the pipeline's canonical Hodge subclass names
  (matching data_input/sn_rosmap_green_hodge/cell_proportions.csv).
- Proportion = (# cells of a type in a donor) / (# QC cells in that donor).
- Restrict to donors present in the ROSMAP genotype (matched donor list),
  so the output lines up with the WGS-derived genotypes.

Outputs
-------
- cell_proportions.csv : individualID + one column per retained cell type
- occupancy_report.tsv : per-cell-type #cells and #donors with nonzero cells
"""
import argparse
import pandas as pd

# Mathys Hodge19 label spelling -> pipeline canonical name
LABEL_MAP = {
    "Astrocyte":      "Astrocyte",
    "Endothelial":    "Endothelial",
    "IT":             "IT",
    "L4 IT":          "L4.IT",
    "L5/6 IT Car3":   "L5.6.IT.Car3",
    "L5/6 NP":        "L5.6.NP",
    "L5 ET":          "L5.ET",
    "L6 CT":          "L6.CT",
    "L6b":            "L6b",
    "LAMP5":          "LAMP5",
    "Microglia":      "Microglia",
    "OPC":            "OPC",
    "Oligodendrocyte":"Oligodendrocyte",
    "PAX6":           "PAX6",
    "PVALB":          "PVALB",
    "Pericyte":       "Pericyte",
    "SST":            "SST",
    "VIP":            "VIP",
    "VLMC":           "VLMC",
}

# Canonical 19-class column order (same as sn_rosmap_green_hodge)
CANON_ORDER = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
    "L5.ET", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", "Oligodendrocyte",
    "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--matched-donors", required=True,
                    help="one R* individualID per line; restrict output to these")
    ap.add_argument("--out-proportions", required=True)
    ap.add_argument("--out-occupancy", required=True)
    ap.add_argument("--min-cells", type=int, default=200,
                    help="drop donors with fewer than this many QC cells")
    ap.add_argument("--drop-celltypes", default="",
                    help="comma-separated canonical names to drop (too rare)")
    args = ap.parse_args()

    matched = set()
    with open(args.matched_donors) as fh:
        for line in fh:
            s = line.strip()
            if s:
                matched.add(s)
    print(f"matched donor list: {len(matched)} donors")

    drop = {c for c in args.drop_celltypes.split(",") if c}

    usecols = ["Donor ID", "passed_QC", "Hodge_subclass_fine_hodge19"]
    print("reading labels (this is a ~200MB csv) ...")
    df = pd.read_csv(args.labels, usecols=usecols, dtype=str)
    n0 = len(df)

    # QC filter (passed_QC stored as 'true'/'false' strings)
    df = df[df["passed_QC"].str.lower() == "true"]
    print(f"cells: {n0} -> {len(df)} after passed_QC==true")

    # restrict to genotype-matched donors
    df = df[df["Donor ID"].isin(matched)]
    print(f"cells after restricting to matched donors: {len(df)}")

    # map labels
    df["ct"] = df["Hodge_subclass_fine_hodge19"].map(LABEL_MAP)
    unmapped = df["ct"].isna().sum()
    if unmapped:
        bad = df.loc[df["ct"].isna(), "Hodge_subclass_fine_hodge19"].unique()
        raise SystemExit(f"ERROR: {unmapped} cells with unmapped labels: {bad}")

    # counts donor x celltype
    counts = df.groupby(["Donor ID", "ct"]).size().unstack(fill_value=0)

    # occupancy report (before donor min-cells filter)
    occ = pd.DataFrame({
        "n_cells": counts.sum(axis=0),
        "n_donors_nonzero": (counts > 0).sum(axis=0),
    }).reindex([c for c in CANON_ORDER if c in counts.columns])
    occ.to_csv(args.out_occupancy, sep="\t")
    print("\n=== occupancy (per cell type) ===")
    print(occ.to_string())

    # donor min-cells filter
    donor_tot = counts.sum(axis=1)
    keep_donors = donor_tot[donor_tot >= args.min_cells].index
    dropped = counts.shape[0] - len(keep_donors)
    print(f"\ndonors: {counts.shape[0]} -> {len(keep_donors)} "
          f"(dropped {dropped} with < {args.min_cells} QC cells)")
    counts = counts.loc[keep_donors]
    donor_tot = donor_tot.loc[keep_donors]

    # proportions
    prop = counts.div(donor_tot, axis=0)

    # column order + drops
    cols = [c for c in CANON_ORDER if c in prop.columns and c not in drop]
    prop = prop[cols]
    prop.index.name = "individualID"
    prop = prop.reset_index()

    prop.to_csv(args.out_proportions, index=False)
    print(f"\nwrote {args.out_proportions}: {prop.shape[0]} donors x "
          f"{len(cols)} cell types")
    print("cell types:", ", ".join(cols))


if __name__ == "__main__":
    main()
