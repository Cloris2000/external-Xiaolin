#!/usr/bin/env python3
"""
Assemble the per-sample PC table behind Figure 1 panels A, B and D.

Inputs are the plink2 --score projections written by figure1_genotype_pca.sh
(one per cohort plus the 1000G reference), the reference .psam (superpopulation
labels), the eigenvalues (variance explained) and the repo's self-reported
ancestry assignments.

Genetic ancestry is assigned by k-nearest-neighbours on the first 6 reference
PCs: each cohort sample takes the majority superpopulation among its k nearest
1000G neighbours.  A sample whose majority vote is below --knn-min-frac is left
UNCERTAIN rather than forced into a group.

Note the distinction this table preserves:
  - `genetic_ancestry`  inferred here from genotypes vs the 1000G reference
  - `reported_ancestry` self-reported race/ethnicity from cohort metadata
    (docs/ancestry_specific/ancestry_sample_assignments.tsv)
These are different quantities and are deliberately kept in separate columns.
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

SUPERPOPS = ["EUR", "AFR", "AMR", "EAS", "SAS"]

# Cohort display order for the figure (matches the meta-analysis config order).
COHORT_ORDER = [
    "ROSMAP", "ROSMAP_array", "Mayo", "MSBB",
    "CMC_MSSM", "CMC_PENN", "CMC_PITT",
    "GTEx_v10", "NABEC",
    "NIMH_HBCC_1M", "NIMH_HBCC_h650", "NIMH_HBCC_Omni5M",
    "GVEX", "AMP_AD_Rush", "AMP_AD_Mayo",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proj-dir", required=True,
                   help="dir holding proj_<name>.sscore files")
    p.add_argument("--kg-psam", required=True)
    p.add_argument("--eigenval", required=True)
    p.add_argument("--ancestry-assignments", required=True)
    p.add_argument("--out-samples", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--n-pcs", type=int, default=10)
    p.add_argument("--knn-k", type=int, default=20)
    p.add_argument("--knn-min-frac", type=float, default=0.6,
                   help="min neighbour fraction to call an ancestry")
    return p.parse_args()


def read_sscore(path, n_pcs):
    """Read a plink2 .sscore into (ids, PC matrix)."""
    df = pd.read_csv(path, sep="\t")
    id_col = "IID" if "IID" in df.columns else df.columns[0]
    pc_cols = [c for c in df.columns if c.endswith("_AVG")][:n_pcs]
    if not pc_cols:
        sys.exit(f"ERROR: no *_AVG score columns in {path}")
    return df[id_col].astype(str).values, df[pc_cols].values.astype(float), pc_cols


def knn_assign(ref_pcs, ref_labels, qry_pcs, k, min_frac, n_dims=6):
    """Majority superpopulation among the k nearest reference samples."""
    R = ref_pcs[:, :n_dims]
    Q = qry_pcs[:, :n_dims]
    labels, fracs = [], []
    # Chunked to keep the distance matrix small.
    for start in range(0, len(Q), 512):
        chunk = Q[start:start + 512]
        d = ((chunk[:, None, :] - R[None, :, :]) ** 2).sum(axis=2)
        idx = np.argpartition(d, kth=k, axis=1)[:, :k]
        for row in idx:
            votes = Counter(ref_labels[row])
            lab, n = votes.most_common(1)[0]
            frac = n / k
            labels.append(lab if frac >= min_frac else "UNCERTAIN")
            fracs.append(frac)
    return np.array(labels), np.array(fracs)


def main():
    args = parse_args()
    proj_dir = Path(args.proj_dir)

    # --- reference ---------------------------------------------------------
    ref_file = proj_dir / "proj_1000G.sscore"
    if not ref_file.exists():
        sys.exit(f"ERROR: missing reference projection {ref_file}")
    ref_ids, ref_pcs, pc_cols = read_sscore(ref_file, args.n_pcs)

    psam = pd.read_csv(args.kg_psam, sep="\t", dtype=str)
    psam.columns = [c.lstrip("#") for c in psam.columns]
    id_col = "IID" if "IID" in psam.columns else psam.columns[0]
    pop_map = dict(zip(psam[id_col].astype(str), psam["SuperPop"]))
    ref_labels = np.array([pop_map.get(i, "NA") for i in ref_ids])

    keep = ref_labels != "NA"
    ref_ids, ref_pcs, ref_labels = ref_ids[keep], ref_pcs[keep], ref_labels[keep]
    print(f"[assemble] reference samples: {len(ref_ids):,} "
          f"({dict(Counter(ref_labels))})", flush=True)

    # --- cohorts -----------------------------------------------------------
    frames = []
    ref_df = pd.DataFrame(ref_pcs, columns=[f"PC{i+1}" for i in range(ref_pcs.shape[1])])
    ref_df.insert(0, "IID", ref_ids)
    ref_df.insert(0, "cohort", "1000G")
    ref_df["dataset"] = "reference"
    ref_df["genetic_ancestry"] = ref_labels
    ref_df["knn_frac"] = np.nan
    frames.append(ref_df)

    for f in sorted(proj_dir.glob("proj_*.sscore")):
        name = f.stem.replace("proj_", "")
        if name == "1000G":
            continue
        ids, pcs, _ = read_sscore(f, args.n_pcs)
        lab, frac = knn_assign(ref_pcs, ref_labels, pcs,
                               args.knn_k, args.knn_min_frac)
        d = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(pcs.shape[1])])
        d.insert(0, "IID", ids)
        d.insert(0, "cohort", name)
        d["dataset"] = "cohort"
        d["genetic_ancestry"] = lab
        d["knn_frac"] = frac
        frames.append(d)
        print(f"[assemble] {name:18s} n={len(ids):5d}  "
              f"{dict(Counter(lab))}", flush=True)

    samples = pd.concat(frames, ignore_index=True)

    # --- self-reported labels, kept separate -------------------------------
    rep = pd.read_csv(args.ancestry_assignments, sep="\t", dtype=str)
    rep = rep.rename(columns={"ancestry": "reported_ancestry"})
    samples = samples.merge(
        rep[["cohort", "IID", "reported_ancestry"]],
        on=["cohort", "IID"], how="left",
    )

    # --- variance explained ------------------------------------------------
    ev = pd.read_csv(args.eigenval, header=None)[0].values.astype(float)
    pct = 100 * ev / ev.sum()
    for i in range(min(args.n_pcs, len(pct))):
        samples[f"PC{i+1}_pct_var"] = round(float(pct[i]), 2)
    print("[assemble] variance explained: " +
          ", ".join(f"PC{i+1}={pct[i]:.1f}%" for i in range(min(5, len(pct)))),
          flush=True)

    cohort_cat = pd.Categorical(
        samples["cohort"],
        categories=["1000G"] + COHORT_ORDER, ordered=True)
    samples = samples.assign(_o=cohort_cat).sort_values(["_o", "IID"]).drop(columns="_o")
    samples.to_csv(args.out_samples, sep="\t", index=False)

    # --- per-cohort ancestry composition (panel D) -------------------------
    coh = samples[samples.dataset == "cohort"]
    summ = (coh.groupby(["cohort", "genetic_ancestry"])
               .size().rename("n").reset_index())
    tot = coh.groupby("cohort").size().rename("n_total").reset_index()
    summ = summ.merge(tot, on="cohort")
    summ["pct"] = (100 * summ.n / summ.n_total).round(2)

    # Concordance with self-report, where a self-reported label exists.
    both = coh.dropna(subset=["reported_ancestry"])
    both = both[both.reported_ancestry.isin(SUPERPOPS)]
    conc = (both.assign(agree=both.genetic_ancestry == both.reported_ancestry)
                .groupby("cohort").agree.agg(["sum", "count"]))
    conc["pct_concordant"] = (100 * conc["sum"] / conc["count"]).round(1)
    summ = summ.merge(
        conc[["pct_concordant"]].reset_index(), on="cohort", how="left")

    summ.to_csv(args.out_summary, sep="\t", index=False)

    print(f"[assemble] wrote {args.out_samples} ({len(samples):,} rows)")
    print(f"[assemble] wrote {args.out_summary}")
    if len(both):
        overall = 100 * (both.genetic_ancestry == both.reported_ancestry).mean()
        print(f"[assemble] genetic vs self-reported concordance: {overall:.1f}% "
              f"(n={len(both):,} with both labels)")


if __name__ == "__main__":
    main()
