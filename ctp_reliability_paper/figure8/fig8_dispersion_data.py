"""Per-supertype layer dispersion vs concordance, and the depth negative control
(PI suggestion 2026-10-23).

Hypothesis: a supertype spread across many layers is less sensitive to where
the block was cut, so it should (a) be more concordant across datasets and
(b) show a flatter relationship between its own cross-dataset disagreement
and the depth mismatch of the two blocks.

Per supertype:
  from cortdepth's Xenium reference (prior_DFC_2026.csv):
    sigma         depth spread of the supertype's cells (IQR / 1.349, >= 0.02)
    p_wm          share of its cells in white matter (grey/white)
    depth_median, n_cells_ref
  from common/data/reliability_summary.tsv (Figure 1 metric):
    median_icc_sym, q25, q75, n_pairs   median symmetric ICC over the 5 pairs
  depth sensitivity, recomputed here on the 511 passing donors of
  data/fig8_depth_pair_donor.tsv:
    for each pair and donor, |residual CLR_A - residual CLR_B| of this
    supertype (CLR within its class set -- inhibitory, glial, or excitatory --
    residual = minus the dataset's across-donor mean, as fig8_data.py), then
    pooled over pairs: slope of that on wasserstein with pair as a fixed
    effect (depth_slope, depth_slope_se, depth_slope_p), its partial R2, and
    the Spearman rho of the pair-demeaned quantities; n_donors.
  mean_share: mean across the five pairs' datasets of the supertype's mean
    share of its class set (abundance, a known correlate of ICC)
Excitatory supertypes are computed too but flagged (class Excitatory): they
define x, so their depth_slope is circular and is drawn hollow, for reference.

Output: data/fig8_supertype_dispersion.tsv, one row per supertype.
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig8_dispersion_data.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as t_dist

CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
HERE = Path(__file__).resolve().parent
COMMON = HERE.parent / "common"
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/"
           "refs/taxonomy_DFC_2026.tsv")
PRIOR = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/cortical-depth/src/cortdepth/"
             "data/prior_DFC_2026.csv")
PAIRS = {
    "Green x Mathys": ("Green_DFC_person", "Mathys_DFC_person"),
    "Green x Multiome": ("Green_DFC_person", "Multiome2025_DLPFC_DFC_person"),
    "Mathys x Multiome": ("Mathys_DFC_person", "Multiome2025_DLPFC_DFC_person"),
    "Multiome x PsychAD": ("Multiome2025_DLPFC_DFC_person", "PsychAD_full_DFC_person"),
    "BrainSCOPE CMC x PsychAD": ("BrainSCOPE_CMC_DFC_person", "PsychAD_full_DFC_person"),
}
PSEUDO = 0.5
MIN_DONORS = 100


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def load_counts(name):
    c = pd.read_parquet(CACHE / f"{name}_cells.parquet", columns=["donor", "supertype"])
    c["donor"] = c.donor.astype(str)
    return pd.crosstab(c.donor, c.supertype)


def per_type_absdiff(tA, tB, donors, types):
    """Long table (donor, supertype, absdiff) for one class set in one pair."""
    types = [t for t in types if t in tA.columns and t in tB.columns]
    dA = tA.reindex(index=donors, columns=types).fillna(0).values
    dB = tB.reindex(index=donors, columns=types).fillna(0).values
    keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
    yA, yB = clr(dA[keep]), clr(dB[keep])
    diff = np.abs((yA - yA.mean(0)) - (yB - yB.mean(0)))
    share = 0.5 * ((dA[keep] / dA[keep].sum(1, keepdims=True)).mean(0)
                   + (dB[keep] / dB[keep].sum(1, keepdims=True)).mean(0))
    long = pd.DataFrame(diff, index=np.array(donors)[keep], columns=types).stack()
    long = long.rename("absdiff").reset_index().rename(columns={"level_0": "donor", "level_1": "supertype"})
    return long, pd.Series(share, index=types)


def fixed_effects_slope(df):
    """absdiff ~ wasserstein + pair (dummies). Returns slope, se, p, partial R2, demeaned rho."""
    pairs = pd.get_dummies(df.pair, drop_first=True).values.astype(float)
    X = np.column_stack([np.ones(len(df)), df.wasserstein.values, pairs])
    y = df.absdiff.values
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    dof = len(y) - X.shape[1]
    s2 = (r @ r) / dof
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(cov[1, 1])
    p = 2 * t_dist.sf(abs(b[1] / se), dof)
    X0 = np.delete(X, 1, axis=1)
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
    r0 = y - X0 @ b0
    pr2 = 1 - (r @ r) / (r0 @ r0)
    g = df.groupby("pair")
    rho = spearmanr(df.wasserstein - g.wasserstein.transform("mean"),
                    df.absdiff - g.absdiff.transform("mean")).correlation
    return b[1], se, p, pr2, rho


def main():
    tax = pd.read_csv(TAX, sep="\t")
    cls_of = dict(zip(tax.supertype_label, tax.class_id))
    sets = {k: set(g.supertype_label) for k, g in tax.groupby("class_id")}
    don = pd.read_csv(HERE / "data/fig8_depth_pair_donor.tsv", sep="\t")
    don = don[don["pass"]]

    cache, long, shares = {}, [], []
    for pair, (a, b) in PAIRS.items():
        for n in (a, b):
            if n not in cache:
                cache[n] = load_counts(n)
        g = don[don.pair == pair]
        for cls, types in sets.items():
            L, sh = per_type_absdiff(cache[a], cache[b], list(g.donor), sorted(types))
            L["pair"] = pair
            long.append(L)
            shares.append(sh.rename(pair))
        print(f"  {pair:26s} {len(g)} donors")
        later = {x for p in list(PAIRS.values())[list(PAIRS).index(pair) + 1:] for x in p}
        for n in list(cache):
            if n not in later:
                del cache[n]
    long = pd.concat(long, ignore_index=True)
    long = long.merge(don[["pair", "donor", "wasserstein"]], on=["pair", "donor"])
    mean_share = pd.concat(shares, axis=1).mean(axis=1).rename("mean_share")

    rows = []
    for st, df in long.groupby("supertype"):
        if len(df) < MIN_DONORS:
            continue
        slope, se, p, pr2, rho = fixed_effects_slope(df)
        rows.append(dict(supertype=st, n_donors=len(df), depth_slope=slope, depth_slope_se=se,
                         depth_slope_p=p, depth_partial_r2=pr2, depth_rho=rho))
    sens = pd.DataFrame(rows)

    prior = pd.read_csv(PRIOR).rename(columns={"n_cells": "n_cells_ref"})
    rel = pd.read_csv(COMMON / "data/reliability_summary.tsv", sep="\t")
    rel = rel[rel.level == "supertype"][["cell_type", "median_icc_sym", "q25", "q75", "n_pairs"]]
    rel = rel.rename(columns={"cell_type": "supertype"})

    out = (tax[["supertype_label", "subclass_label", "class_id", "reactive_state"]]
           .rename(columns={"supertype_label": "supertype", "subclass_label": "subclass", "class_id": "cls"})
           .merge(prior[["supertype", "n_cells_ref", "depth_median", "sigma", "p_wm"]], on="supertype", how="left")
           .merge(rel, on="supertype", how="left")
           .merge(sens, on="supertype", how="left")
           .merge(mean_share, left_on="supertype", right_index=True, how="left"))
    out.to_csv(HERE / "data/fig8_supertype_dispersion.tsv", sep="\t", index=False)
    ok = out.dropna(subset=["sigma", "median_icc_sym"])
    print(f"\n{len(out)} supertypes; {len(ok)} with reference depth and ICC; "
          f"{out.depth_slope.notna().sum()} with a depth slope (>= {MIN_DONORS} donors)")
    for cls, g in ok.groupby("cls"):
        r1 = spearmanr(g.sigma, g.median_icc_sym)
        gg = g.dropna(subset=["depth_slope"])
        r2 = spearmanr(gg.sigma, gg.depth_slope) if len(gg) > 3 else None
        print(f"  {cls:11s} n={len(g):3d}  sigma vs ICC rho={r1.correlation:+.2f} p={r1.pvalue:.2g}"
              + (f"   sigma vs depth slope rho={r2.correlation:+.2f} p={r2.pvalue:.2g} (n={len(gg)})" if r2 else ""))
    print(f"wrote data/fig8_supertype_dispersion.tsv")


if __name__ == "__main__":
    main()
