"""Figure 8 (simple layer-marker version): sampled-depth mismatch from two
subclass layer markers, no depth reference (user request 2026-10-23; made the
Figure 8 design 2026-11-03, replacing the cortdepth Wasserstein axis).

Among a sample's excitatory nuclei, L2/3 IT marks the upper layers and L6 IT the
deep layers:

    upper_frac = n(L2/3 IT) / (n(L2/3 IT) + n(L6 IT))      1 = all upper, 0 = all deep

The depth mismatch between the two samples of a donor is |upper_frac_A - upper_frac_B|.

Pairs and donors (2026-11-03): the eight pairs of common/reliability_stats.py,
from data/fig8_depth_pair_donor.tsv (fig8_data.py): registry person_ids, the
PI's `pass` QC flags (frontal, not a duplicate, no excitatory label collapse,
sex-concordant, identity not assumed) and disc_inh. Marker counts come from the
registry composition subclass tables (reliability_stats.COUNTS), the same
counts behind every pair statistic.

Two donor-level y measures over INHIBITORY supertypes (excitatory types define x):
  disc_inh   mean |residual CLR_A - residual CLR_B| (fig8_data.py; a distance,
             higher = more discordant); residual = CLR within the inhibitory set
             minus the type's mean over the pair's shared donors in that dataset
  agree_inh  Pearson r between the donor's two inhibitory profiles, as Figure 3c
             (fig3_data.py): CLR of inhibitory types with mean share >= 0.5% in
             both datasets (among donors with >= 50 inhibitory nuclei in both),
             centred and scaled per type within each dataset, then correlated
             across types per donor (higher = more concordant)
  tv_inh     % of the donor's inhibitory neurons that differ between the two
             datasets (2026-11-03, the Figure 8 y): half the sum over inhibitory
             supertypes of |proportion_A - proportion_B|, proportions within
             inhibitory nuclei (total variation distance; 0 = identical, 100 =
             no shared types). Plain proportions, so common types dominate; the
             part shared by all donors of a pair (lab/protocol) is a per-pair
             baseline, absorbed by the pair term of the pooled model
  icc_inh    the same profiles as agree_inh, scored per donor with the paper's
             symmetric ICC across types, 2 cov(a, b) / (var a + var b)
             (2026-11-03): equals agree_inh when the donor's two profiles have
             equal spread and is lower when one dataset shows a stronger or
             weaker version of the same pattern (higher = more concordant)

Output: data/fig8_simple_depth.tsv = fig8_depth_pair_donor.tsv plus, per side,
n_l23, n_l6it, upper_frac; abs_d_upper; tv_inh; n_inh_types, agree_inh and icc_inh.
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig8_simple_check_data.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import reliability_stats as rs  # noqa: E402

UPPER, DEEP = "L2_3_IT", "L6_IT"          # subclass_id columns of the count tables
MIN_CLASS_CELLS, MIN_SHARE = 50, 0.005     # as fig3_data.py


def marker_counts(name):
    t = pd.read_parquet(rs.COUNTS / f"{name}__subclass.parquet")
    t.index = t.index.astype(str)
    t = t[[UPPER, DEEP]].rename(columns={UPPER: "n_l23", DEEP: "n_l6it"})
    t["upper_frac"] = t.n_l23 / (t.n_l23 + t.n_l6it)
    return t


def tv_inh(ta, tb, donors, inh_types):
    """% of inhibitory neurons that differ: half the summed |share_A - share_B|."""
    inh = [c for c in inh_types if c in ta.columns or c in tb.columns]
    A = ta.reindex(index=donors, columns=inh, fill_value=0).values.astype(float)
    B = tb.reindex(index=donors, columns=inh, fill_value=0).values.astype(float)
    pa, pb = A / A.sum(1, keepdims=True), B / B.sum(1, keepdims=True)
    return 50 * np.abs(pa - pb).sum(1)


def agreement_inh(ta, tb, inh_types):
    """Figure 3c's per-donor agreement of inhibitory composition (fig3_data.py)."""
    shared = sorted(set(ta.index) & set(tb.index))
    inh = [c for c in inh_types if c in ta.columns and c in tb.columns]
    A = ta.loc[shared, inh].values.astype(float)
    B = tb.loc[shared, inh].values.astype(float)
    ok = (A.sum(1) >= MIN_CLASS_CELLS) & (B.sum(1) >= MIN_CLASS_CELLS)
    keep = ((A[ok] / A[ok].sum(1, keepdims=True)).mean(0) >= MIN_SHARE) & \
           ((B[ok] / B[ok].sum(1, keepdims=True)).mean(0) >= MIN_SHARE)
    za, zb = rs.clr(A[ok][:, keep]), rs.clr(B[ok][:, keep])
    za = (za - za.mean(0)) / za.std(0, ddof=1)
    zb = (zb - zb.mean(0)) / zb.std(0, ddof=1)
    r = [np.corrcoef(za[i], zb[i])[0, 1] for i in range(len(za))]
    icc = [2 * np.cov(za[i], zb[i], ddof=1)[0, 1] / (za[i].var(ddof=1) + zb[i].var(ddof=1))
           for i in range(len(za))]
    return pd.DataFrame({"donor": np.array(shared)[ok], "n_inh_types": int(keep.sum()),
                         "agree_inh": r, "icc_inh": icc})


def main():
    tax = rs.taxonomy()
    inh_types = list(tax.supertype_label[tax.class_id == "Inhibitory"])
    d = pd.read_csv(HERE / "data/fig8_depth_pair_donor.tsv", sep="\t")
    d["donor"] = d.donor.astype(str)
    marks, sup, out = {}, {}, []
    for a, b, pair in rs.PAIRS:
        for n in (a, b):
            if n not in marks:
                marks[n], sup[n] = marker_counts(n), rs.load(n)
        g = d[d.pair == pair].copy()
        assert len(g), f"{pair}: no rows in fig8_depth_pair_donor.tsv"
        for side, name in (("a", a), ("b", b)):
            m = marks[name].reindex(g.donor)
            assert m.n_l23.notna().all(), f"{pair} {name}: donors without marker counts"
            for col in ("n_l23", "n_l6it", "upper_frac"):
                g[f"{col}_{side}"] = m[col].values
        g["abs_d_upper"] = (g.upper_frac_a - g.upper_frac_b).abs()
        g["tv_inh"] = tv_inh(sup[a], sup[b], g.donor.tolist(), inh_types)
        g = g.merge(agreement_inh(sup[a], sup[b], inh_types), on="donor", how="left")
        print(f"  {pair:30s} n={len(g)}  median upper_frac A {g.upper_frac_a.median():.2f} "
              f"B {g.upper_frac_b.median():.2f}  median |diff| {g.abs_d_upper.median():.3f}  "
              f"agree_inh for {g.agree_inh.notna().sum()} ({int(g.n_inh_types.max())} inh types)")
        out.append(g)
    res = pd.concat(out, ignore_index=True)
    res.to_csv(HERE / "data/fig8_simple_depth.tsv", sep="\t", index=False)
    print(f"wrote data/fig8_simple_depth.tsv ({len(res)} rows)")


if __name__ == "__main__":
    main()
