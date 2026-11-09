"""Explainer for the Figure 8 y-axis ("% of inhibitory neurons that differ"),
for discussion with the PI, not a paper panel (2026-11-03).

Two real Green x Mathys donors from data/fig8_simple_depth.tsv, both passing the
Figure 8 filters:
  BDR000612  sampled at the same depth in both datasets (upper-layer share of
             excitatory neurons 0.91 vs 0.91); 17% of inhibitory neurons differ
  BDR000842  sampled at different depths (0.90 vs 0.33); 42% differ
Neither is the most extreme of its kind (matched donors range 15-19%, mismatched
34-57%).

Output: data/fig8_y_explainer.tsv, one row per donor x dataset x inhibitory
supertype: share of the donor's inhibitory nuclei, subclass, and the donor's
upper_frac and tv_inh, which are checked against fig8_simple_depth.tsv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import reliability_stats as rs  # noqa: E402

PAIR, A, B = "Green x Mathys", "Green", "Mathys"
DONORS = {"BDR000612": "Same sampled depth", "BDR000842": "Different sampled depth"}


def main():
    tax = rs.taxonomy()
    inh = tax[tax.class_id == "Inhibitory"]
    sc = dict(zip(inh.supertype_label, inh.subclass_label))
    d = pd.read_csv(HERE / "data/fig8_simple_depth.tsv", sep="\t")
    d = d[d.pair == PAIR].set_index("donor")
    rows = []
    for ds, side in ((A, "a"), (B, "b")):
        c = rs.load(ds)
        for donor, case in DONORS.items():
            x = c.loc[donor, [t for t in sc if t in c.columns]].astype(float)
            share = 100 * x / x.sum()
            rows.append(pd.DataFrame({"case": case, "donor": donor, "dataset": ds,
                                      "supertype": share.index,
                                      "subclass": share.index.map(sc), "share": share.values,
                                      "upper_frac": d.loc[donor, f"upper_frac_{side}"],
                                      "tv_inh": d.loc[donor, "tv_inh"]}))
    out = pd.concat(rows, ignore_index=True)
    # the bars must give back the plotted % exactly
    for donor, g in out.groupby("donor"):
        w = g.pivot(index="supertype", columns="dataset", values="share").fillna(0)
        tv = 0.5 * np.abs(w[A] - w[B]).sum()
        assert abs(tv - g.tv_inh.iloc[0]) < 1e-6, (donor, tv, g.tv_inh.iloc[0])
        print(f"{donor}: upper-layer share {g[g.dataset == A].upper_frac.iloc[0]:.2f} vs "
              f"{g[g.dataset == B].upper_frac.iloc[0]:.2f}; {tv:.1f}% of inhibitory neurons differ")
    out.to_csv(HERE / "data/fig8_y_explainer.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
