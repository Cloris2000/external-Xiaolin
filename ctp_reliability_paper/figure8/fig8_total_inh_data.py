"""Figure 8 variant (PI suggestion 2026-10-23): y = the TOTAL inhibitory share
of neurons, not the per-supertype average.

For each donor and dataset: f = n_inh / (n_inh + n_exc) over all neuronal
nuclei (class_id Inhibitory / Excitatory from taxonomy_DFC_2026.tsv). On the
logit scale with a 0.5 pseudocount on both counts (the two-part CLR up to a
factor of 2), each dataset's across-donor mean is removed, and

    disc_inh_total = |residual_A - residual_B|

Why: the inh:exc ratio varies far less across layers than the interneuron
subtype mix does, so it is a composition measure that is not mainly a depth
readout. If depth mismatch still predicts disagreement in it, the Figure 8
effect is not just "laminar subtypes measured twice". Same donors, pass filter
and x (wasserstein) as the main figure, from data/fig8_depth_pair_donor.tsv.

Output: data/fig8_total_inh.tsv = fig8_depth_pair_donor.tsv plus n_inh_*,
n_exc_*, inh_frac_* per side and disc_inh_total.
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig8_total_inh_data.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/"
           "refs/taxonomy_DFC_2026.tsv")
HERE = Path(__file__).resolve().parent
PAIRS = {
    "Green x Mathys": ("Green_DFC_person", "Mathys_DFC_person"),
    "Green x Multiome": ("Green_DFC_person", "Multiome2025_DLPFC_DFC_person"),
    "Mathys x Multiome": ("Mathys_DFC_person", "Multiome2025_DLPFC_DFC_person"),
    "Multiome x PsychAD": ("Multiome2025_DLPFC_DFC_person", "PsychAD_full_DFC_person"),
    "BrainSCOPE CMC x PsychAD": ("BrainSCOPE_CMC_DFC_person", "PsychAD_full_DFC_person"),
}
PSEUDO = 0.5


def class_counts(name, cls_of):
    c = pd.read_parquet(CACHE / f"{name}_cells.parquet", columns=["donor", "supertype"])
    c["donor"] = c.donor.astype(str)
    c["cls"] = c.supertype.map(cls_of)
    t = pd.crosstab(c.donor, c.cls).reindex(columns=["Inhibitory", "Excitatory"], fill_value=0)
    t.columns = ["n_inh", "n_exc"]
    t["inh_frac"] = t.n_inh / (t.n_inh + t.n_exc)
    t["logit"] = np.log((t.n_inh + PSEUDO) / (t.n_exc + PSEUDO))
    return t


def main():
    tax = pd.read_csv(TAX, sep="\t")
    cls_of = dict(zip(tax.supertype_label, tax.class_id))
    d = pd.read_csv(HERE / "data/fig8_depth_pair_donor.tsv", sep="\t")
    cache, out = {}, []
    for pair, (a, b) in PAIRS.items():
        for n in (a, b):
            if n not in cache:
                cache[n] = class_counts(n, cls_of)
        g = d[d.pair == pair].copy()
        shared = g.donor
        res = {}
        for side, name in (("a", a), ("b", b)):
            m = cache[name].reindex(shared)
            for col in ("n_inh", "n_exc", "inh_frac"):
                g[f"{col}_{side}"] = m[col].values
            # residual: this dataset's across-donor mean logit (over the pair's
            # shared donors) removed, as the per-supertype version does
            res[side] = m.logit.values - np.nanmean(m.logit.values)
        g["disc_inh_total"] = np.abs(res["a"] - res["b"])
        print(f"  {pair:26s} n={len(g)}  median inh_frac A {g.inh_frac_a.median():.2f} "
              f"B {g.inh_frac_b.median():.2f}  median disc_inh_total {g.disc_inh_total.median():.3f}")
        out.append(g)
    r = pd.concat(out, ignore_index=True)
    r.to_csv(HERE / "data/fig8_total_inh.tsv", sep="\t", index=False)
    print(f"wrote data/fig8_total_inh.tsv ({len(r)} rows)")


if __name__ == "__main__":
    main()
