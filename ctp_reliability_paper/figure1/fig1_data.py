"""Data for Figure 1: supertype-level proportions are more reliable than subclass.

Three lineages, DFC 2026 taxonomy: the subclass that contains microglia
("Immune"), Sst, and L2/3 IT.

Metric (decision 2026-10-05): symmetric ICC = 2 cov / (var_a + var_b) of
donor-level CLR between the two datasets of a pair, from
common/reliability_stats.py (rho, r and icc_person_a carried alongside).
Pair set and data source (decision 2026-10-29): the eight pairs of
reliability_stats.PAIRS, from the registry composition counts.

  (a) per-pair icc_sym (plus rho, r, icc_person_a) for each subclass and each of
      its supertypes across the eight matched-donor pairs, from
      common/data/pair_stats_percelltype.tsv.
  (c, d) donor-level CLR for Green x Mathys (panel c), BrainSCOPE CMC x PsychAD
      MSBB (panel d), plus Multiome x PsychAD ROSMAP (kept for reference), at
      both levels, from the same donors and CLR as reliability_stats
      (reliability_stats.load / clr_blocks):
        - registry count tables, shared donors = person_id intersection
        - a donor is kept in a compartment if it has >0 nuclei of that
          compartment in both datasets
        - CLR within compartment (neuronal vs non-neuronal), pseudocount 0.5
      icc_sym and n_donors are recomputed here and checked against
      common/data/pair_stats_percelltype.tsv; the script stops on any difference.
      (PsychAD x Ruzicka, formerly kept for reference, is gone: Ruzicka is not
      carried forward in the registry and has no count table.)

Reads /scratch/shreejoy (read-only). Writes data/ next to this script.
Run (login node is fine, ~1 min):
  PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig1_data.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import reliability_stats as rs  # noqa: E402

TAX = rs.TAX
OUT = Path(__file__).resolve().parent / "data"
COMMON = Path(__file__).resolve().parent.parent / "common" / "data"

SUBCLASSES = ["Immune", "Sst", "L2/3 IT"]
CLR_PAIRS = ["Green x Mathys", "Multiome x PsychAD ROSMAP", "BrainSCOPE CMC x PsychAD MSBB"]
PAIRS = [p for p in rs.PAIRS if p[2] in CLR_PAIRS]   # (dataset a, dataset b, pair label)


def pair_clr(ca, cb, level, pair, keep_types):
    long, icc = [], []
    for grp, types, donors, yA, yB in rs.clr_blocks(ca, cb, level):
        for j, ct in enumerate(types):
            a, b = yA[:, j], yB[:, j]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < rs.MIN_DONORS:
                continue
            ac, bc = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
            icc.append(dict(pair=pair, level=level, cell_type=ct, n_donors=int(ok.sum()),
                            icc_sym=float(2 * np.cov(ac, bc, ddof=1)[0, 1]
                                          / (ac.var(ddof=1) + bc.var(ddof=1)))))
            if ct in keep_types:
                long.append(pd.DataFrame(dict(pair=pair, level=level, cell_type=ct,
                                              group=grp, donor=donors[ok],
                                              clr_a=a[ok], clr_b=b[ok])))
    return pd.concat(long), pd.DataFrame(icc)


def check_against_common(icc, level, pair):
    ref = pd.read_csv(COMMON / "pair_stats_percelltype.tsv", sep="\t")
    ref = ref[(ref.pair == pair) & (ref.level == level)][["cell_type", "icc_sym", "n_donors"]]
    m = icc.merge(ref, on="cell_type", how="outer", indicator=True, suffixes=("", "_ref"))
    if (m._merge != "both").any():
        raise SystemExit(f"{pair} {level}: type sets differ from pair_stats:\n"
                         f"{m[m._merge != 'both']}")
    d = (m.icc_sym - m.icc_sym_ref).abs().max()
    if d > 1e-9 or (m.n_donors != m.n_donors_ref).any():
        raise SystemExit(f"{pair} {level}: ICC does not reproduce pair_stats (max |diff| = {d})")
    print(f"  {pair} {level}: {len(icc)} types, ICC reproduces pair_stats (max |diff| = {d:.1e})")


def main():
    OUT.mkdir(exist_ok=True)
    tax = pd.read_csv(TAX, sep="\t")
    tax = tax[tax.subclass_label.isin(SUBCLASSES)]
    supertypes = set(tax.supertype_label)
    st2sc = dict(zip(tax.supertype_label, tax.subclass_label))

    # ---- (a) per-pair statistics from the shared table (independent pairs) ---
    a = pd.read_csv(COMMON / "pair_stats_percelltype.tsv", sep="\t")
    a = a[a.subclass.isin(SUBCLASSES)]
    a = a[["level", "subclass", "cell_type", "pair", "n_donors", "icc_sym", "rho", "r",
           "icc_person_a"]]
    n_pairs = a.pair.nunique()
    print(f"(a) {n_pairs} matched-donor pairs: {sorted(a.pair.unique())}")
    print(a.groupby(["level", "subclass"]).cell_type.nunique().to_string())
    a.to_csv(OUT / "fig1a_pairs.tsv", sep="\t", index=False)

    # ---- (c, d) donor-level CLR ------------------------------------------------
    keep_types = set(SUBCLASSES) | supertypes
    cache = {}
    for name in sorted({x for p in PAIRS for x in p[:2]}):
        cache[name] = rs.load(name)
        print(f"{name}: {len(cache[name])} donors, {int(cache[name].values.sum()):,} nuclei")

    long = []
    for a_name, b_name, pair in PAIRS:
        ca, cb = cache[a_name], cache[b_name]
        print(f"{pair}: shared donors {len(set(ca.index) & set(cb.index))}")
        for level in ("subclass", "supertype"):
            lg, icc = pair_clr(ca, cb, level, pair, keep_types)
            check_against_common(icc, level, pair)
            long.append(lg)
    long = pd.concat(long)
    long.to_csv(OUT / "fig1bc_clr.tsv", sep="\t", index=False)
    print(f"wrote {OUT}/fig1a_pairs.tsv, fig1bc_clr.tsv "
          f"({long.groupby(['pair', 'level']).donor.nunique().to_dict()})")


if __name__ == "__main__":
    main()
