"""Data for Figure 1: supertype-level proportions are more reliable than subclass.

Three lineages, DFC 2026 taxonomy: the subclass that contains microglia
("Immune"), Sst, and L2/3 IT.

Metric (decision 2026-10-05): symmetric ICC = 2 cov / (var_a + var_b) of
donor-level CLR between the two datasets of a pair, from
common/reliability_stats.py (recomputed from the donor caches; rho, r and
icc_person_a are carried alongside and r / icc_sym / icc_person_a were checked
against 114 there).

  (a) per-pair icc_sym (plus rho, r, icc_person_a) for each subclass and each of
      its supertypes across the five independent matched-donor pairs, from
      common/data/pair_stats_percelltype.tsv.
  (c, d) donor-level CLR for Green x Mathys (panel c), BrainSCOPE CMC x PsychAD
      (panel d), plus Multiome x PsychAD and PsychAD x Ruzicka (kept for
      reference; Ruzicka's ~500 non-neuronal nuclei per donor make rare
      supertypes unmeasurable and the pair is excluded from the five),
      at both levels, reproducing wp3lib.decompose exactly as fig2a_clr_data.py:
        - person-keyed DFC caches, shared donors = ID intersection
        - a donor is kept in a compartment if it has >0 nuclei of that
          compartment in both datasets
        - CLR within compartment (neuronal vs non-neuronal), pseudocount 0.5
        - r = Pearson correlation of the two CLR vectors (r_cross_study of 114)
      Both r and icc_person_a are recomputed and checked against 114; the
      script stops on any difference.

Reads /scratch/shreejoy (read-only). Writes data/ next to this script.
Run (login node is fine, ~1 min):
  PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig1_data.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
RESULTS = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/"
           "refs/taxonomy_DFC_2026.tsv")
OUT = Path(__file__).resolve().parent / "data"
COMMON = Path(__file__).resolve().parent.parent / "common" / "data"

SUBCLASSES = ["Immune", "Sst", "L2/3 IT"]
PAIRS = [  # (dataset a, dataset b, label in 114) -- same order as 114_pairs_dfc.py
    ("Green_DFC_person", "Mathys_DFC_person", "Green x Mathys"),
    ("Multiome2025_DLPFC_DFC_person", "PsychAD_full_DFC_person", "Multiome x PsychAD"),
    ("BrainSCOPE_CMC_DFC_person", "PsychAD_full_DFC_person", "BrainSCOPE CMC x PsychAD"),
    ("PsychAD_full_DFC_person", "Ruzicka_full_DFC_person",
     "PsychAD x Ruzicka (MSSM2 x MSSM1)"),
]
PSEUDO = 0.5
NEURONAL_PREFIX = "Neuronal"
MIN_DONORS = 40  # MIN_SHARED in 114_pairs_dfc.py


def count_grid(cells, level):
    grp = np.where(cells.cls.str.startswith(NEURONAL_PREFIX), "neu", "glia")
    t = pd.crosstab(cells.donor, [grp, cells[level]])
    return {g: t[g] for g in ("neu", "glia") if g in t.columns.get_level_values(0)}


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def pair_clr(ca, cb, level, pair, keep_types):
    shared = sorted(set(ca.donor) & set(cb.donor))
    gA, gB = count_grid(ca, level), count_grid(cb, level)
    long, icc = [], []
    for grp in sorted(set(gA) & set(gB)):
        types = [c for c in gA[grp].columns if c in gB[grp].columns]
        dA = gA[grp].reindex(index=shared, columns=types).fillna(0).values
        dB = gB[grp].reindex(index=shared, columns=types).fillna(0).values
        keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
        dA, dB = dA[keep], dB[keep]
        donors = np.array(shared)[keep]
        yA, yB = clr(dA), clr(dB)
        for j, ct in enumerate(types):
            a, b = yA[:, j], yB[:, j]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < MIN_DONORS:
                continue
            ac, bc = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
            icc.append(dict(pair=pair, level=level, cell_type=ct, group=grp,
                            n_donors=int(ok.sum()),
                            r=float(np.corrcoef(ac, bc)[0, 1]),
                            icc=float(np.cov(ac, bc, ddof=1)[0, 1] / ac.var(ddof=1))))
            if ct in keep_types:
                long.append(pd.DataFrame(dict(pair=pair, level=level, cell_type=ct,
                                              group=grp, donor=donors[ok],
                                              clr_a=a[ok], clr_b=b[ok])))
    return pd.concat(long), pd.DataFrame(icc)


def check_against_114(icc, level, pair):
    ref = pd.read_csv(RESULTS / f"114_pair_percelltype_{level}.csv")
    ref = ref[ref.pair == pair][["cell_type", "r_cross_study", "icc_person_a", "n_donors"]]
    m = icc.merge(ref, on="cell_type", how="outer", indicator=True,
                  suffixes=("", "_114"))
    if (m._merge != "both").any():
        raise SystemExit(f"{pair} {level}: type sets differ from 114:\n"
                         f"{m[m._merge != 'both']}")
    d = max((m.r - m.r_cross_study).abs().max(), (m.icc - m.icc_person_a).abs().max())
    if d > 1e-9 or (m.n_donors != m.n_donors_114).any():
        raise SystemExit(f"{pair} {level}: r/ICC do not reproduce 114 (max |diff| = {d})")
    print(f"  {pair} {level}: {len(icc)} types, r and ICC reproduce 114 (max |diff| = {d:.1e})")


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
    print(f"(a) {n_pairs} independent pairs: {sorted(a.pair.unique())}")
    print(a.groupby(["level", "subclass"]).cell_type.nunique().to_string())
    a.to_csv(OUT / "fig1a_pairs.tsv", sep="\t", index=False)

    # ---- (b, c) donor-level CLR ------------------------------------------------
    keep_types = set(SUBCLASSES) | supertypes
    cols = ["donor", "supertype", "subclass", "cls"]
    cache = {}
    for name in sorted({x for p in PAIRS for x in p[:2]}):
        c = pd.read_parquet(CACHE / f"{name}_cells.parquet", columns=cols)
        c["donor"] = c["donor"].astype(str)
        cache[name] = c
        print(f"{name}: {c.donor.nunique()} donors, {len(c):,} nuclei")

    long = []
    for a_name, b_name, pair in PAIRS:
        ca, cb = cache[a_name], cache[b_name]
        print(f"{pair}: shared donors {len(set(ca.donor) & set(cb.donor))}")
        for level in ("subclass", "supertype"):
            lg, icc = pair_clr(ca, cb, level, pair, keep_types)
            check_against_114(icc, level, pair)
            long.append(lg)
    long = pd.concat(long)
    long.to_csv(OUT / "fig1bc_clr.tsv", sep="\t", index=False)
    print(f"wrote {OUT}/fig1a_pairs.tsv, fig1bc_clr.tsv "
          f"({long.groupby(['pair', 'level']).donor.nunique().to_dict()})")


if __name__ == "__main__":
    main()
