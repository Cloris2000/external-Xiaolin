"""Per-pair reliability statistics recomputed from the donor-level CLR values,
for the independent matched-donor pairs, at subclass and supertype level.

Pair set (decision 2026-10-05): five pairs. PsychAD x Ruzicka is excluded:
Ruzicka has ~500 non-neuronal nuclei per donor, so rare supertypes sit on the
pseudocount floor, and its 51 donors are a subset of the 54 in BrainSCOPE CMC x
PsychAD, so it was not independent of that pair either.

Main-text metric (decision 2026-10-05): the symmetric ICC
    icc_sym = 2 cov(a, b) / (var(a) + var(b))
i.e. the agreement (one-way) ICC of two interchangeable measurements of the same
donors, identical to Lin's concordance coefficient on centred data. Unlike
Shreejoy's icc_person_a = cov / var(a) it does not depend on which dataset was
listed first in the pair. Pearson r, Spearman rho and icc_person_a are carried
alongside for the supplement. All four are computed from the same CLR vectors,
rebuilt exactly as wp3lib.decompose does (see fig2a_clr_data.py); r, icc_sym
and icc_person_a are checked against 114 for every pair and level, which proves
the donors and the CLR are the ones behind the published numbers.

Inputs (read-only): /scratch/shreejoy/cell_type_bias/wp3/
  cache/<dataset>_cells.parquet      per-nucleus donor / supertype / subclass / class
  results/114_pair_percelltype_{subclass,supertype}.csv   for the check
Outputs (data/ next to this script):
  pair_stats_percelltype.tsv   per pair x level x cell type: n_donors, icc_sym, rho, r,
                               icc_person_a
  reliability_summary.tsv      per level x cell type: n_pairs; median/IQR/min of icc_sym;
                               medians of rho, r, icc_person_a
                               (types measured in >= MIN_PAIRS pairs, as 117_whitelist)
Run (login node is fine, ~1 min wall, ~6 GB RAM for the PsychAD cache):
  PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python reliability_stats.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
RESULTS = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/"
           "refs/taxonomy_DFC_2026.tsv")
OUT = Path(__file__).resolve().parent / "data"

# (dataset a, dataset b, label) in the order of 114_pairs_dfc.py; a is the 114 reference
PAIRS = [
    ("Green_DFC_person", "Mathys_DFC_person", "Green x Mathys"),
    ("Green_DFC_person", "Multiome2025_DLPFC_DFC_person", "Green x Multiome"),
    ("Mathys_DFC_person", "Multiome2025_DLPFC_DFC_person", "Mathys x Multiome"),
    ("Multiome2025_DLPFC_DFC_person", "PsychAD_full_DFC_person", "Multiome x PsychAD"),
    ("BrainSCOPE_CMC_DFC_person", "PsychAD_full_DFC_person", "BrainSCOPE CMC x PsychAD"),
    # ("PsychAD_full_DFC_person", "Ruzicka_full_DFC_person",
    #  "PsychAD x Ruzicka (MSSM2 x MSSM1)"),   # excluded 2026-10-05, see header
]
PSEUDO = 0.5
NEURONAL_PREFIX = "Neuronal"
MIN_DONORS = 40   # MIN_SHARED in 114_pairs_dfc.py
MIN_PAIRS = 2     # as 117_whitelist_dfc.py


def count_grid(cells, level):
    grp = np.where(cells.cls.str.startswith(NEURONAL_PREFIX), "neu", "glia")
    t = pd.crosstab(cells.donor, [grp, cells[level]])
    return {g: t[g] for g in ("neu", "glia") if g in t.columns.get_level_values(0)}


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def pair_stats(ca, cb, level, pair):
    shared = sorted(set(ca.donor) & set(cb.donor))
    gA, gB = count_grid(ca, level), count_grid(cb, level)
    rows = []
    for grp in sorted(set(gA) & set(gB)):
        types = [c for c in gA[grp].columns if c in gB[grp].columns]
        dA = gA[grp].reindex(index=shared, columns=types).fillna(0).values
        dB = gB[grp].reindex(index=shared, columns=types).fillna(0).values
        keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
        yA, yB = clr(dA[keep]), clr(dB[keep])
        for j, ct in enumerate(types):
            a, b = yA[:, j], yB[:, j]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < MIN_DONORS:
                continue
            ac, bc = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
            cov = float(np.cov(ac, bc, ddof=1)[0, 1])
            va, vb = float(ac.var(ddof=1)), float(bc.var(ddof=1))
            rows.append(dict(pair=pair, level=level, cell_type=ct, group=grp,
                             n_donors=int(ok.sum()),
                             icc_sym=2 * cov / (va + vb),
                             rho=float(spearmanr(ac, bc).correlation),
                             r=float(np.corrcoef(ac, bc)[0, 1]),
                             icc_person_a=cov / va))
    return pd.DataFrame(rows)


def check_against_114(st, level, pair):
    ref = pd.read_csv(RESULTS / f"114_pair_percelltype_{level}.csv")
    ref = ref[ref.pair == pair][["cell_type", "r_cross_study", "icc_person_a", "n_donors",
                                 "v_person", "v_total_a", "v_total_b"]]
    ref["icc_sym_114"] = 2 * ref.v_person / (ref.v_total_a + ref.v_total_b)
    m = st.merge(ref, on="cell_type", how="outer", indicator=True, suffixes=("", "_114"))
    if (m._merge != "both").any():
        raise SystemExit(f"{pair} {level}: type sets differ from 114:\n{m[m._merge != 'both']}")
    d = max((m.r - m.r_cross_study).abs().max(),
            (m.icc_person_a - m.icc_person_a_114).abs().max(),
            (m.icc_sym - m.icc_sym_114).abs().max())
    if d > 1e-9 or (m.n_donors != m.n_donors_114).any():
        raise SystemExit(f"{pair} {level}: r/ICC do not reproduce 114 (max |diff| = {d})")
    print(f"  {pair} {level}: {len(st)} types, r, icc_sym and icc_person_a reproduce 114 "
          f"(max |diff| = {d:.1e})")


def load(name):
    c = pd.read_parquet(CACHE / f"{name}_cells.parquet",
                        columns=["donor", "supertype", "subclass", "cls"])
    c["donor"] = c["donor"].astype(str)
    return c


def main():
    OUT.mkdir(exist_ok=True)
    tax = pd.read_csv(TAX, sep="\t")
    st2sc = dict(zip(tax.supertype_label, tax.subclass_label))

    cache, out = {}, []
    for a_name, b_name, pair in PAIRS:
        for n in (a_name, b_name):
            if n not in cache:
                cache[n] = load(n)
                print(f"{n}: {cache[n].donor.nunique()} donors, {len(cache[n]):,} nuclei")
        ca, cb = cache[a_name], cache[b_name]
        print(f"{pair}: shared donors {len(set(ca.donor) & set(cb.donor))}")
        for level in ("subclass", "supertype"):
            st = pair_stats(ca, cb, level, pair)
            check_against_114(st, level, pair)
            out.append(st)
        # free datasets no later pair needs
        later = {x for p in PAIRS[PAIRS.index((a_name, b_name, pair)) + 1:] for x in p[:2]}
        for n in list(cache):
            if n not in later:
                del cache[n]

    per = pd.concat(out, ignore_index=True)
    per["subclass"] = np.where(per.level == "subclass", per.cell_type, per.cell_type.map(st2sc))
    per = per[["pair", "level", "subclass", "cell_type", "group", "n_donors",
               "icc_sym", "rho", "r", "icc_person_a"]]
    per.to_csv(OUT / "pair_stats_percelltype.tsv", sep="\t", index=False)

    keys = ["level", "subclass", "cell_type", "group"]
    g = per.groupby(keys).icc_sym
    s = pd.DataFrame({"n_pairs": g.size(), "median_icc_sym": g.median(),
                      "q25": g.quantile(0.25), "q75": g.quantile(0.75), "min_icc_sym": g.min(),
                      "n_pairs_ge_030": g.apply(lambda x: int((x >= 0.30).sum()))})
    for c in ("rho", "r", "icc_person_a"):
        s[f"median_{c}"] = per.groupby(keys)[c].median()
    s = s.reset_index()
    s = s[s.n_pairs >= MIN_PAIRS].sort_values(["level", "median_icc_sym"], ascending=[True, False])
    s.to_csv(OUT / "reliability_summary.tsv", sep="\t", index=False)

    for level in ("subclass", "supertype"):
        x = s[s.level == level]
        print(f"{level}: {len(x)} types in >= {MIN_PAIRS} pairs; median ICC {x.median_icc_sym.median():.3f}; "
              f"{(x.median_icc_sym >= 0.30).sum()} with median ICC >= 0.30")
    print("top 10 supertypes by median symmetric ICC:")
    print(s[s.level == "supertype"].head(10)[["cell_type", "n_pairs", "median_icc_sym", "q25", "q75",
                                              "median_rho", "median_r"]].round(3).to_string(index=False))
    print(f"wrote {OUT}/pair_stats_percelltype.tsv, reliability_summary.tsv")


if __name__ == "__main__":
    main()
