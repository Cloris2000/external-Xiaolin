"""Donor-level CLR for Figure 2a: Green vs Mathys, DFC 2026 taxonomy.

Reproduces wp3lib.decompose (114_pairs_dfc.py) for one pair so the scatter
shows exactly the data behind the published ICC:
  - person-keyed DFC caches (Green_DFC_person, Mathys_DFC_person)
  - shared donors only; a donor is kept in a compartment if it has >0 nuclei
    of that compartment in both datasets
  - CLR within compartment (neuronal vs non-neuronal), pseudocount 0.5
  - ICC = cov(a, b) / var(a), the icc_person_a of 114

The ICC recomputed here is checked against 114_pair_percelltype_supertype.csv
and the script stops if any supertype differs.

Reads /scratch/shreejoy (read-only). Writes data/ next to this script.
"""
from pathlib import Path

import numpy as np
import pandas as pd

CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
RESULTS = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
OUT = Path(__file__).resolve().parent / "data"

A, B, PAIR = "Green_DFC_person", "Mathys_DFC_person", "Green x Mathys"
PSEUDO = 0.5
NEURONAL_PREFIX = "Neuronal"
MIN_DONORS = 40  # MIN_SHARED in 114_pairs_dfc.py


def count_grid(cells):
    grp = np.where(cells.cls.str.startswith(NEURONAL_PREFIX), "neu", "glia")
    t = pd.crosstab(cells.donor, [grp, cells.supertype])
    return {g: t[g] for g in ("neu", "glia") if g in t.columns.get_level_values(0)}


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def main():
    OUT.mkdir(exist_ok=True)
    cols = ["donor", "supertype", "cls"]
    ca = pd.read_parquet(CACHE / f"{A}_cells.parquet", columns=cols)
    cb = pd.read_parquet(CACHE / f"{B}_cells.parquet", columns=cols)
    for c in (ca, cb):
        c["donor"] = c["donor"].astype(str)
    shared = sorted(set(ca.donor) & set(cb.donor))
    print(f"{A}: {ca.donor.nunique()} donors; {B}: {cb.donor.nunique()} donors; "
          f"shared {len(shared)}")

    gA, gB = count_grid(ca), count_grid(cb)
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
            icc.append(dict(cell_type=ct, group=grp, n_donors=int(ok.sum()),
                            icc=float(np.cov(ac, bc, ddof=1)[0, 1] / ac.var(ddof=1))))
            long.append(pd.DataFrame(dict(pair=PAIR, cell_type=ct, group=grp,
                                          donor=donors[ok], clr_a=a[ok], clr_b=b[ok])))

    icc = pd.DataFrame(icc)
    ref = pd.read_csv(RESULTS / "114_pair_percelltype_supertype.csv")
    ref = ref[ref.pair == PAIR][["cell_type", "icc_person_a", "n_donors"]]
    m = icc.merge(ref, on="cell_type", how="outer", indicator=True,
                  suffixes=("", "_114"))
    if (m._merge != "both").any():
        raise SystemExit(f"supertype sets differ from 114:\n{m[m._merge != 'both']}")
    d = (m.icc - m.icc_person_a).abs().max()
    if d > 1e-9 or (m.n_donors != m.n_donors_114).any():
        raise SystemExit(f"ICC does not reproduce 114 (max |diff| = {d})")
    print(f"{len(icc)} supertypes; ICC reproduces 114 (max |diff| = {d:.1e})")

    pd.concat(long).to_csv(OUT / "fig2a_clr_green_mathys.tsv", sep="\t", index=False)
    icc.to_csv(OUT / "fig2a_icc_green_mathys.tsv", sep="\t", index=False)
    print(f"wrote {OUT}/fig2a_clr_green_mathys.tsv, fig2a_icc_green_mathys.tsv")


if __name__ == "__main__":
    main()
