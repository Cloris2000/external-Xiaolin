"""Donor-level CLR for Figure 2a: Green vs Mathys, DFC 2026 taxonomy.

Same donors and CLR as common/reliability_stats.py (reliability_stats.load /
clr_blocks; registry composition counts since 2026-10-29), so the scatter
shows exactly the data behind the per-pair ICC:
  - registry count tables, shared donors = person_id intersection
  - a donor is kept in a compartment if it has >0 nuclei of that compartment
    in both datasets
  - CLR within compartment (neuronal vs non-neuronal), pseudocount 0.5
  - r = Pearson correlation of the two CLR vectors; icc = cov(a, b) / var(a)
    (icc_person_a), kept for reference; the paper's metric is icc_sym

icc_sym and n_donors are recomputed here and checked against
common/data/pair_stats_percelltype.tsv; the script stops if any supertype differs.

Reads /scratch/shreejoy (read-only). Writes data/ next to this script.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import reliability_stats as rs  # noqa: E402

OUT = Path(__file__).resolve().parent / "data"
COMMON = Path(__file__).resolve().parent.parent / "common" / "data"
A, B, PAIR = "Green", "Mathys", "Green x Mathys"


def main():
    OUT.mkdir(exist_ok=True)
    ca, cb = rs.load(A), rs.load(B)
    print(f"{A}: {len(ca)} donors; {B}: {len(cb)} donors; "
          f"shared {len(set(ca.index) & set(cb.index))}")

    long, icc = [], []
    for grp, types, donors, yA, yB in rs.clr_blocks(ca, cb, "supertype"):
        for j, ct in enumerate(types):
            a, b = yA[:, j], yB[:, j]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < rs.MIN_DONORS:
                continue
            ac, bc = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
            cov = np.cov(ac, bc, ddof=1)[0, 1]
            icc.append(dict(cell_type=ct, group=grp, n_donors=int(ok.sum()),
                            r=float(np.corrcoef(ac, bc)[0, 1]),
                            icc=float(cov / ac.var(ddof=1)),
                            icc_sym=float(2 * cov / (ac.var(ddof=1) + bc.var(ddof=1)))))
            long.append(pd.DataFrame(dict(pair=PAIR, cell_type=ct, group=grp,
                                          donor=donors[ok], clr_a=a[ok], clr_b=b[ok])))

    icc = pd.DataFrame(icc)
    ref = pd.read_csv(COMMON / "pair_stats_percelltype.tsv", sep="\t")
    ref = ref[(ref.pair == PAIR) & (ref.level == "supertype")][["cell_type", "icc_sym", "n_donors"]]
    m = icc.merge(ref, on="cell_type", how="outer", indicator=True, suffixes=("", "_ref"))
    if (m._merge != "both").any():
        raise SystemExit(f"supertype sets differ from pair_stats:\n{m[m._merge != 'both']}")
    d = (m.icc_sym - m.icc_sym_ref).abs().max()
    if d > 1e-9 or (m.n_donors != m.n_donors_ref).any():
        raise SystemExit(f"ICC does not reproduce pair_stats (max |diff| = {d})")
    print(f"{len(icc)} supertypes; ICC reproduces pair_stats (max |diff| = {d:.1e})")

    pd.concat(long).to_csv(OUT / "fig2a_clr_green_mathys.tsv", sep="\t", index=False)
    icc.to_csv(OUT / "fig2a_icc_green_mathys.tsv", sep="\t", index=False)
    print(f"wrote {OUT}/fig2a_clr_green_mathys.tsv, fig2a_icc_green_mathys.tsv")


if __name__ == "__main__":
    main()
