"""Donor-bootstrap confidence intervals for the symmetric ICC (Fig 4 'Core Supertype').

Same donors, CLR vectors and icc_sym as reliability_stats.py (imported, not
duplicated; point estimates are checked against data/pair_stats_percelltype.tsv).
For every pair x level x cell type the matched donors are resampled with
replacement N_BOOT times (the same donor draw is used for all cell types of a
pair, so the cell types stay jointly resampled), icc_sym is recomputed, and
  - per-pair 95% CI = 2.5/97.5% quantiles
  - pooled 95% CI of the median ICC across pairs = quantiles of the per-resample
    median over the pairs in which the type is measured (same point estimate as
    Fig 1a / reliability_summary.tsv). Pairs are resampled independently, so the
    pooled CI ignores donor overlap between pairs (slightly narrow).
Outputs (RELIABILITY_CI_OUT, default data/ next to this script; /project is
read-only on Trillium compute nodes, so the sbatch points it at scratch):
  ci_pair_percelltype.tsv      pair, level, cell_type, n_donors, icc_sym, ci_lo, ci_hi
  reliability_ci_summary.tsv   level, cell_type, n_pairs, median_icc_sym, ci_lo, ci_hi,
                               n_pairs_ci_gt0
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

import reliability_stats as rs

N_BOOT = 1000
CHUNK = 100
SEED = 20261027
OUT = Path(os.environ.get("RELIABILITY_CI_OUT", Path(__file__).resolve().parent / "data"))


def boot_icc(yA, yB, rng):
    """Bootstrap icc_sym for all columns at once; returns (N_BOOT, T)."""
    n = yA.shape[0]
    aa, bb, ab = yA * yA, yB * yB, yA * yB
    out = []
    for _ in range(N_BOOT // CHUNK):
        w = rng.multinomial(n, np.full(n, 1.0 / n), size=CHUNK) / n   # (CHUNK, n) resample weights
        ma, mb = w @ yA, w @ yB
        cov = w @ ab - ma * mb
        va, vb = w @ aa - ma ** 2, w @ bb - mb ** 2
        out.append(2 * cov / (va + vb))
    return np.vstack(out)


def pair_boot(ca, cb, level, pair, rng):
    shared = sorted(set(ca.index) & set(cb.index))
    gA, gB = rs.count_grid(ca, level), rs.count_grid(cb, level)
    rows, boots = [], []
    for grp in sorted(set(gA) & set(gB)):
        types = [c for c in gA[grp].columns if c in gB[grp].columns]
        dA = gA[grp].reindex(index=shared, columns=types).fillna(0).values
        dB = gB[grp].reindex(index=shared, columns=types).fillna(0).values
        keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
        yA, yB = rs.clr(dA[keep]), rs.clr(dB[keep])
        if keep.sum() < rs.MIN_DONORS:
            continue
        bt = boot_icc(yA, yB, rng)
        for j, ct in enumerate(types):
            a, b = yA[:, j], yB[:, j]
            ac, bc = a - a.mean(), b - b.mean()
            icc = 2 * np.cov(ac, bc, ddof=1)[0, 1] / (ac.var(ddof=1) + bc.var(ddof=1))
            rows.append(dict(pair=pair, level=level, cell_type=ct, n_donors=int(keep.sum()),
                             icc_sym=icc, ci_lo=np.quantile(bt[:, j], 0.025),
                             ci_hi=np.quantile(bt[:, j], 0.975)))
            boots.append(bt[:, j])
    return pd.DataFrame(rows), np.array(boots)   # boots: (types, N_BOOT)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    cache, per, boot = {}, [], {}
    for a_name, b_name, pair in rs.PAIRS:
        for n in (a_name, b_name):
            if n not in cache:
                cache[n] = rs.load(n)
        for level in ("subclass", "supertype"):
            df, bt = pair_boot(cache[a_name], cache[b_name], level, pair, rng)
            per.append(df)
            for ct, b in zip(df.cell_type, bt):
                boot.setdefault((level, ct), []).append(b)
            print(f"{pair} {level}: {len(df)} types", flush=True)
        later = {x for p in rs.PAIRS[rs.PAIRS.index((a_name, b_name, pair)) + 1:] for x in p[:2]}
        for n in list(cache):
            if n not in later:
                del cache[n]
    per = pd.concat(per, ignore_index=True)

    # point estimates must equal the published per-pair table
    ref = pd.read_csv(Path(__file__).resolve().parent / "data" / "pair_stats_percelltype.tsv", sep="\t")
    m = per.merge(ref[["pair", "level", "cell_type", "icc_sym", "n_donors"]],
                  on=["pair", "level", "cell_type"], how="outer", indicator=True, suffixes=("", "_ref"))
    assert (m._merge == "both").all(), m[m._merge != "both"]
    assert np.abs(m.icc_sym - m.icc_sym_ref).max() < 1e-9 and (m.n_donors == m.n_donors_ref).all()
    print("point estimates reproduce pair_stats_percelltype.tsv")
    per.to_csv(OUT / "ci_pair_percelltype.tsv", sep="\t", index=False)

    rows = []
    for (level, ct), bl in boot.items():
        if len(bl) < rs.MIN_PAIRS:
            continue
        med = np.median(np.vstack(bl), axis=0)
        p = per[(per.level == level) & (per.cell_type == ct)]
        rows.append(dict(level=level, cell_type=ct, n_pairs=len(bl),
                         median_icc_sym=p.icc_sym.median(),
                         ci_lo=np.quantile(med, 0.025), ci_hi=np.quantile(med, 0.975),
                         n_pairs_ci_gt0=int((p.ci_lo > 0).sum())))
    s = pd.DataFrame(rows).sort_values(["level", "median_icc_sym"], ascending=[True, False])
    s.to_csv(OUT / "reliability_ci_summary.tsv", sep="\t", index=False)
    for level in ("subclass", "supertype"):
        x = s[s.level == level]
        print(f"{level}: {len(x)} types; pooled CI lower bound > 0 in {(x.ci_lo > 0).sum()}; "
              f"CI > 0 in all pairs in {(x.n_pairs_ci_gt0 == x.n_pairs).sum()}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
