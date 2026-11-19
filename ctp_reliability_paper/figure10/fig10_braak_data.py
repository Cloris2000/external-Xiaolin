"""Sst_25 vs Braak stage (AD severity), Mount Sinai (MSBB) cohorts, bulk and snRNA-seq.

Does the Sst_25 decrease track AD neuropathological severity, not only diagnosis?

Donors: AD continuum only -- registry dx in {Control, MCI, AD}; other neurodegenerative
diagnoses (FTD, PD, OtherDementia, ...) are left out because their Braak stage reflects a
different process. All ages (as the Fig 10 AD analysis). Braak = registry covariates.tsv
`braak` (0-6, harmonised from PsychAD/MSBB). Complete age, sex, PMI required.

Cohorts:
  bulk   MSBB_FP       production CelMod Sst_25 prediction (no CelMod training donors here)
  snRNA  PsychAD_MSBB  within-compartment CLR (fig10_data.snrna_matrix)
  same   people in both, each modality z-scored over the shared donors

Model per cohort/modality (outcome z-scored over its analysis donors, SD units):
  z ~ braak + age + sex + PMI     beta = change in Sst_25 per Braak stage (linear trend)
Same-donor difference test: (z_bulk - z_snRNA) ~ braak + age + sex + PMI, whose braak
coefficient = beta_bulk - beta_snRNA (do the two modalities track Braak equally?).

Writes data/fig10_braak_donors.tsv and data/fig10_braak_effects.tsv.
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig10_braak_data.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fig10_data as F   # noqa: E402

TYPE = "Sst_25"
DX = ["Control", "MCI", "AD"]
BULK, SNRNA = "MSBB_FP", "PsychAD_MSBB"


def zs(v):
    return (v - v.mean()) / v.std(ddof=1)


def fit(y, m, adj_dx=False):
    """coefficient on column 1 = braak. adj_dx adds AD and MCI indicators (vs Control), so the
    Braak slope is severity beyond diagnosis."""
    extra = [(m.dx == "AD").astype(float).values, (m.dx == "MCI").astype(float).values] if adj_dx else []
    M = np.column_stack([np.ones(len(m)), m.braak.values] + [m[c].values for c in F.COVARS] + extra)
    return F.ols_case(y, M)


def main():
    reg = pd.read_csv(F.COV, sep="\t", low_memory=False, dtype={"person_id": str}).set_index("person_id")
    C = F.covariates()
    C["braak"] = pd.to_numeric(reg.braak, errors="coerce")
    assert C.braak.dropna().isin(range(7)).all(), "Braak outside 0-6"
    training = set(pd.read_csv(F.DECONV / "prod_labels.csv").donor.astype(str))

    Xb, n_tr = F.bulk_matrix(BULK, training)
    Xs = F.snrna_matrix(SNRNA)
    print(f"{BULK}: {n_tr} CelMod training donors removed")

    def donors(X):
        m = C.reindex(X.index)
        m = m[m.dx.isin(DX)].dropna(subset=["braak"] + F.COVARS)
        y = X.loc[m.index, TYPE].astype(float)
        return m[y.notna().values], y[y.notna()]

    don, eff = [], []
    sets = {}
    for mod, X, cohort in (("bulk", Xb, BULK), ("snRNA", Xs, SNRNA)):
        m, y = donors(X)
        sets[mod] = (m, X)
        b, se, p = fit(zs(y.values), m)
        eff.append(dict(set="all", modality=mod, cohort=cohort, n=len(m), beta=b, se=se, p=p,
                        **{f"n_{d}": int((m.dx == d).sum()) for d in DX}))
        b2, se2, p2 = fit(zs(y.values), m, adj_dx=True)
        eff[-1].update(beta_adj_dx=b2, se_adj_dx=se2, p_adj_dx=p2)
        print(f"      {mod:5s} + dx adjustment: beta/Braak stage {b2:+.3f} (SE {se2:.3f}) p {p2:.2g}")
        don.append(pd.DataFrame(dict(set="all", modality=mod, cohort=cohort, person_id=m.index,
                                     dx=m.dx.values, braak=m.braak.values, z=zs(y.values))))
        print(f"all   {mod:5s} {cohort:13s} n={len(m):4d} "
              f"({', '.join(f'{d} {(m.dx == d).sum()}' for d in DX)})  "
              f"beta/Braak stage {b:+.3f} (SE {se:.3f}) p {p:.2g}")

    # same donors
    ids = sorted(set(sets["bulk"][0].index) & set(sets["snRNA"][0].index))
    m = C.loc[ids]
    z = {mod: zs(sets[mod][1].loc[ids, TYPE].astype(float).values) for mod in ("bulk", "snRNA")}
    for mod, cohort in (("bulk", BULK), ("snRNA", SNRNA)):
        b, se, p = fit(z[mod], m)
        eff.append(dict(set="same_donor", modality=mod, cohort=cohort, n=len(ids), beta=b, se=se, p=p,
                        **{f"n_{d}": int((m.dx == d).sum()) for d in DX}))
        b2, se2, p2 = fit(z[mod], m, adj_dx=True)
        eff[-1].update(beta_adj_dx=b2, se_adj_dx=se2, p_adj_dx=p2)
        don.append(pd.DataFrame(dict(set="same_donor", modality=mod, cohort=cohort, person_id=ids,
                                     dx=m.dx.values, braak=m.braak.values, z=z[mod])))
        print(f"same  {mod:5s} {cohort:13s} n={len(ids):4d}  beta/Braak stage {b:+.3f} (SE {se:.3f}) p {p:.2g}"
              f"  | + dx {b2:+.3f} p {p2:.2g}")
    b, se, p = fit(z["bulk"] - z["snRNA"], m)
    assert np.isclose(b, eff[-2]["beta"] - eff[-1]["beta"])
    eff.append(dict(set="same_donor", modality="bulk_minus_snRNA", cohort=f"{BULK} x {SNRNA}",
                    n=len(ids), beta=b, se=se, p=p, **{f"n_{d}": int((m.dx == d).sum()) for d in DX}))
    r = float(pd.Series(z["bulk"]).corr(pd.Series(z["snRNA"]), method="spearman"))
    print(f"same  bulk - snRNA difference in Braak slope {b:+.3f} (SE {se:.3f}) p {p:.2g}; "
          f"per-donor Spearman rho {r:.2f}")

    out = HERE / "data"
    pd.concat(don).to_csv(out / "fig10_braak_donors.tsv", sep="\t", index=False)
    pd.DataFrame(eff).to_csv(out / "fig10_braak_effects.tsv", sep="\t", index=False)
    print("wrote data/fig10_braak_donors.tsv, data/fig10_braak_effects.tsv")


if __name__ == "__main__":
    main()
