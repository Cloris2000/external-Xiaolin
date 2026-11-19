"""Sst_25 supplement to Figure 10, same-donor comparison: people measured by BOTH a bulk and an
snRNA-seq cohort in the Fig 10 analysis sets (data/fig10_donors.tsv), so bulk and snRNA are
compared on identical donors, identical case/control labels and identical covariates.

Per (disease, bulk cohort x snRNA cohort) pair with >= MIN_PER_GROUP matched cases and controls:
  * outcome per modality = Sst_25 within-compartment CLR (bulk: production CelMod prediction;
    snRNA: count CLR), z-scored over the MATCHED donors (fig10_data.fit_cohort);
  * case beta from the Fig 10 primary model (z ~ case + age + sex + PMI);
  * Pearson r between the two modalities' z across the matched donors (per-donor agreement).
Bulk persons are unique across bulk cohorts (fig10_data dedup), so pairs never share a person.

Writes data/fig10_sst25_matched_donors.tsv and data/fig10_sst25_matched_effects.tsv.
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig10_sst25_matched.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fig10_data as F   # noqa: E402

TYPE = "Sst_25"
DISEASES = ["SCZ", "AD"]


def main():
    D = pd.read_csv(HERE / "data" / "fig10_donors.tsv", sep="\t", dtype={"person_id": str})
    D = D[D.disease.isin(DISEASES)]
    C = F.covariates()
    training = set(pd.read_csv(F.DECONV / "prod_labels.csv").donor.astype(str))
    mats = {}

    def mat(mod, name):
        if (mod, name) not in mats:
            mats[(mod, name)] = F.snrna_matrix(name) if mod == "snRNA" else F.bulk_matrix(name, training)[0]
        return mats[(mod, name)]

    don, eff = [], []
    for disease in DISEASES:
        d = D[D.disease == disease]
        bulk = d[d.modality == "bulk"].groupby("cohort").person_id.apply(set)
        sn = d[d.modality == "snRNA"].groupby("cohort").person_id.apply(set)
        for bc in bulk.index:
            for sc in sn.index:
                ids = sorted(bulk[bc] & sn[sc])
                if not ids:
                    continue
                m = C.reindex(ids).copy()
                case = d[(d.cohort == bc) & d.person_id.isin(ids)].set_index("person_id").case
                case_sn = d[(d.cohort == sc) & d.person_id.isin(ids)].set_index("person_id").case
                assert (case.reindex(ids) == case_sn.reindex(ids)).all(), (bc, sc, "dx differs")
                m["case"] = case.reindex(ids).values
                n_case, n_ctrl = int(m.case.sum()), int((1 - m.case).sum())
                if min(n_case, n_ctrl) < F.MIN_PER_GROUP:
                    print(f"{disease} {bc} x {sc}: {len(ids)} shared ({n_case}/{n_ctrl}), below floor, skipped")
                    continue
                z = {}
                for mod, name in (("bulk", bc), ("snRNA", sc)):
                    X = mat(mod, name)
                    rows, _, _ = F.fit_cohort(X, m, [TYPE], F.COVARS)
                    assert len(rows) == 1 and rows[0]["n"] == len(ids), (mod, name)
                    y = X.loc[ids, TYPE].astype(float).values
                    z[mod] = (y - y.mean()) / y.std(ddof=1)
                    eff.append(dict(disease=disease, bulk_cohort=bc, snrna_cohort=sc, modality=mod,
                                    n=len(ids), n_case=n_case, n_ctrl=n_ctrl,
                                    beta=rows[0]["beta"], se=rows[0]["se"], p=rows[0]["p"]))
                r = float(np.corrcoef(z["bulk"], z["snRNA"])[0, 1])
                for e in eff[-2:]:
                    e["r_bulk_snrna"] = r
                for mod in ("bulk", "snRNA"):
                    don.append(pd.DataFrame(dict(disease=disease, bulk_cohort=bc, snrna_cohort=sc,
                                                 modality=mod, person_id=ids, case=m.case.values,
                                                 z=z[mod])))
                b, s = eff[-2], eff[-1]
                print(f"{disease} {bc:20s} x {sc:15s} n={len(ids):3d} ({n_case}/{n_ctrl})  "
                      f"bulk beta {b['beta']:+.2f} p {b['p']:.3g} | snRNA beta {s['beta']:+.2f} "
                      f"p {s['p']:.3g} | r(bulk, snRNA) {r:+.2f}")
    pd.concat(don).to_csv(HERE / "data" / "fig10_sst25_matched_donors.tsv", sep="\t", index=False)
    pd.DataFrame(eff).to_csv(HERE / "data" / "fig10_sst25_matched_effects.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
