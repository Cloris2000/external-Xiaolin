"""Sst_25 supplement to Figure 10: donor-level Sst_25 proportion (% of all neuronal nuclei)
in each snRNA-seq cohort, on exactly the donors that entered the Fig 10 primary model
(data/fig10_donors.tsv: case/control, age window, one person per modality x disease).

Neuronal = taxonomy class_label starting 'Neuronal' (Glutamatergic + GABAergic), the
compartment split of common/reliability_stats.py. Bulk has no proportion (CelMod predicts
a CLR-scale score), so it appears only in the forest plot (fig10_effects_by_cohort.tsv).

Bulk: Sst_25 production CelMod prediction (within-compartment CLR scale), z-scored within
cohort on the Fig 10 bulk analysis donors.

Writes data/fig10_sst25_props.tsv (snRNA) and data/fig10_sst25_bulk.tsv (bulk).
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig10_sst25_data.py
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import reliability_stats as RS   # noqa: E402
sys.path.insert(0, str(HERE))
import fig10_data as F           # noqa: E402  (bulk_matrix: name mapping + training-donor removal)

TYPE = "Sst_25"
DISEASES = ["SCZ", "AD"]


def main():
    D = pd.read_csv(HERE / "data" / "fig10_donors.tsv", sep="\t", dtype={"person_id": str})
    D = D[(D.modality == "snRNA") & D.disease.isin(DISEASES)]
    tax = RS.taxonomy().set_index("supertype_label")
    rows = []
    for cohort, g in D.groupby("cohort", sort=False):
        c = RS.load(cohort)
        neu = c.columns[tax.loc[c.columns, "class_label"].str.startswith(RS.NEURONAL_PREFIX)]
        assert TYPE in c.columns, cohort
        missing = set(g.person_id) - set(c.index)
        assert not missing, (cohort, len(missing))
        n_neu = c.loc[g.person_id, neu].sum(axis=1)
        n_sst = c.loc[g.person_id, TYPE]
        # within-compartment CLR (fig10_data.snrna_matrix), z-scored over these donors per
        # disease = the snRNA outcome of the Fig 10 regression
        clr = F.snrna_matrix(cohort).loc[g.person_id, TYPE].astype(float)
        assert clr.notna().all(), cohort
        out = g.assign(n_neuronal=n_neu.values, n_sst25=n_sst.values,
                       pct_of_neurons=100 * n_sst.values / n_neu.values, clr=clr.values)
        out["z"] = out.groupby("disease").clr.transform(lambda v: (v - v.mean()) / v.std(ddof=1))
        rows.append(out)
        for d, h in out.groupby("disease"):
            print(f"{d:3s} {cohort:20s} n={len(h):4d}  median neurons/donor {h.n_neuronal.median():7.0f}  "
                  f"median Sst_25 % ctrl {h[h.case == 0].pct_of_neurons.median():.3f} "
                  f"case {h[h.case == 1].pct_of_neurons.median():.3f}  "
                  f"zero-neuron donors {(h.n_neuronal == 0).sum()}")
    P = pd.concat(rows)
    P.to_csv(HERE / "data" / "fig10_sst25_props.tsv", sep="\t", index=False)
    print(f"wrote {len(P)} snRNA donor rows")

    # bulk: production CelMod prediction of the within-compartment CLR (prod__consensus,
    # trained on the z-scored Green/Mathys/Multiome consensus CLR), z-scored within cohort on
    # the analysis donors -- exactly the outcome of the Fig 10 regression (fit_cohort yz)
    B = pd.read_csv(HERE / "data" / "fig10_donors.tsv", sep="\t", dtype={"person_id": str})
    B = B[(B.modality == "bulk") & B.disease.isin(DISEASES)]
    training = set(pd.read_csv(F.DECONV / "prod_labels.csv").donor.astype(str))
    rows = []
    for (disease, cohort), g in B.groupby(["disease", "cohort"], sort=False):
        X, _ = F.bulk_matrix(cohort, training)
        y = X.loc[g.person_id, TYPE].astype(float)
        assert y.notna().all(), cohort
        out = g.assign(pred_clr=y.values, z=((y - y.mean()) / y.std(ddof=1)).values)
        rows.append(out)
        print(f"{disease:3s} bulk {cohort:24s} n={len(out):4d}  mean z ctrl {out[out.case == 0].z.mean():+.3f} "
              f"case {out[out.case == 1].z.mean():+.3f}")
    Q = pd.concat(rows)
    Q.to_csv(HERE / "data" / "fig10_sst25_bulk.tsv", sep="\t", index=False)
    print(f"wrote {len(Q)} bulk donor rows")


if __name__ == "__main__":
    main()
