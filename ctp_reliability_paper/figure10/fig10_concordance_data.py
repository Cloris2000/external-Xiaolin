"""Concordance of bulk-deconvolved and snRNA-seq composition, Fig 10 supplement.

Two levels, both on the 86 Fig 4 Core Supertypes:

  MEASUREMENT (per donor). On the same-donor pairs of fig10_sst25_matched.py (people in both
  a bulk and an snRNA analysis set, >= 10 cases and controls), Pearson r between the bulk
  CelMod prediction and the snRNA within-compartment CLR across donors, per supertype.

  DISEASE EFFECT. For each pair, the case beta of the Fig 10 primary model (z ~ case + age +
  sex + PMI, outcome z-scored over the matched donors) per supertype and modality; Spearman
  rho between the 86 bulk and 86 snRNA betas, and the principal-axis (TLS) slope of bulk on
  snRNA. Because both betas come from the SAME people, they correlate even with no disease
  signal (the measurements correlate, see above). The p value is therefore from permuting
  the case labels (N_PERM times, one permutation shared by both modalities), which keeps the
  donors, covariates, between-type correlation and measurement agreement and removes only
  the disease. The meta-level (all cohorts, different donors) betas are in
  fig10_effects_meta.tsv and are summarised by the R script without a p value.

Writes data/fig10_concordance_donor_r.tsv, fig10_concordance_betas.tsv,
fig10_concordance_summary.tsv, fig10_concordance_sst25_donors.tsv, fig10_concordance_null.tsv.
TLS slope is carried but unstable at weak rho; the SD ratio of the betas is the scale summary.
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig10_concordance_data.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fig10_data as F   # noqa: E402

N_PERM = 2000
SEED = 20261117
TYPE = "Sst_25"


def tls_slope(x, y):
    """principal-axis slope of y on x (both have error)."""
    c = np.cov(x, y)
    w, v = np.linalg.eigh(c)
    a = v[:, np.argmax(w)]
    return float(a[1] / a[0])


def betas(M_cov, case, Yz):
    """case beta for every column of Yz (columns already z-scored), one shared design."""
    M = np.column_stack([np.ones(len(case)), case, M_cov])
    return np.linalg.lstsq(M, Yz, rcond=None)[0][1]


def zcols(Y):
    return (Y - Y.mean(0)) / Y.std(0, ddof=1)


def main():
    C = F.covariates()
    core = F.core_types()
    training = set(pd.read_csv(F.DECONV / "prod_labels.csv").donor.astype(str))
    E = pd.read_csv(HERE / "data" / "fig10_sst25_matched_effects.tsv", sep="\t")
    MD = pd.read_csv(HERE / "data" / "fig10_sst25_matched_donors.tsv", sep="\t", dtype={"person_id": str})
    pairs = E[["disease", "bulk_cohort", "snrna_cohort"]].drop_duplicates().itertuples(index=False)
    rng = np.random.default_rng(SEED)

    rrows, brows, srows, drows, nrows = [], [], [], [], []
    for disease, bc, sc in pairs:
        g = MD[(MD.disease == disease) & (MD.bulk_cohort == bc) & (MD.snrna_cohort == sc) &
               (MD.modality == "bulk")]
        ids = g.person_id.tolist()
        case = g.case.values.astype(float)
        Xb = F.bulk_matrix(bc, training)[0].loc[ids, core].astype(float).values
        Xs = F.snrna_matrix(sc).loc[ids, core].astype(float).values
        assert np.isfinite(Xb).all() and np.isfinite(Xs).all(), (bc, sc, "NaN in a Core type")
        cov = C.reindex(ids)[F.COVARS].values.astype(float)
        assert np.isfinite(cov).all()
        Zb, Zs = zcols(Xb), zcols(Xs)
        label = f"{bc} x {sc}"

        # measurement agreement per supertype
        r = np.array([np.corrcoef(Xb[:, j], Xs[:, j])[0, 1] for j in range(len(core))])
        for j, ct in enumerate(core):
            rrows.append(dict(disease=disease, pair=label, bulk_cohort=bc, snrna_cohort=sc,
                              cell_type=ct, n=len(ids), r=r[j]))

        # disease-effect concordance, observed and label-permuted
        bb, bs = betas(cov, case, Zb), betas(cov, case, Zs)
        j25 = core.index(TYPE)
        ref = E[(E.disease == disease) & (E.bulk_cohort == bc) & (E.snrna_cohort == sc)].set_index("modality")
        assert np.isclose(bb[j25], ref.loc["bulk", "beta"]) and np.isclose(bs[j25], ref.loc["snRNA", "beta"]), \
            "vectorised betas differ from fig10_sst25_matched.py"
        rho = spearmanr(bb, bs).correlation
        null = np.empty(N_PERM)
        for i in range(N_PERM):
            pc = rng.permutation(case)
            null[i] = spearmanr(betas(cov, pc, Zb), betas(cov, pc, Zs)).correlation
        p_perm = (1 + (null >= rho).sum()) / (N_PERM + 1)
        for j, ct in enumerate(core):
            brows.append(dict(disease=disease, pair=label, bulk_cohort=bc, snrna_cohort=sc,
                              cell_type=ct, beta_bulk=bb[j], beta_snrna=bs[j]))
        srows.append(dict(disease=disease, pair=label, bulk_cohort=bc, snrna_cohort=sc,
                          n=len(ids), n_case=int(case.sum()), n_ctrl=int((1 - case).sum()),
                          median_donor_r=float(np.median(r)), sst25_donor_r=float(r[j25]),
                          sst25_donor_r_rank=int((r > r[j25]).sum() + 1),
                          rho_beta=float(rho), null_rho_median=float(np.median(null)),
                          null_rho_q95=float(np.quantile(null, 0.95)), p_perm=float(p_perm),
                          sign_agree=float(np.mean(np.sign(bb) == np.sign(bs))),
                          tls_slope=tls_slope(bs, bb),
                          sd_ratio_bulk_over_snrna=float(bb.std(ddof=1) / bs.std(ddof=1)),
                          sd_beta_bulk=float(bb.std(ddof=1)), sd_beta_snrna=float(bs.std(ddof=1))))
        nrows.append(pd.DataFrame(dict(disease=disease, pair=label, null_rho=null)))
        drows.append(pd.DataFrame(dict(disease=disease, pair=label, person_id=ids, case=case,
                                       z_bulk=Zb[:, j25], z_snrna=Zs[:, j25])))
        s = srows[-1]
        print(f"{disease} {label:38s} n={len(ids):3d}  donor r median {s['median_donor_r']:+.2f} "
              f"(Sst_25 {s['sst25_donor_r']:+.2f}, rank {s['sst25_donor_r_rank']}/86)  "
              f"beta rho {rho:+.2f} (null median {s['null_rho_median']:+.2f}, 95% {s['null_rho_q95']:+.2f}, "
              f"p_perm {p_perm:.4f})  sign agree {s['sign_agree']:.0%}  SD(beta) bulk/sn {s['sd_ratio_bulk_over_snrna']:.2f}")

    out = HERE / "data"
    pd.DataFrame(rrows).to_csv(out / "fig10_concordance_donor_r.tsv", sep="\t", index=False)
    pd.DataFrame(brows).to_csv(out / "fig10_concordance_betas.tsv", sep="\t", index=False)
    pd.DataFrame(srows).to_csv(out / "fig10_concordance_summary.tsv", sep="\t", index=False)
    pd.concat(drows).to_csv(out / "fig10_concordance_sst25_donors.tsv", sep="\t", index=False)
    pd.concat(nrows).to_csv(out / "fig10_concordance_null.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
