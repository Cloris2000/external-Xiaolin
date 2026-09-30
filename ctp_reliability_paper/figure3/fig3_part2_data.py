"""Figure 3 part 2 inputs: why sampled depth cannot simply be covaried for AD.

Design follows the PI's 68_ad_depth_explains.py / 67_tilt_is_it_real.py, on the
DFC 2026 caches:
  CLR within the neuronal compartment; per supertype s
    beta_s    phenotype effect, OLS  CLR ~ phenotype + age + sex
    gamma_s   depth sensitivity, slope of CLR on the donor's sampled depth in
              donors WITHOUT dementia (fallback: lowest phenotype quartile)
    delta     slope of sampled depth on the phenotype, all donors
    predicted gamma_s * delta, the effect a pure sampling shift would make
    beta_adj  phenotype effect with sampled depth added as a covariate
  Sampled depth is leave-one-class-out, as in 68: an inhibitory supertype is
  modelled with depth from excitatory neurons and vice versa, so a supertype
  never predicts itself.
  Gradient ("tilt") = slope of beta_s on the supertype's cortical depth across
  supertypes; share = predicted gradient / observed gradient, with a bootstrap
  over supertypes.

Depth prior and taxonomy as in fig3_data.py. Writes data/.
"""
import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr

import fig3_data as F

DATASETS = {"Green_DFC_person": "Green", "Mathys_DFC_person": "Mathys",
            "Multiome2025_DLPFC_DFC_person": "Multiome 2025",
            "PsychAD_full_DFC_person": "PsychAD"}
PHENOS = ["cerad_ad", "cogdx", "dementia"]
N_BOOT = 2000
SEED = 20260929


def ols(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    cov = float(r @ r) / (len(y) - X.shape[1]) * np.linalg.pinv(X.T @ X)
    return float(b[1]), float(np.sqrt(max(cov[1, 1], 0)))


def main():
    rng = np.random.default_rng(SEED)
    tax = pd.read_csv(F.TAX, sep="\t").set_index("supertype_label")
    cls = tax.class_id
    pri = pd.read_csv(F.RESULTS / "12_supertype_depth_wm.csv").set_index("supertype")
    mu = pri.depth_median[pri.index.isin(tax.index)]
    ph = pd.read_csv(F.RESULTS / "65_person_phenotypes.csv")
    ph["person_id"] = ph.person_id.astype(str)
    ph = ph.set_index("person_id")

    eff, summ = [], []
    for cache, label in DATASETS.items():
        cells = pd.read_parquet(F.CACHE / f"{cache}_cells.parquet", columns=["donor", "supertype"])
        t = pd.crosstab(cells.donor.astype(str), cells.supertype)
        neu = [c for c in t.columns if cls[c] in ("Excitatory", "Inhibitory")]
        t = t.loc[t[neu].sum(axis=1) >= F.MIN_CELLS, neu]
        exc = [c for c in neu if cls[c] == "Excitatory" and c in mu.index]
        inh = [c for c in neu if cls[c] == "Inhibitory" and c in mu.index]
        for pheno in PHENOS:
            M = ph.reindex(t.index)
            ok = (M[pheno].notna() & M.sex.notna()).values
            if ok.sum() < F.MIN_DONORS:
                print(f"  {label:13s} {pheno:9s} skipped: {ok.sum()} donors"); continue
            T, M = t[ok], M[ok].copy()
            M["age"] = M.age.fillna(M.age.median())
            Y = F.clr(T.values.astype(float))
            dep_exc, dep_inh = F.class_depth(T, exc, mu), F.class_depth(T, inh, mu)
            healthy = (M.dementia == 0).values
            if healthy.sum() < 25:
                healthy = (M[pheno] <= M[pheno].quantile(0.25)).values
            if healthy.sum() < 25:
                print(f"  {label:13s} {pheno:9s} skipped: {healthy.sum()} unaffected"); continue
            ph_v, age, sex = (M[c].values.astype(float) for c in (pheno, "age", "sex"))
            rows = []
            for j, ct in enumerate(T.columns):
                if ct not in mu.index:
                    continue
                dep = dep_exc if cls[ct] == "Inhibitory" else dep_inh
                fin = np.isfinite(dep)
                hh = fin & healthy
                if hh.sum() < 25 or fin.sum() < F.MIN_DONORS:
                    continue
                y = Y[:, j]
                beta, se = ols(y, np.column_stack([np.ones(len(y)), ph_v, age, sex]))
                Xa = np.column_stack([np.ones(fin.sum()), ph_v[fin], age[fin], sex[fin], dep[fin]])
                beta_adj, _ = ols(y[fin], Xa)
                g = linregress(dep[hh], y[hh]).slope
                d = linregress(ph_v[fin], dep[fin]).slope
                rows.append(dict(dataset=label, phenotype=pheno, cell_type=ct, cls=cls[ct],
                                 depth=mu[ct], beta=beta, se=se, beta_adj=beta_adj,
                                 gamma=g, delta=d, predicted=g * d))
            O = pd.DataFrame(rows)
            eff.append(O)
            lo = linregress(O.depth, O.beta).slope
            lp = linregress(O.depth, O.predicted).slope
            la = linregress(O.depth, O.beta_adj).slope
            boot = []
            for _ in range(N_BOOT):
                i = rng.integers(0, len(O), len(O))
                a = linregress(O.depth.values[i], O.beta.values[i]).slope
                if abs(a) > 1e-6:
                    boot.append(linregress(O.depth.values[i], O.predicted.values[i]).slope / a)
            r_obs, p_obs = spearmanr(O.depth, O.beta)
            r_adj, p_adj = spearmanr(O.depth, O.beta_adj)
            s = dict(dataset=label, phenotype=pheno, n_donors=int(ok.sum()),
                     n_healthy=int(healthy.sum()), n_types=len(O),
                     slope_obs=lo, slope_pred=lp, slope_adj=la,
                     share=lp / lo if abs(lo) > 1e-6 else np.nan,
                     share_lo=np.percentile(boot, 2.5), share_hi=np.percentile(boot, 97.5),
                     rho_obs=r_obs, p_obs=p_obs, rho_adj=r_adj, p_adj=p_adj,
                     rho_gamma_beta=spearmanr(O.gamma, O.beta)[0])
            summ.append(s)
            print(f"  {label:13s} {pheno:9s} n={s['n_donors']:4d} healthy={s['n_healthy']:4d} "
                  f"tilt rho {r_obs:+.2f} (p {p_obs:.1e}) -> adjusted {r_adj:+.2f}; "
                  f"share {100 * s['share']:6.0f}% [{100 * s['share_lo']:.0f}, {100 * s['share_hi']:.0f}]")
    pd.concat(eff).to_csv(F.OUT / "fig3_part2_effects.tsv", sep="\t", index=False)
    pd.DataFrame(summ).to_csv(F.OUT / "fig3_part2_summary.tsv", sep="\t", index=False)
    print(f"wrote fig3_part2_effects.tsv, fig3_part2_summary.tsv to {F.OUT}")


if __name__ == "__main__":
    main()
