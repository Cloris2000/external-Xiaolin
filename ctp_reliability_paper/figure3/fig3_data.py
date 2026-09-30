"""Figure 3 inputs: sampled cortical depth, DFC 2026 taxonomy (DRAFT).

Question: does inconsistent sampling depth between tissue blocks contribute to
poor agreement between datasets, and can depth simply be covaried for AD?

Depth prior: results/12_supertype_depth_wm.csv (SCZ Xenium, DLPFC; the prior
used by 25/68). It is keyed by supertype NAME and was built before DFC 2026, so
it is joined by name; DFC supertypes absent from it (15 reactive states and
Pvalb_14) are left out of every depth score, and the share of nuclei this drops
is reported per dataset.

Depth scores per donor (first moment of the sampled depth distribution):
  depth_all  composition-weighted mean depth over all labelled nuclei (as 25)
  depth_exc  excitatory neurons only; depth_inh  inhibitory neurons only (as 68)
  wm_frac    composition-weighted P(white matter), as 25
A class-restricted score needs >= 50 nuclei of that class (68's rule).

Outputs (data/):
  fig3_donor_depth.tsv       one row per (dataset, donor)
  fig3_dataset_profiles.tsv  mean depth density per dataset on a 0-1 grid
  fig3_pair_donor.tsv        matched donors in the 6 independent pairs: depth of
                             the same person in both datasets, and agreement of
                             that person's INHIBITORY composition between them
                             (depth taken from EXCITATORY cells, so the two
                             measures share no nuclei)
  fig3_ad_effects.tsv        per (dataset, neuronal supertype): CERAD effect on
                             CLR and depth sensitivity gamma (68's design)
  fig3_ad_donor.tsv          donor depth and CERAD for the AD datasets
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress

CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
RESULTS = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv")
OUT = Path(__file__).resolve().parent / "data"

DATASETS = [
    "SEAAD_DFC_2026",
    "Batiuk_full_DFC_person", "BrainSCOPE_CMC_DFC_person",
    "BrainSCOPE_Multiome_DFC_person", "BrainSCOPE_UCLA_DFC_person",
    "Cain_2023_DFC_person", "Frohlich_full_DFC_person", "Green_DFC_person",
    "Lau_2020_DFC_person", "Leng_2021_DFC_person", "Mathys_2019_DFC_person",
    "Mathys_DFC_person", "Morabito_2021_DFC_person",
    "Multiome2025_DLPFC_DFC_person", "Olah_2020_DFC_person",
    "PsychAD_full_DFC_person", "Ruzicka_full_DFC_person", "Yang_2022_DFC_person",
    "Zhou_2020_DFC_person",
]
# the six independent pairs of 114/117
PAIRS = [
    ("Green_DFC_person", "Mathys_DFC_person", "Green x Mathys"),
    ("Green_DFC_person", "Multiome2025_DLPFC_DFC_person", "Green x Multiome"),
    ("Mathys_DFC_person", "Multiome2025_DLPFC_DFC_person", "Mathys x Multiome"),
    ("Multiome2025_DLPFC_DFC_person", "PsychAD_full_DFC_person", "Multiome x PsychAD"),
    ("BrainSCOPE_CMC_DFC_person", "PsychAD_full_DFC_person", "BrainSCOPE CMC x PsychAD"),
    ("PsychAD_full_DFC_person", "Ruzicka_full_DFC_person", "PsychAD x Ruzicka"),
]
AD_DATASETS = ["Green_DFC_person", "Mathys_DFC_person"]
PHENO = "cerad_ad"

GRID = np.linspace(0.0, 1.0, 201)
MIN_CELLS = 100            # donor floor for a profile (25)
MIN_CLASS_CELLS = 50       # class-restricted depth score (68)
MIN_DONORS = 40
PSEUDO = 0.5


def short(name):
    return name.replace("_DFC_person", "").replace("_full", "")


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def ols_beta(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    dof = len(y) - X.shape[1]
    cov = float(r @ r) / dof * np.linalg.pinv(X.T @ X)
    return float(b[1]), float(np.sqrt(max(cov[1, 1], 0)))


def class_depth(tab, cols, mu):
    C = tab[cols].values.astype(float)
    tot = C.sum(axis=1)
    out = np.full(len(C), np.nan)
    ok = tot >= MIN_CLASS_CELLS
    out[ok] = (C[ok] / tot[ok, None]) @ mu[cols].values
    return out


def main():
    OUT.mkdir(exist_ok=True)
    tax = pd.read_csv(TAX, sep="\t").set_index("supertype_label")
    cls = tax.class_id
    pri = pd.read_csv(RESULTS / "12_supertype_depth_wm.csv").set_index("supertype")
    pri = pri[pri.index.isin(tax.index)]
    mu, pwm = pri.depth_median, pri.p_wm.fillna(0)
    sig = pri.depth_sd.clip(lower=0.02)
    print(f"depth prior covers {len(pri)} of {len(tax)} DFC supertypes")

    tabs, donor_rows, prof_rows = {}, [], []
    for name in DATASETS:
        f = CACHE / f"{name}_cells.parquet"
        if not f.exists():
            print(f"  {name}: no cache, skipped"); continue
        cells = pd.read_parquet(f, columns=["donor", "supertype"])
        cells["donor"] = cells.donor.astype(str)
        unk = ~cells.supertype.isin(tax.index)
        if unk.any():
            raise SystemExit(f"{name}: {unk.sum()} nuclei with non-DFC labels, e.g. "
                             f"{cells.supertype[unk].unique()[:5]}")
        t = pd.crosstab(cells.donor, cells.supertype)
        tabs[name] = t
        have = [c for c in t.columns if c in mu.index]
        dropped = 1 - t[have].values.sum() / t.values.sum()
        T = t[have][t[have].sum(axis=1) >= MIN_CELLS]
        P = T.div(T.sum(axis=1), axis=0)
        exc = [c for c in have if cls[c] == "Excitatory"]
        inh = [c for c in have if cls[c] == "Inhibitory"]
        d = pd.DataFrame({
            "dataset": short(name), "donor": T.index, "n_nuclei": T.sum(axis=1).values,
            "depth_all": P.values @ mu[have].values,
            "depth_exc": class_depth(T, exc, mu), "depth_inh": class_depth(T, inh, mu),
            "wm_frac": P.values @ pwm[have].values,
        })
        donor_rows.append(d)
        B = np.exp(-((GRID[None, :] - mu[have].values[:, None]) ** 2)
                   / (2 * sig[have].values[:, None] ** 2))
        B /= sig[have].values[:, None] * np.sqrt(2 * np.pi)
        prof_rows.append(pd.DataFrame({"dataset": short(name), "depth": GRID,
                                       "density": (P.values @ B).mean(0)}))
        print(f"  {short(name):22s} {len(T):5d} donors  depth {d.depth_all.median():.3f}  "
              f"WM {100 * d.wm_frac.median():4.1f}%  nuclei without prior {100 * dropped:.1f}%")
    D = pd.concat(donor_rows)
    D.to_csv(OUT / "fig3_donor_depth.tsv", sep="\t", index=False)
    pd.concat(prof_rows).to_csv(OUT / "fig3_dataset_profiles.tsv", sep="\t", index=False)

    # ---- matched donors: depth mismatch vs composition agreement -------------
    Dk = D.set_index(["dataset", "donor"])
    pair_rows = []
    for a, b, label in PAIRS:
        ta, tb = tabs[a], tabs[b]
        shared = sorted(set(ta.index) & set(tb.index))
        inh = [c for c in ta.columns if c in tb.columns and cls[c] == "Inhibitory"]
        A, Bm = ta.loc[shared, inh].values.astype(float), tb.loc[shared, inh].values.astype(float)
        ok = (A.sum(axis=1) >= MIN_CLASS_CELLS) & (Bm.sum(axis=1) >= MIN_CLASS_CELLS)
        # keep inhibitory types with a mean share >= 0.5% in both, so the
        # per-donor correlation is not driven by near-empty columns
        keep = ((A[ok] / A[ok].sum(1, keepdims=True)).mean(0) >= 0.005) & \
               ((Bm[ok] / Bm[ok].sum(1, keepdims=True)).mean(0) >= 0.005)
        za, zb = clr(A[ok][:, keep]), clr(Bm[ok][:, keep])
        za = (za - za.mean(0)) / za.std(0, ddof=1)       # remove dataset offsets
        zb = (zb - zb.mean(0)) / zb.std(0, ddof=1)
        agree = np.array([np.corrcoef(za[i], zb[i])[0, 1] for i in range(len(za))])
        don = np.array(shared)[ok]
        g = lambda ds, col: Dk[col].reindex(list(zip([short(ds)] * len(don), don))).values
        pair_rows.append(pd.DataFrame({
            "pair": label, "donor": don, "n_inh_types": int(keep.sum()),
            "agree_inh": agree,
            "depth_exc_a": g(a, "depth_exc"), "depth_exc_b": g(b, "depth_exc"),
            "depth_all_a": g(a, "depth_all"), "depth_all_b": g(b, "depth_all"),
            "wm_a": g(a, "wm_frac"), "wm_b": g(b, "wm_frac")}))
        print(f"  {label:26s} {ok.sum():4d} donors, {keep.sum()} inhibitory types")
    pd.concat(pair_rows).to_csv(OUT / "fig3_pair_donor.tsv", sep="\t", index=False)

    # ---- AD: CERAD effect vs depth sensitivity (68's design, DFC) ------------
    ph = pd.read_csv(RESULTS / "65_person_phenotypes.csv")
    ph["person_id"] = ph.person_id.astype(str)
    ph = ph.set_index("person_id")
    eff_rows, addon_rows = [], []
    for name in AD_DATASETS:
        t = tabs[name]
        neu = [c for c in t.columns if cls[c] in ("Excitatory", "Inhibitory")]
        M = ph.reindex(t.index)
        ok = M[PHENO].notna() & M.sex.notna() & (t[neu].sum(axis=1) >= MIN_CELLS)
        T, M = t.loc[ok, neu], M[ok].copy()
        M["age"] = M.age.fillna(M.age.median())
        Y = clr(T.values.astype(float))
        exc = [c for c in neu if cls[c] == "Excitatory" and c in mu.index]
        inh = [c for c in neu if cls[c] == "Inhibitory" and c in mu.index]
        dep_exc, dep_inh = class_depth(T, exc, mu), class_depth(T, inh, mu)
        healthy = (M.dementia == 0).values
        X = np.column_stack([np.ones(len(M)), M[[PHENO, "age", "sex"]].values.astype(float)])
        for j, ct in enumerate(neu):
            if ct not in mu.index:
                continue
            dep = dep_exc if cls[ct] == "Inhibitory" else dep_inh   # never predicts itself
            fin = np.isfinite(dep)
            hh = fin & healthy
            if hh.sum() < 25:
                continue
            beta, se = ols_beta(Y[:, j], X)
            gm = linregress(dep[hh], Y[hh, j])
            dl = linregress(M[PHENO].values[fin], dep[fin])
            eff_rows.append(dict(dataset=short(name), cell_type=ct, cls=cls[ct],
                                 depth=mu[ct], beta_cerad=beta, se_cerad=se,
                                 gamma=gm.slope, delta=dl.slope,
                                 predicted=gm.slope * dl.slope))
        addon_rows.append(pd.DataFrame({
            "dataset": short(name), "donor": T.index, PHENO: M[PHENO].values,
            "dementia": M.dementia.values, "age": M.age.values, "sex": M.sex.values,
            "depth_exc": dep_exc, "depth_inh": dep_inh}))
        print(f"  AD {short(name):10s} {len(M)} donors with CERAD, "
              f"{int(healthy.sum())} without dementia")
    pd.DataFrame(eff_rows).to_csv(OUT / "fig3_ad_effects.tsv", sep="\t", index=False)
    pd.concat(addon_rows).to_csv(OUT / "fig3_ad_donor.tsv", sep="\t", index=False)
    print(f"wrote 5 tables to {OUT}")


if __name__ == "__main__":
    main()
