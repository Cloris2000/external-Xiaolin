"""Figure 10 inputs: disease effects on the Core Supertypes, bulk-deconvolved vs snRNA-seq.

For every Fig 4 Core Supertype (ci_lo of the median symmetric ICC > 0; 86 types) and each
disease (SCZ, BD, AD), the case-vs-control effect is estimated twice:
  bulk    production CelMod predictions (results_2026-10-21c deconv, prod__consensus)
  snRNA   CLR of supertype nucleus counts, within compartment (neuronal / non-neuronal),
          pseudocount 0.5 -- the CLR of common/reliability_stats.py, imported, not copied.

Design (user decisions 2026-11-12):
  * all cohorts with the disease in each modality; per-cohort OLS, then DerSimonian-Laird
    random-effects meta-analysis per modality (wp3 121/197/200 convention).
  * outcome z-scored within cohort (on the analysis donors) so beta = SD units, comparable
    between modalities.  y_z ~ case + age + sex + PMI.  PRIMARY MODEL HAS NO RIN: the paper's
    standing rule is "do not use RIN as a covariate" (2026-09-28); a +RIN model is carried
    as a sensitivity column only.
  * SCZ and BD: adults 18 <= age < 70 (Shreejoy's SCZ convention, 197/200); AD: all ages.
  * cases = registry harmonised dx == disease; controls = dx == 'Control'; every other dx
    (MCI, OtherDementia, PathAging, Schizoaffective, ...) is left out of that contrast.
  * each person counted once per modality x disease: cohorts are listed in priority order
    and a person already used by an earlier cohort is dropped from later ones (e.g. Green x
    Mathys share 227 ROSMAP donors).  The two modalities may share people.
  * bulk donors that trained the production CelMod (deconv/prod_labels.csv, 886 ROSMAP
    people) are excluded, so every bulk prediction is out of sample; diagnosis never enters
    CelMod.
  * >= MIN_PER_GROUP cases and controls per cohort after all filters.
  * NeuN-sorted datasets (Batiuk, Bast_2025) are not used: sorting fixes the neuron fraction,
    so their composition is not comparable.  Ruzicka is not used (duplicate of BrainSCOPE CMC
    / SZBD, analysis_set = no).

Inputs (read-only):
  /project/rrg-shreejoy/brain-data-registry/releases/2026-10-21c/covariates.tsv   dx/age/sex/pmi/rin
  /scratch/shreejoy/composition/results_2026-10-21c/counts/<ds>__supertype.parquet
  /scratch/shreejoy/composition/results_2026-10-21c/deconv/<cohort>_DFC_supertype_prod__consensus.parquet
  ../figure4/data/fig4a_icc_vs_deconv.tsv                                        Core set
  celltype-composition/refs/taxonomy_DFC_2026.tsv                                class labels

Writes data/fig10_*.tsv.
Run: PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig10_data.py
"""
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import reliability_stats as RS   # noqa: E402  (load, count_grid, clr, taxonomy)

warnings.filterwarnings("ignore", category=FutureWarning)

COV = Path("/project/rrg-shreejoy/brain-data-registry/releases/2026-10-21c/covariates.tsv")
DECONV = Path("/scratch/shreejoy/composition/results_2026-10-21c/deconv")
CORE_TSV = HERE.parent / "figure4" / "data" / "fig4a_icc_vs_deconv.tsv"
OUT = HERE / "data"

DISEASES = ["SCZ", "BD", "AD"]
AGE_RANGE = {"SCZ": (18, 70), "BD": (18, 70), "AD": (None, None)}   # [lo, hi)
MIN_PER_GROUP = 10
COVARS = ["age", "sex", "pmi"]
MODELS = {"primary": COVARS, "plus_rin": COVARS + ["rin"]}

# priority order within each list = which cohort keeps a person seen in several
SNRNA = {
    "SCZ": ["PsychAD_MSBB", "PsychAD_HBCC", "BrainSCOPE_CMC", "BrainSCOPE_SZBD", "Frohlich"],
    "BD":  ["PsychAD_HBCC", "BrainSCOPE_SZBD", "PsychAD_MSBB", "BrainSCOPE_Multiome"],
    "AD":  ["PsychAD_MSBB", "Green", "Mathys", "Multiome2025_DLPFC", "CUIMC2",
            "PsychAD_ROSMAP", "SEAAD"],
}
BULK = {   # frontal-cortex bulk only; MSBB uses BA10 (FP), one region per person
    "SCZ": ["CommonMind_MSSM", "GVEX", "NIMH_HBCC", "LIBD_szControl_LIBD",
            "LIBD_szControl_HBCC", "CommonMind_Pitt", "CommonMind_Penn"],
    "BD":  ["GVEX", "NIMH_HBCC", "CommonMind_Pitt", "CommonMind_MSSM"],
    "AD":  ["ROSMAP", "MSBB_FP", "Diverse_Rush2_DLPFC", "Diverse_MayoEmory_DLPFC",
            "Diverse_MSSM_DLPFC", "Diverse_Columbia_DLPFC"],
}


def covariates():
    c = pd.read_csv(COV, sep="\t", low_memory=False)
    assert c.person_id.is_unique
    c = c.set_index("person_id")
    out = pd.DataFrame(index=c.index)
    out["dx"] = c.dx.astype("string")
    out["age"] = pd.to_numeric(c.age, errors="coerce")
    out["sex"] = c.sex.map({"M": 1.0, "F": 0.0})
    out["pmi"] = pd.to_numeric(c.pmi, errors="coerce")
    out["rin"] = pd.to_numeric(c.rin, errors="coerce")
    return out


def core_types():
    a = pd.read_csv(CORE_TSV, sep="\t")
    core = a[a.concordant].cell_type.tolist()
    assert len(core) == 86, len(core)
    return core


def snrna_matrix(ds):
    """donor x supertype within-compartment CLR, as reliability_stats (all donors of ds)."""
    counts = RS.load(ds)
    blocks = []
    for grp, t in RS.count_grid(counts, "supertype").items():
        keep = t.sum(axis=1) > 0
        blocks.append(pd.DataFrame(RS.clr(t[keep].values.astype(float)),
                                   index=t.index[keep], columns=t.columns))
    return pd.concat(blocks, axis=1)   # a donor lacking one compartment gets NaN there


def bulk_matrix(cohort, training):
    """Registry deconv columns use underscore names ('L2_3_IT_13'); map them back to the
    taxonomy labels with the normalisation of figure4/fig4_data.py add_extra_legs."""
    X = pd.read_parquet(DECONV / f"{cohort}_DFC_supertype_prod__consensus.parquet")
    X.index = X.index.astype(str)
    norm = lambda x: re.sub(r"[^0-9A-Za-z-]+", "_", x).strip("_")
    labels = RS.taxonomy().supertype_label
    lab = {norm(x): x for x in labels.unique()}
    assert len(lab) == labels.nunique(), "name normalisation collides"
    miss = set(X.columns) - set(lab)
    assert not miss, (cohort, sorted(miss)[:5])
    X = X.rename(columns=lab)
    n0 = len(X)
    X = X[~X.index.isin(training)]
    return X, n0 - len(X)


def ols_case(y, M):
    """beta, se, p of column 1 (case indicator)."""
    b, *_ = np.linalg.lstsq(M, y, rcond=None)
    r = y - M @ b
    dof = len(y) - M.shape[1]
    cov = float(r @ r) / dof * np.linalg.pinv(M.T @ M)
    se = float(np.sqrt(max(cov[1, 1], 0)))
    return float(b[1]), se, float(2 * tdist.sf(abs(b[1] / se), dof))


def dl_meta(b, se):
    """DerSimonian-Laird random effects (as wp3 200_scz_by_bank.py)."""
    b, se = np.asarray(b, float), np.asarray(se, float)
    w = 1 / se ** 2
    bf = (b * w).sum() / w.sum()
    if len(b) == 1:
        return float(b[0]), float(se[0]), 0.0, np.nan
    Q = float((w * (b - bf) ** 2).sum()); df = len(b) - 1
    Cq = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / Cq) if Cq > 0 else 0.0
    wr = 1 / (se ** 2 + tau2)
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    return float((b * wr).sum() / wr.sum()), float(np.sqrt(1 / wr.sum())), tau2, I2


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def donors_for(X, C, disease, used):
    """analysis donors of one cohort: case/control, age window, not used by an earlier cohort."""
    m = C.reindex(X.index)
    m = m[m.dx.isin([disease, "Control"])].copy()
    lo, hi = AGE_RANGE[disease]
    if lo is not None:
        m = m[(m.age >= lo) & (m.age < hi)]
    m = m[~m.index.isin(used)]
    m["case"] = (m.dx == disease).astype(float)
    return m


def fit_cohort(X, m, types, cols):
    mm = m.dropna(subset=["case"] + cols)
    n_case, n_ctrl = int(mm.case.sum()), int((1 - mm.case).sum())
    rows = []
    if min(n_case, n_ctrl) < MIN_PER_GROUP:
        return rows, n_case, n_ctrl
    for ct in types:
        if ct not in X.columns:
            continue
        y = X.loc[mm.index, ct].astype(float)
        ok = y.notna().values
        yy = y.values[ok]
        if ok.sum() < 2 * MIN_PER_GROUP or np.std(yy) < 1e-9:
            continue
        yz = (yy - yy.mean()) / yy.std(ddof=1)
        M = np.column_stack([np.ones(ok.sum()), mm.case.values[ok]] +
                            [mm[c].values[ok] for c in cols])
        b, se, p = ols_case(yz, M)
        rows.append(dict(cell_type=ct, n=int(ok.sum()),
                         n_case=int(mm.case.values[ok].sum()),
                         n_ctrl=int((1 - mm.case.values[ok]).sum()),
                         beta=b, se=se, p=p))
    return rows, n_case, n_ctrl


def main():
    OUT.mkdir(exist_ok=True)
    C = covariates()
    core = core_types()
    tax = RS.taxonomy().set_index("supertype_label")
    training = set(pd.read_csv(DECONV / "prod_labels.csv").donor.astype(str))

    mats = {}
    def get(mod, name):
        if (mod, name) not in mats:
            if mod == "snRNA":
                mats[(mod, name)] = (snrna_matrix(name), 0)
            else:
                mats[(mod, name)] = bulk_matrix(name, training)
            X = mats[(mod, name)][0]
            miss = [c for c in core if c not in X.columns]
            assert len(miss) < len(core), (mod, name, "no Core supertype column matches")
            if miss:
                print(f"  note: {mod} {name} lacks {len(miss)} Core types: {miss[:5]}")
        return mats[(mod, name)]

    per, cohorts, donors = [], [], []
    for disease in DISEASES:
        for mod, plan in (("bulk", BULK), ("snRNA", SNRNA)):
            for model, cols in MODELS.items():
                used = set()
                for name in plan[disease]:
                    X, n_train = get(mod, name)
                    m = donors_for(X, C, disease, used)
                    rows, n_case, n_ctrl = fit_cohort(X, m, core, cols)
                    entered = bool(rows)
                    if model == "primary":
                        cohorts.append(dict(disease=disease, modality=mod, cohort=name,
                                            n_case=n_case, n_ctrl=n_ctrl, entered=entered,
                                            n_training_excluded=n_train,
                                            age_mean=float(m.age.mean()) if len(m) else np.nan))
                        print(f"{disease:3s} {mod:5s} {name:24s} cases {n_case:4d} ctrl {n_ctrl:4d}"
                              f"{'' if entered else '   -> below floor, not used'}"
                              f"{f'   ({n_train} CelMod training donors removed)' if n_train else ''}")
                    if entered:
                        mm = m.dropna(subset=["case"] + cols)
                        used |= set(mm.index)
                        if model == "primary":
                            donors.append(pd.DataFrame({"disease": disease, "modality": mod,
                                                        "cohort": name, "person_id": mm.index,
                                                        "case": mm.case.values,
                                                        "age": mm.age.values}))
                    for r in rows:
                        r.update(disease=disease, modality=mod, model=model, cohort=name)
                    per.extend(rows)
    P = pd.DataFrame(per)

    meta = []
    for (disease, mod, model, ct), g in P.groupby(["disease", "modality", "model", "cell_type"]):
        b, se, tau2, I2 = dl_meta(g.beta, g.se)
        meta.append(dict(disease=disease, modality=mod, model=model, cell_type=ct,
                         beta=b, se=se, p=float(2 * norm.sf(abs(b / se))), tau2=tau2, I2=I2,
                         n_cohorts=len(g), n_donors=int(g.n.sum()),
                         n_case=int(g.n_case.sum()), n_ctrl=int(g.n_ctrl.sum()),
                         n_same_sign=int((np.sign(g.beta) == np.sign(b)).sum()),
                         cohorts=",".join(g.cohort)))
    M = pd.DataFrame(meta)
    M["fdr"] = M.groupby(["disease", "modality", "model"]).p.transform(lambda p: bh(p.values))
    M["class_label"] = M.cell_type.map(tax.class_label)
    M["subclass_label"] = M.cell_type.map(tax.subclass_label)
    assert M.class_label.notna().all()

    P.to_csv(OUT / "fig10_effects_by_cohort.tsv", sep="\t", index=False)
    M.to_csv(OUT / "fig10_effects_meta.tsv", sep="\t", index=False)
    pd.DataFrame(cohorts).to_csv(OUT / "fig10_cohorts.tsv", sep="\t", index=False)
    pd.concat(donors).to_csv(OUT / "fig10_donors.tsv", sep="\t", index=False)   # primary-model donors

    pm = M[M.model == "primary"]
    print("\nprimary model: Core types tested / FDR < 0.05, and bulk-vs-snRNA agreement")
    for disease in DISEASES:
        w = pm[pm.disease == disease].pivot(index="cell_type", columns="modality", values="beta")
        line = []
        for mod in ("bulk", "snRNA"):
            s = pm[(pm.disease == disease) & (pm.modality == mod)]
            line.append(f"{mod} {len(s)} types, {int((s.fdr < 0.05).sum())} FDR<0.05")
        if {"bulk", "snRNA"} <= set(w.columns):
            w = w.dropna()
            rho = w.corr(method="spearman").iloc[0, 1]
            line.append(f"Spearman(beta_bulk, beta_snRNA) = {rho:+.2f} over {len(w)} types")
        print(f"  {disease}: " + "; ".join(line))
    rin = M.pivot_table(index=["disease", "modality", "cell_type"], columns="model", values="beta").dropna()
    print(f"\n+RIN sensitivity: Pearson(beta primary, beta +RIN) = "
          f"{np.corrcoef(rin.primary, rin.plus_rin)[0, 1]:.3f} over {len(rin)} estimates")
    print(f"wrote 4 tables to {OUT}")


if __name__ == "__main__":
    main()
