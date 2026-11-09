"""Spearman version of the Fig 9 trainer x test heatmap (user request 2026-11-04).

The matrix run (retention_matrix.score) stores Pearson r only. This recomputes, from the
SAVED predictions, Spearman rho for every heatmap cell on exactly the donor set the
Pearson score used:
  * off-diagonal: pred/<trainer>__<bulk>.parquet; donors = prediction ∩ target's own
    out-of-fold baseline ∩ test snRNA labels, minus every person the trainer trained on
    (retention_matrix.score's rule); >= 50 finite donors per cell type
  * diagonal: the same prediction file (its own-cohort rows are the out-of-fold values)
  * consensus x member cells: cons_oof/<trainer>.parquet (5-fold out-of-fold),
    donors = prediction ∩ member labels (cons_oof.py's rule)
Summary per cell = median over the 86 Core Supertypes, as for the Pearson heatmap.
No fitting; login-node safe.

Output: ../figure9/data/fig9_heatmap_spearman.tsv
Run (env as strategy_scores.py): $PY heatmap_spearman.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
MAT = Path(os.environ.get("STRATEGY_OUT",
                          "/scratch/zhoux156/ctp_reliability_paper/celmod_strategies/matrix"))
FIG = HERE.parent / "figure9" / "data"
REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "stages" / "composition"))
from ctcomp.bulk import to_persons, matrix_stem, bank_of, PART_OF  # noqa: E402
import retention_matrix as RM  # noqa: E402

MIN_SCORE = 50
REP_BULK = ["ROSMAP", "MSBB_FP", "Diverse_Rush2_DLPFC", "CommonMind_MSSM", "NIMH_HBCC", "NABEC"]
PAPER_SN = {"Green", "Mathys", "Multiome2025_DLPFC", "CUIMC2", "PsychAD_MSBB",
            "PsychAD_HBCC", "PsychAD_ROSMAP", "BrainSCOPE_CMC", "Catching_2026_HBCC",
            "Catching_2026_NABEC"}
CONSENSUS = ["ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2",
             "ROSMAP__Green+Mathys+Multiome2025_DLPFC",
             "NIMH_HBCC__PsychAD_HBCC+Catching_2026_HBCC"]


def norm(ct):
    return ct.replace("/", "_").replace(" ", "_")


ci = pd.read_csv(HERE.parent / "common" / "data" / "reliability_ci_summary.tsv", sep="\t")
ci = ci[ci.level == "supertype"]
CORE = set(ci.loc[ci.ci_lo > 0, "cell_type"].map(norm))

pu = pd.read_csv(MAT / "pairs_used.csv")
singles = pu[(pu.bulk.isin(REP_BULK)) & (pu.snrna.isin(PAPER_SN))]
trainers = sorted(set(singles.loc[singles.role == "trainer", "pair"]) | set(CONSENSUS))
targets = sorted(set(singles.loc[singles.role == "target", "pair"]))

_light = {}


def src_persons(b):
    """persons of b's light (source) matrix, as retention_matrix resolves them"""
    if b not in _light:
        idx = pd.read_parquet(Path(os.environ["CTCOMP_BULK"]) / f"{matrix_stem(b)}_light.parquet",
                              columns=[]).index.astype(str)
        p, _ = to_persons(idx, matrix_stem(b))
        s = set(p.dropna())
        if b in PART_OF:
            banks = bank_of(matrix_stem(b))
            s = {x for x in s if banks.get(x) == PART_OF[b][1]}
        _light[b] = s
    return _light[b]


members = {t: t.split("__", 1)[1].split("+") for t in trainers}
labels = RM.load_counts(sorted({m for ms in members.values() for m in ms}
                               | {t.split("__", 1)[1] for t in targets}))
train_people = {t: src_persons(t.split("__")[0]) &
                set().union(*[set(labels[m].index) for m in members[t]]) for t in trainers}


def med_rho(P, Y, D):
    rh = []
    for ct in P.columns.intersection(Y.columns):
        if ct not in CORE:
            continue
        a = P[ct].reindex(D).values.astype(float)
        y = Y[ct].reindex(D).values.astype(float)
        ok = np.isfinite(a) & np.isfinite(y)
        if ok.sum() < MIN_SCORE or a[ok].std() < 1e-9 or y[ok].std() < 1e-9:
            continue
        rho = stats.spearmanr(a[ok], y[ok]).correlation
        rh.append(rho)
        PER_TYPE.append(dict(cell_type=ct, rho=rho, n=int(ok.sum())))
    return (float(np.median(rh)) if rh else np.nan), len(rh)


rows, per_rows = [], []
PER_TYPE = []
for tau in trainers:
    tb = tau.split("__")[0]
    for u in targets:
        ub, us = u.split("__", 1)
        Y = labels[us]
        diag = tau == u
        cons_cell = tau in CONSENSUS and tb == ub and us in members[tau]
        if cons_cell:
            f = MAT / "cons_oof" / f"{tau}.parquet"
            if not f.exists():
                continue
            P = pd.read_parquet(f)
            D = sorted(set(P.index) & set(Y.index))
        else:
            f, fb = MAT / "pred" / f"{tau}__{ub}.parquet", MAT / "oof" / f"{u}.parquet"
            if not (f.exists() and fb.exists()):
                continue
            P = pd.read_parquet(f)
            B = pd.read_parquet(fb)
            D = set(P.index) & set(B.index) & set(Y.index)
            if not diag:
                D -= train_people[tau]
            D = sorted(D)
        if len(D) < MIN_SCORE:
            continue
        PER_TYPE.clear()
        m, k = med_rho(P, Y, D)
        per_rows += [dict(trainer=tau, target=u, diag=diag, cv=cons_cell, **r) for r in PER_TYPE]
        if k:
            rows.append(dict(trainer=tau, target=u, diag=diag, cv=cons_cell,
                             med_core=m, n=len(D), n_types=k))
R = pd.DataFrame(rows)
R.to_csv(FIG / "fig9_heatmap_spearman.tsv", sep="\t", index=False)
# per (cell, Core Supertype) rho, for the supertype x cohort panel (user 2026-11-09)
pd.DataFrame(per_rows).to_csv(FIG / "fig9_heatmap_spearman_celltype.tsv", sep="\t", index=False)
print(f"wrote {len(R)} cells; median of cell medians {R.med_core.median():.3f}")

# agreement check with the Pearson heatmap on the same cells
pr = pd.read_csv(FIG / "fig9_trainer_target.tsv", sep="\t")
pr["cv"] = pr.get("cv", False)
m = R.merge(pr[["trainer", "target", "med_core"]], on=["trainer", "target"],
            suffixes=("_rho", "_r"))
print(f"{len(m)} cells matched to the Pearson heatmap; corr(rho, r) = "
      f"{m.med_core_rho.corr(m.med_core_r):.3f}; median rho - r = "
      f"{(m.med_core_rho - m.med_core_r).median():+.3f}")
