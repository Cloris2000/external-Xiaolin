"""Spearman version of the matrix run's scores.csv (user request 2026-11-10: strategy
selection on Spearman).

retention_matrix.score() stores Pearson r only. This recomputes, for EVERY row of
$STRATEGY_OUT/scores.csv, Spearman rho on exactly the donors and values that row used,
from the saved predictions (no fitting):
  P = pred/<trainer>__<target_bulk>.parquet   (own-cohort rows are out-of-fold values)
  B = oof/<target>.parquet                    (target's own out-of-fold baseline)
  Y = CLR labels of the target snRNA          (RM.load_counts)
  donors = P ∩ B ∩ Y, minus the trainer's training people unless diagonal, where
  training people = persons of the trainer bulk's SOURCE (light) matrix ∩ persons of the
  trainer's label (a '+'-joined consensus label for consensus trainers) -- main()'s rule.
  per cell type: finite in P, B and Y; >= 50 donors; non-constant -- score()'s rule.
Pearson is recomputed alongside on the same values and must reproduce scores.csv's
r_transfer (the check that the donor sets are identical); the script stops if not.

Output: $STRATEGY_OUT/scores_spearman.csv -- scores.csv's key columns with r_transfer
replaced by Spearman rho (r_pearson kept for reference). New file; scores.csv untouched.
Run (env as strategy_scores.py): $PY scores_spearman.py   (single-threaded, a few minutes)
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

MAT = Path(os.environ.get("STRATEGY_OUT",
                          "/scratch/zhoux156/ctp_reliability_paper/celmod_strategies/matrix"))
REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "stages" / "composition"))
from ctcomp.bulk import to_persons, matrix_stem, bank_of, PART_OF  # noqa: E402
import retention_matrix as RM  # noqa: E402

MIN_SCORE = RM.MIN_SCORE
KEEP = ["trainer", "trainer_bulk", "trainer_snrna", "target", "target_bulk", "target_snrna",
        "cell_type", "in_core", "n", "diag", "same_snrna", "same_bulk", "r_transfer"]
S = pd.read_csv(MAT / "scores.csv", usecols=KEEP)
print(f"scores.csv: {len(S)} rows, {S.groupby(['trainer', 'target']).ngroups} trainer x target cells")

labels = RM.load_counts(sorted(set(S.trainer_snrna) | set(S.target_snrna)))
_light = {}


def src_persons(b):
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


train_people = {t: src_persons(tb) & set(labels[tn].index)
                for t, tb, tn in S[["trainer", "trainer_bulk", "trainer_snrna"]].drop_duplicates().itertuples(index=False)}
oof = {u: pd.read_parquet(MAT / "oof" / f"{u}.parquet") for u in S.target.unique()}

out = []
for (tau, ub), g in S.groupby(["trainer", "target_bulk"], sort=True):
    P = pd.read_parquet(MAT / "pred" / f"{tau}__{ub}.parquet")
    for (u, un, diag), h in g.groupby(["target", "target_snrna", "diag"]):
        B, Y = oof[u], labels[un]
        donors = sorted(set(P.index) & set(B.index) & set(Y.index))
        if not diag:
            donors = sorted(set(donors) - train_people[tau])
        Pd, Bd, Yd = P.loc[donors], B.loc[donors], Y.loc[donors]
        for row in h.itertuples(index=False):
            ct = row.cell_type
            a = Pd[ct].values.astype(float)
            b = Bd[ct].values.astype(float)
            y = Yd[ct].values.astype(float)
            ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(y)
            a, y = a[ok], y[ok]
            rec = row._asdict()
            rec["r_pearson"] = float(np.corrcoef(a, y)[0, 1])
            rec["n_check"] = int(ok.sum())
            rec["r_transfer"] = float(stats.spearmanr(a, y).correlation)
            out.append(rec)

R = pd.DataFrame(out)
S_r = S.set_index(["trainer", "target", "cell_type"]).r_transfer
d = (R.set_index(["trainer", "target", "cell_type"]).r_pearson - S_r).abs()
bad_n = (R.n != R.n_check).sum()
print(f"Pearson check: max |r - scores.csv r_transfer| = {d.max():.2e}; donor-count mismatches = {bad_n}")
if d.max() > 1e-8 or bad_n:
    raise SystemExit("donor sets do not reproduce scores.csv -- not writing")
R = R.drop(columns="n_check")
R.to_csv(MAT / "scores_spearman.csv", index=False)
c = R[R.in_core & ~R.diag]
print(f"wrote {len(R)} rows -> {MAT / 'scores_spearman.csv'}")
print(f"Core off-diagonal: median rho {c.r_transfer.median():.3f} vs r {c.r_pearson.median():.3f}; "
      f"corr(rho, r) over rows = {c.r_transfer.corr(c.r_pearson):.3f}")
