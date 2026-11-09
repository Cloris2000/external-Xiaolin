"""5-fold cross-validated (out-of-fold) scores for the consensus models, to fill the
consensus x member-dataset cells of the Fig 9 matrix (user decision 2026-11-04).

Why: with the person hold-out those cells have 0 test donors (the consensus trained on
every donor of each member dataset). They are filled the way the circled diagonal is:
5-fold out-of-fold on the model's own bulk cohort, then scored against each member
dataset on that member's donors. Like the diagonal, the test dataset is part of the
training label, so these cells are drawn circled and read as same-dataset scores.

How: mirrors retention_matrix.main()'s diagonal block line for line (hvg on the fit
donors, zscore source=light / target=raw in each matrix's own space, out_of_fold with
seed 0 / 5 folds), using Shreejoy's functions unchanged from
/scratch/shreejoy/composition/wt_gtex/stages/composition/retention_matrix.py. The gene
space is the full run's (common_genes over the cohorts in matrix/pairs_used.csv), so the
models match the matrix run. Writes only to $STRATEGY_OUT (/scratch/zhoux156/...).

Output: $STRATEGY_OUT/cons_oof_scores.csv  (trainer, target, cell_type, in_core, n, r_cv)
        $STRATEGY_OUT/cons_oof/<trainer>.parquet  out-of-fold predictions
Run on a compute node: sbatch cons_oof.sbatch
"""
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(os.environ["STRATEGY_OUT"]).resolve()
SAFE = Path("/scratch/zhoux156")
if SAFE not in OUT.parents:
    raise SystemExit(f"refusing: STRATEGY_OUT {OUT} is not under {SAFE}")
os.environ["LABCORE_RESULTS"] = str(OUT / "_labcore")
os.environ.pop("LABCORE_READS", None)
REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "stages" / "composition"))
import retention_matrix as RM  # noqa: E402

CONSENSUS = ["ROSMAP__Green+Mathys+Multiome2025_DLPFC+CUIMC2",
             "ROSMAP__Green+Mathys+Multiome2025_DLPFC",
             "NIMH_HBCC__PsychAD_HBCC+Catching_2026_HBCC"]
T0 = time.time()


def log(m):
    print(f"[cons_oof {time.time() - T0:6.0f}s] {m}", flush=True)


def main():
    (OUT / "cons_oof").mkdir(parents=True, exist_ok=True)
    pu = pd.read_csv(OUT / "pairs_used.csv")
    genes = RM.common_genes(sorted(set(pu.bulk)))
    log(f"{len(genes)} common genes (full-run gene space)")
    core = set(RM.T.core())
    rows = []
    for tr in CONSENSUS:
        b, recipe = tr.split("__", 1)
        members = recipe.split("+")
        src = RM.load_bulk([b], RM.BULK_SRC, genes)[b]
        tgt = RM.load_bulk([b], RM.BULK_TGT, genes)[b]
        clr = RM.load_counts([recipe] + members)
        donors = sorted(set(src.index) & set(tgt.index) & set(clr[recipe].index))
        hv = RM.hvg(src, donors, genes)
        Xs = RM.zscore(src, donors, hv)
        Xt = RM.zscore(tgt, donors, hv)
        P, sc = RM.out_of_fold(Xs, Xt, clr[recipe].loc[donors], donors)
        P.to_parquet(OUT / "cons_oof" / f"{tr}.parquet")
        log(f"{tr}: {len(donors)} donors, {P.shape[1]} types, median r_oof vs own "
            f"consensus label {sc.r_oof.median():.3f}")
        for m in members:
            Y = clr[m]
            D = sorted(set(P.index) & set(Y.index))
            for ct in P.columns.intersection(Y.columns):
                a = P.loc[D, ct].values.astype(float)
                y = Y.loc[D, ct].values.astype(float)
                ok = np.isfinite(a) & np.isfinite(y)
                if ok.sum() < RM.MIN_SCORE or a[ok].std() < 1e-9 or y[ok].std() < 1e-9:
                    continue
                rows.append(dict(trainer=tr, target=f"{b}__{m}", cell_type=ct,
                                 in_core=ct in core, n=int(ok.sum()),
                                 r_cv=float(np.corrcoef(a[ok], y[ok])[0, 1])))
            log(f"  scored vs {m}: {len(D)} donors")
    S = pd.DataFrame(rows)
    S.to_csv(OUT / "cons_oof_scores.csv", index=False)
    log(f"wrote {len(S)} rows -> {OUT/'cons_oof_scores.csv'}")


if __name__ == "__main__":
    main()
