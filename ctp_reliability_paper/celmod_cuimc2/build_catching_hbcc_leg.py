"""EXPLORATORY CelMod out-of-fold leg NIMH_HBCC bulk x Catching_2026_HBCC snRNA.

Same code path as build_cuimc2_leg.py (Shreejoy's build_legs.build_one(), unchanged),
except MIN_PAIRED is lowered from 100 to MIN_PAIRED_HERE in THIS process only (module
attribute set after import; his checkout is not touched). The pair has 89 shared
donors (results_2026-10-21c/reliability/pair_overlap.csv), so build_one() refuses it
at the stock minimum. User decision 2026-11-10: build it to see the result before
deciding whether it goes in Fig 4b. Every output carries n_train; treat as a
sensitivity result, not a pipeline leg. Writes to scratch only.

Outputs (OUT): NIMH_HBCC_DFC_<level>_oof__Catching_2026_HBCC.parquet, oof_vs_insample__Catching_2026_HBCC.csv
"""
import os
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "stages" / "deconvolution"))
import build_legs as BL  # noqa: E402
from build_reference import compartments  # noqa: E402

MIN_PAIRED_HERE = 80
BL.MIN_PAIRED = MIN_PAIRED_HERE
OUT = Path(os.environ.get("CATCHING_HBCC_LEG_OUT", "/scratch/zhoux156/ctp_reliability_paper/celmod_catching_hbcc"))
COHORT, TRAINER = "NIMH_HBCC", "Catching_2026_HBCC"
T0 = time.time()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"MIN_PAIRED lowered 100 -> {BL.MIN_PAIRED} (exploratory); BULK_SRC={BL.BULK_SRC}", flush=True)
    scores = []
    for level in ("supertype", "subclass"):
        BL.LEVEL = level
        pred, sc = BL.build_one(COHORT, TRAINER, compartments(level))
        if pred is None:
            raise SystemExit(f"{COHORT} x {TRAINER} {level}: not built")
        pred.to_parquet(OUT / f"{COHORT}_DFC_{level}_oof__{TRAINER}.parquet")
        scores.append(sc)
        print(f"[{time.time() - T0:5.0f}s] {level}: n_train={int(sc.n_train.iloc[0])}, {sc.shape[0]} cell types, "
              f"median r_oof {sc.r_oof.median():.3f} (in-sample {sc.r_insample.median():.3f})", flush=True)
    pd.concat(scores, ignore_index=True).to_csv(OUT / f"oof_vs_insample__{TRAINER}.csv", index=False)


if __name__ == "__main__":
    main()
