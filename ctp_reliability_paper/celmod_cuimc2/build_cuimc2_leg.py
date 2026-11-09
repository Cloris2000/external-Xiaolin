"""CelMod out-of-fold leg ROSMAP bulk x CUIMC2 snRNA, built with Shreejoy's own code.

Why: CUIMC2 joins the paper's snRNA pairs (2026-10-29) and has no CelMod leg in any
composition results root (checked results_2026-10-12b ... 2026-10-21c). The handoff
(2026-10-28) lists CUIMC2 x bulk_ROSMAP at 239 shared people (DLPFC), above
build_legs.MIN_PAIRED = 100.

How: imports stages/deconvolution/build_legs.py from Shreejoy's composition checkout
(/scratch/shreejoy/composition/wt_gtex, the checkout that built the 2026-10-21c legs,
per results_2026-10-21c/deconv_legs/DONE.json) and calls its build_one() unchanged:
  - bulk:   ROSMAP DLPFC, BULK_SRC=light (as the 2026-10-21c legs), from CTCOMP_BULK
  - labels: results_2026-10-21c/counts/CUIMC2__<level>.parquet (registry person_id)
  - CelMod: N_HVG 15000, 5-fold OOF, SEED 0, CLR within compartment, as every leg
Nothing is written into Shreejoy's results root: build_one() returns data frames and
this script writes them to OUT (scratch). LABCORE_READS is left unset, so labcore
writes no read-tracking files either.

Steps
  1. validation: rebuild the existing ROSMAP x Multiome2025_DLPFC supertype leg (the
     smallest ROSMAP leg, 269 paired) and compare with Shreejoy's saved parquet and
     oof_vs_insample.csv; stop if it does not reproduce (tolerance VAL_TOL).
  2. ROSMAP x CUIMC2, supertype.
  3. ROSMAP x CUIMC2, subclass (Fig 4b uses subclass out-of-fold r as well).

Outputs (OUT):
  ROSMAP_DFC_<level>_oof__CUIMC2.parquet   donor x cell type predictions (OOF for the
                                            paired donors, out-of-sample for the rest)
  oof_vs_insample__CUIMC2.csv              per level x cell type: n_train, r_oof, r_insample
  validation_ROSMAP_x_Multiome2025_DLPFC.csv  max |diff| of the rebuilt vs saved leg

Run on a compute node only (CelMod exhausts the login node's thread budget):
  sbatch build_cuimc2_leg.sbatch
"""
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "stages" / "deconvolution"))
import build_legs as BL  # noqa: E402
from build_reference import compartments  # noqa: E402

OUT = Path(os.environ.get("CUIMC2_LEG_OUT", "/scratch/zhoux156/ctp_reliability_paper/celmod_cuimc2"))
SAVED = Path(os.environ["CTCOMP_RESULTS"]) / "deconv"
COHORT, TRAINER = "ROSMAP", "CUIMC2"
VAL_TRAINER = "Multiome2025_DLPFC"
VAL_TOL = 1e-6   # threaded BLAS may reorder sums; anything larger is a real difference
T0 = time.time()


def log(m):
    print(f"[cuimc2 {time.time() - T0:6.0f}s] {m}", flush=True)


def leg(trainer, level):
    BL.LEVEL = level   # build_one reads the module-level LEVEL for the counts file name
    pred, sc = BL.build_one(COHORT, trainer, compartments(level))
    if pred is None:
        raise SystemExit(f"{COHORT} x {trainer} {level}: not built (see log above)")
    return pred, sc


def validate():
    pred, sc = leg(VAL_TRAINER, "supertype")
    ref_p = pd.read_parquet(SAVED / f"{COHORT}_DFC_supertype_oof__{VAL_TRAINER}.parquet")
    ref_s = pd.read_csv(SAVED / "oof_vs_insample.csv")
    ref_s = ref_s[(ref_s.cohort == COHORT) & (ref_s.trainer == VAL_TRAINER)]
    assert set(pred.index) == set(ref_p.index) and set(pred.columns) == set(ref_p.columns), \
        "validation: donor or cell-type sets differ from the saved leg"
    d_pred = np.nanmax(np.abs(pred.loc[ref_p.index, ref_p.columns].values - ref_p.values))
    m = sc.merge(ref_s, on="cell_type", suffixes=("", "_saved"))
    d_r = np.nanmax(np.abs(m.r_oof - m.r_oof_saved))
    res = pd.DataFrame([dict(leg=f"{COHORT} x {VAL_TRAINER}", n_train=int(sc.n_train.iloc[0]),
                             n_train_saved=int(ref_s.n_train.iloc[0]),
                             max_abs_diff_pred=d_pred, max_abs_diff_r_oof=d_r)])
    res.to_csv(OUT / f"validation_{COHORT}_x_{VAL_TRAINER}.csv", index=False)
    log(f"validation {COHORT} x {VAL_TRAINER}: max |diff| predictions {d_pred:.2e}, r_oof {d_r:.2e}")
    if max(d_pred, d_r) > VAL_TOL or res.n_train[0] != res.n_train_saved[0]:
        raise SystemExit("validation failed: the rebuilt leg does not reproduce Shreejoy's saved leg")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"CTCOMP_REPO={REPO}  CTCOMP_RESULTS={os.environ['CTCOMP_RESULTS']}  "
        f"BULK_SRC={BL.BULK_SRC}  OUT={OUT}")
    validate()
    scores = []
    for level in ("supertype", "subclass"):
        pred, sc = leg(TRAINER, level)
        pred.to_parquet(OUT / f"{COHORT}_DFC_{level}_oof__{TRAINER}.parquet")
        scores.append(sc)
        n_tr = int(sc.n_train.iloc[0])
        log(f"{COHORT} x {TRAINER} {level}: {n_tr} paired donors, {sc.shape[0]} cell types, "
            f"median r_oof {sc.r_oof.median():.3f} (in-sample {sc.r_insample.median():.3f})")
    pd.concat(scores, ignore_index=True).to_csv(OUT / f"oof_vs_insample__{TRAINER}.csv", index=False)
    log(f"wrote {OUT}")


if __name__ == "__main__":
    main()
