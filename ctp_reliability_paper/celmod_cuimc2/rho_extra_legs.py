"""Spearman rho (out-of-fold prediction vs observed snRNA-seq CLR) for the extra Fig 4 legs.

The legs built/used so far saved only Pearson r_oof. This recomputes, from the saved
out-of-fold prediction parquets, the same r_oof AND the Spearman rho on the same matched
donors and the same CLR targets (build_legs.clr_within, unchanged). It asserts that the
recomputed Pearson r equals the saved r_oof, so rho is on exactly the donors/targets of r.
ROSMAP x Green is included as a cross-check against the wp3 115 table's rho column.
Writes OUT/rho_extra_legs.csv (scratch). No CelMod is fitted.
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

REPO = Path(os.environ["CTCOMP_REPO"]); sys.path.insert(0, str(REPO / "stages" / "deconvolution"))
import build_legs as BL  # noqa: E402
from build_reference import compartments  # noqa: E402

RES = Path(os.environ["CTCOMP_RESULTS"])
OUT = Path(os.environ.get("RHO_OUT", "/scratch/zhoux156/ctp_reliability_paper/celmod_cuimc2"))
CUIMC2 = Path("/scratch/zhoux156/ctp_reliability_paper/celmod_cuimc2")
# (cohort, trainer, level, prediction parquet, saved r_oof table)
LEGS = [("ROSMAP", "CUIMC2", "supertype", CUIMC2 / "ROSMAP_DFC_supertype_oof__CUIMC2.parquet", CUIMC2 / "oof_vs_insample__CUIMC2.csv"),
        ("ROSMAP", "CUIMC2", "subclass", CUIMC2 / "ROSMAP_DFC_subclass_oof__CUIMC2.parquet", CUIMC2 / "oof_vs_insample__CUIMC2.csv"),
        ("NABEC", "Catching_2026_NABEC", "supertype", RES / "deconv/NABEC_DFC_supertype_oof__Catching_2026_NABEC.parquet", RES / "deconv/oof_vs_insample.csv"),
        ("ROSMAP", "Green", "supertype", RES / "deconv/ROSMAP_DFC_supertype_oof__Green.parquet", RES / "deconv/oof_vs_insample.csv")]
rows = []
for cohort, trainer, level, pp, sp in LEGS:
    BL.LEVEL = level
    B = BL.read_target(os.environ["CTCOMP_BULK"], cohort, BL.BULK_SRC); B.index = B.index.astype(str)
    C = pd.read_parquet(BL.COUNTS / f"{trainer}__{level}.parquet"); C.index = C.index.astype(str)
    train = sorted(set(B.index) & set(C.index))
    Y = BL.clr_within(C.loc[train], compartments(level))
    P = pd.read_parquet(pp); P.index = P.index.astype(str)
    saved = pd.read_csv(sp).query("cohort == @cohort and trainer == @trainer and level == @level").set_index("cell_type")
    for ct in P.columns:
        y, p = Y[ct].values.astype(float), P.loc[train, ct].values.astype(float)
        r = float(np.corrcoef(p, y)[0, 1]); rho = float(spearmanr(p, y).correlation)
        assert abs(r - saved.loc[ct, "r_oof"]) < 1e-6, (cohort, trainer, level, ct, r, saved.loc[ct, "r_oof"])
        rows.append(dict(cohort=cohort, trainer=trainer, level=level, cell_type=ct, n_train=len(train), r_oof=r, rho_oof=rho))
    print(f"{cohort} x {trainer} {level}: {len(train)} donors, {len(P.columns)} cell types; Pearson r reproduces saved r_oof", flush=True)
pd.DataFrame(rows).to_csv(OUT / "rho_extra_legs.csv", index=False)
