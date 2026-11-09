"""Figure 5 inputs: transferring CelMod across cohorts.

PI's DFC 2026 tables (read-only):
  192_retention_clean.csv    3x3 (+pooled) retention matrix: median transferred r /
                             median same-cohort r, CORE and all scope; in-sample
                             diagonal excluded by the PI
  192_retention_matrix.csv   per (target, source, supertype) r_transfer, r_same
  191_arm0_vs_celmod.csv     reference-profile NNLS vs CelMod, ROSMAP triad
  176_training_cohort_bias.csv  two bulk legs sharing a training cohort: inflation
  185_portable_consensus.csv CORE flag

Writes data/fig5_*.tsv.
"""
from pathlib import Path
import numpy as np
import pandas as pd

R = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
OUT = Path(__file__).resolve().parent / "data"
# donor-bootstrap 95% CI of the median symmetric ICC (common/reliability_ci.py), as in figure4/fig4_data.py
RELIABILITY_CI = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/common/data/reliability_ci_summary.tsv")
ICC_THRESH = 0.0
SET_HI, SET_LO = "Core Supertype", "Non-core Supertype"
N_BOOT = 2000
SRC ={"srcROSMAP": "ROSMAP", "srcMSBB": "MSBB", "srcHBCC": "NIMH HBCC", "POOL_MSBB_HBCC": "MSBB + HBCC"}
TGT = {"ROSMAP": "ROSMAP", "MSBB_IFG": "MSBB IFG", "MSBB_STG": "MSBB STG", "NIMH_HBCC": "NIMH HBCC"}


def main():
    OUT.mkdir(exist_ok=True)
    core = set(pd.read_csv(R / "185_portable_consensus.csv").query("CORE == True").cell_type)

    a = pd.read_csv(R / "192_retention_clean.csv")
    a["source_lab"] = a.source.map(SRC)
    a["target_lab"] = a.target.map(TGT)
    assert a.source_lab.notna().all() and a.target_lab.notna().all()
    a.to_csv(OUT / "fig5a_retention.tsv", sep="\t", index=False)
    print(a.pivot_table(index=["scope", "source_lab"], columns="target_lab", values="retention").round(2))

    m = pd.read_csv(R / "192_retention_matrix.csv")
    m["source_lab"] = m.source.map(SRC)
    m["target_lab"] = m.target.map(TGT)
    m["core"] = m.cell_type.isin(core)
    # Fig 4 grouping (reliability-only, no bulk data): Core Supertype = lower bound of the
    # donor-bootstrap 95% CI of the median symmetric ICC > 0. Supertypes without an ICC -> NA.
    ci = pd.read_csv(RELIABILITY_CI, sep="\t").query("level == 'supertype'")
    core_sup = dict(zip(ci.cell_type, ci.ci_lo > ICC_THRESH))
    m["core_supertype"] = m.cell_type.map(core_sup)
    print(f"b: {m.core_supertype.notna().groupby(m.cell_type).first().sum()} of {m.cell_type.nunique()} supertypes have an ICC; "
          f"no ICC: {sorted(m[m.core_supertype.isna()].cell_type.unique())}")
    m.to_csv(OUT / "fig5b_retention_by_type.tsv", sep="\t", index=False)

    # retention by Fig 4 group, same pairs as the PI's clean matrix (in-sample pairs excluded).
    # Check first that the ratio of medians reproduces the PI's 'all' retention.
    pairs = a[a.scope == "all"][["source", "target", "retention"]]
    mm = m.merge(pairs[["source", "target"]], on=["source", "target"]).dropna(subset=["r_transfer", "r_same"])
    chk = mm.groupby(["source", "target"]).apply(lambda g: g.r_transfer.median() / g.r_same.median(), include_groups=False)
    chk = pairs.set_index(["source", "target"]).retention.to_frame().join(chk.rename("recomputed"))
    print(chk.round(3))
    assert (chk.retention - chk.recomputed).abs().max() < 0.01, "recomputed retention does not match the PI's"
    g = mm.dropna(subset=["core_supertype"]).groupby(["source", "target", "core_supertype"]).agg(
        n_types=("cell_type", "nunique"), med_transfer=("r_transfer", "median"), med_same=("r_same", "median")).reset_index()
    g["retention"] = g.med_transfer / g.med_same
    g["scope"] = g.core_supertype.map({True: SET_HI, False: SET_LO})
    g["source_lab"] = g.source.map(SRC)
    g["target_lab"] = g.target.map(TGT)
    g.to_csv(OUT / "fig5_retention_by_set.tsv", sep="\t", index=False)
    print(g.pivot_table(index=["scope", "source_lab"], columns="target_lab", values="retention").round(2))

    # 95% CI of Core Supertype retention per pair: resample supertypes with replacement,
    # ratio of medians recomputed on each resample (same supertypes in numerator and denominator)
    rng = np.random.default_rng(0)
    rows = []
    for (s, t), d in mm[mm.core_supertype == True].groupby(["source", "target"]):
        rt, rs = d.r_transfer.to_numpy(), d.r_same.to_numpy()
        idx = rng.integers(0, len(d), size=(N_BOOT, len(d)))
        boot = np.median(rt[idx], axis=1) / np.median(rs[idx], axis=1)
        rows.append({"source": s, "target": t, "n_types": len(d),
                     "retention": np.median(rt) / np.median(rs),
                     "ci_lo": np.quantile(boot, 0.025), "ci_hi": np.quantile(boot, 0.975)})
    ci_ret = pd.DataFrame(rows)
    ci_ret["source_lab"] = ci_ret.source.map(SRC)
    ci_ret["target_lab"] = ci_ret.target.map(TGT)
    ci_ret.to_csv(OUT / "fig5_core_retention_ci.tsv", sep="\t", index=False)
    print(ci_ret.round(2))

    c = pd.read_csv(R / "191_arm0_vs_celmod.csv")
    c.to_csv(OUT / "fig5c_profile_vs_celmod.tsv", sep="\t", index=False)
    print(c.round(3))

    d = pd.read_csv(R / "176_training_cohort_bias.csv").rename(columns={"Unnamed: 0": "cell_type"})
    d["core"] = d.cell_type.isin(core)
    d.to_csv(OUT / "fig5d_trainer_inflation.tsv", sep="\t", index=False)
    print(f"d: n={len(d)} median inflation {d.inflation.median():+.3f}, "
          f"CORE {d[d.core].inflation.median():+.3f}")
    print(f"wrote 4 tables to {OUT}")


if __name__ == "__main__":
    main()
