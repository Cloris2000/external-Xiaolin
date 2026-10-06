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
import pandas as pd

R = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
OUT = Path(__file__).resolve().parent / "data"
SRC = {"srcROSMAP": "ROSMAP", "srcMSBB": "MSBB", "srcHBCC": "NIMH HBCC", "POOL_MSBB_HBCC": "MSBB + HBCC"}
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
    m.to_csv(OUT / "fig5b_retention_by_type.tsv", sep="\t", index=False)

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
