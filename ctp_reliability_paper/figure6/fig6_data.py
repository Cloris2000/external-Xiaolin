"""Figure 6 inputs: the production classifier.

CelMod trained on ROSMAP bulk with a consensus label over Green + Mathys +
AMP-AD Diverse multiome (the PI's "Multiome/Diverse"), tag `prod`
(BULK_SRC=batch, BULK_TGT=raw; decision record 2026-09-19). Every table below
is from that build, not the superseded `light` build of 2026-09-18.

PI's tables (read-only):
  195_prod_labels.csv          which snRNA legs each ROSMAP donor's label used
  190_arm_comparison.csv       consensus vs single label at matched n, scored on a
                               held-out snRNA dataset (disattenuated r), 3 families
  196_production_scorecard.csv MSBB harness: production vs single-label references
  197_demoA_cross_region.csv   bulk cross-region agreement vs SEA-AD snRNA prediction
  deconv/*_DFC_supertype_prod__consensus.parquet  the deliverable (donor counts)

Writes data/fig6_*.tsv.
"""
import glob
import os
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

R = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
D = Path("/scratch/shreejoy/cell_type_bias/wp3/deconv")
OUT = Path(__file__).resolve().parent / "data"
GROUPS = {  # gene-space groups of 195_production_model.py
    "rnaseq": ["ROSMAP", "ROSMAP_PCC", "MSBB_FP", "MSBB_IFG", "MSBB_STG", "NIMH_HBCC", "CommonMind", "GVEX"],
    "microarray": ["HBTRC"],
    "gtex_regions": ["GTEx_BA9_frz", "GTEx_BA24", "GTEx_CTX", "GTEx_CBL"],
    "wp1_bulk": ["Diverse_MSSM_DLPFC", "Diverse_MSSM_STG", "Diverse_MayoEmory_DLPFC",
                 "Diverse_MayoEmory_STG", "Diverse_Rush2_DLPFC", "Diverse_Rush2_STG",
                 "Diverse_Columbia_DLPFC", "Mayo_TCX", "NABEC"],
}


def main():
    OUT.mkdir(exist_ok=True)
    L = pd.read_csv(R / "195_prod_labels.csv")
    legs = ["Green_DFC", "Mathys_DFC", "Multiome_DFC"]
    L["combo"] = L[legs].apply(lambda r: " + ".join(l.replace("_DFC", "") for l in legs if r[l]), axis=1)
    a = L.groupby(["n_legs", "combo"]).size().reset_index(name="n_donors")
    a.to_csv(OUT / "fig6a_label_legs.tsv", sep="\t", index=False)
    print(f"a: {len(L)} labelled ROSMAP donors\n{a}")

    g = pd.read_csv(R / "190_arm_comparison.csv")
    g = g[~g.arm.str.contains("consTau")]  # production uses the plain mean
    g["family_lab"] = g.family.str.replace("exclude_", "held out: ")
    g["arm_lab"] = g.arm.str.split("__").str[1].str.replace("cons_", "consensus: ").str.replace("solo_", "single: ").str.replace("_", " + ")
    g.to_csv(OUT / "fig6b_consensus_vs_single.tsv", sep="\t", index=False)
    print(g[g.vs.isna() & (g.scope == "CORE")][["family_lab", "arm_lab", "median_r_disatt"]].round(3))

    c = pd.read_csv(R / "196_production_scorecard.csv")
    c = c[c.arm.notna()].copy()
    c.to_csv(OUT / "fig6c_msbb_harness.tsv", sep="\t", index=False)
    print(c[["arm", "scope", "median_reliability", "frac_sign_ok"]].round(3))

    x = pd.read_csv(R / "197_demoA_cross_region.csv")
    x.to_csv(OUT / "fig6d_cross_region.tsv", sep="\t", index=False)

    rows = []
    for f in sorted(glob.glob(str(D / "*_DFC_supertype_prod__consensus.parquet"))):
        coh = os.path.basename(f).split("_DFC_")[0]
        grp = next((k for k, v in GROUPS.items() if coh in v), "unassigned")
        rows.append(dict(cohort=coh, group=grp, n_donors=pq.ParquetFile(f).metadata.num_rows))
    e = pd.DataFrame(rows)
    assert (e.group != "unassigned").all(), e[e.group == "unassigned"]
    e.to_csv(OUT / "fig6e_deliverable.tsv", sep="\t", index=False)
    print(f"e: {len(e)} cohorts, {e.n_donors.sum()} donors\n{e.groupby('group').n_donors.agg(['size', 'sum'])}")
    print(f"wrote 5 tables to {OUT}")


if __name__ == "__main__":
    main()
