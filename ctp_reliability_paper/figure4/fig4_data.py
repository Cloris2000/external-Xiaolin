"""Figure 4 inputs: deconvolution resolution is better at supertype level.

All inputs are the PI's DFC 2026 tables (read-only):
  117_whitelist_dfc.csv            median snRNA ICC (6 independent pairs) and median
                                   out-of-fold CelMod r (6 bulk x snRNA pairs) per supertype
  196_subclass_vs_supertype.csv    triangulation reliability in the ROSMAP triad, both levels
  115_celmod_dfc_{level}.csv       per (bulk x snRNA pair, cell type) out-of-fold CelMod r
  185_portable_consensus.csv       CORE flag (34 portable supertypes)
  taxonomy_DFC_2026.tsv            supertype -> subclass, compartment, reactive_state

Writes data/fig4_*.tsv.
"""
from pathlib import Path
import pandas as pd

R = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv")
OUT = Path(__file__).resolve().parent / "data"


def main():
    OUT.mkdir(exist_ok=True)
    tax = pd.read_csv(TAX, sep="\t")
    sub_of = dict(zip(tax.supertype_label, tax.subclass_label))
    comp = dict(zip(tax.supertype_label, tax.compartment))
    react = dict(zip(tax.supertype_label, tax.reactive_state))
    core = set(pd.read_csv(R / "185_portable_consensus.csv").query("CORE == True").cell_type)

    # a: deconvolution accuracy vs snRNA reliability, one row per supertype
    w = pd.read_csv(R / "117_whitelist_dfc.csv")
    a = w[["cell_type", "median_icc", "n_pairs", "bulk_r", "bulk_n_pairs"]].copy()
    a["compartment"] = a.cell_type.map(comp)
    a["reactive"] = a.cell_type.map(react)
    a["core"] = a.cell_type.isin(core)
    assert a.compartment.notna().all()
    a.to_csv(OUT / "fig4a_icc_vs_deconv.tsv", sep="\t", index=False)

    # b: subclass reliability vs its own child supertypes (ROSMAP triad, bulk leg)
    s = pd.read_csv(R / "196_subclass_vs_supertype.csv")
    s = s[s.leg == "bulk"] if "bulk" in set(s.leg) else s
    s = s.drop_duplicates(["level", "cell_type"])
    sup = s[s.level == "supertype"].copy()
    sup["subclass"] = sup.cell_type.map(sub_of)
    sub = s[s.level == "subclass"].rename(columns={"cell_type": "subclass"})
    b = sup.merge(sub[["subclass", "reliability", "portability", "sign_ok"]],
                  on="subclass", suffixes=("", "_subclass"))
    b["resolved"] = (b.portability >= 0.975) & b.sign_ok
    b["resolved_subclass"] = (b.portability_subclass >= 0.975) & b.sign_ok_subclass
    b["core"] = b.cell_type.isin(core)
    b.to_csv(OUT / "fig4b_subclass_vs_children.tsv", sep="\t", index=False)
    n_sub = b.subclass.nunique()
    med_child = b.groupby("subclass").reliability.median()
    rel_sub = b.groupby("subclass").reliability_subclass.first()
    print(f"b: {n_sub} subclasses; subclass below median child in {(rel_sub < med_child).sum()}/{n_sub}; "
          f"resolved subclasses {b.groupby('subclass').resolved_subclass.first().sum()}/{n_sub}, "
          f"resolved supertypes {b.resolved.sum()}/{len(b)}")

    # c: out-of-fold r by resolution, per bulk x snRNA pair
    S = pd.read_csv(R / "115_celmod_dfc_supertype.csv")
    B = pd.read_csv(R / "115_celmod_dfc_subclass.csv")
    S["set"] = ["CORE supertypes" if c in core else "Other supertypes" for c in S.cell_type]
    B["set"] = "Subclasses"
    c = pd.concat([S, B])[["pair", "bulk", "snRNA", "level", "cell_type", "set", "n_donors", "r"]]
    c.to_csv(OUT / "fig4c_oof_by_resolution.tsv", sep="\t", index=False)
    print(c.groupby(["pair", "set"]).r.median().unstack().round(3))
    print(f"wrote 3 tables to {OUT}")


if __name__ == "__main__":
    main()
