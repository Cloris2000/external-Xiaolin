"""Figure 4 inputs: snRNA-seq reliability predicts deconvolution accuracy, and a
reliability-only threshold (no bulk data used in selection) picks out a
deconvolution-ready set of supertypes.

Redesign (2026-10-14, user-approved after independent verification of a pasted
analysis proposal): panel a now uses the paper's decided reliability metric
(symmetric ICC, 5 matched-donor pairs, common/data/reliability_summary.tsv)
instead of the old 6-pair icc_person_a-based median_icc from 117_whitelist_dfc.
Old panel b (ROSMAP-triad triangulation reliability) is dropped. New panel b
groups supertypes by snRNA-seq reproducibility alone (median ICC > 0, i.e.
reproducible at all, chosen with no bulk data, unlike CORE; user decision
2026-10-21, was >= 0.3 on 2026-10-14) vs ICC <= 0 supertypes vs subclasses,
and shows out-of-fold CelMod Spearman rho by group per bulk x snRNA pair -- the direct
test of the panel a relationship. NB at ICC > 0 the split is 134 vs 15
supertypes and the Mann-Whitney contrasts are much weaker than at >= 0.3
(12 vs 137, all p <= 5e-5): vs subclasses most pairs are n.s. A supplementary table (fig4_supp_rho_by_pair.tsv; was
panel c briefly) shows the panel a relationship (ICC vs bulk_r) holds per
bulk x snRNA pair, not just in the pooled median (Spearman rho per pair).
The 12 concordant supertypes are 11/12 of CORE (one new: Micro-PVM_3-SEAAD) --
CORE remains the production/portability call (also requires cross-cohort bulk
portability); this figure's "concordant" set is reliability-only.

All inputs are the PI's DFC 2026 tables (read-only):
  common/data/reliability_summary.tsv   symmetric ICC (8 pairs since 2026-10-29) per supertype/subclass
  115_celmod_dfc_{level}.csv            per (bulk x snRNA pair, cell type) out-of-fold CelMod r
  185_portable_consensus.csv            CORE flag (34 portable supertypes) -- reported, not used for grouping
  taxonomy_DFC_2026.tsv                 supertype -> subclass, compartment, reactive_state, class_id

Metric (user decision 2026-11-11): out-of-fold SPEARMAN rho between held-out CelMod predictions
and observed snRNA-seq CLR across matched donors (was Pearson r; r is still carried in the tables).
wp3 115 legs use their saved `rho` column; the two extra legs use celmod_cuimc2/data/rho_extra_legs.csv
(rho_extra_legs.py, asserted to sit on the saved r_oof).

Writes data/fig4_*.tsv.
"""
import re
import numpy as np
from pathlib import Path
import pandas as pd
from scipy import stats

R = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv")
RELIABILITY = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/common/data/reliability_summary.tsv")
# donor-bootstrap 95% CI of the median ICC across the 8 pairs (common/reliability_ci.py)
RELIABILITY_CI = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/common/data/reliability_ci_summary.tsv")
# Two extra bulk x snRNA out-of-fold legs from the 2026-10-21c registry build (user decision
# 2026-11-11): NABEC x Catching (Shreejoy's saved oof_vs_insample.csv) and ROSMAP x CUIMC2
# (celmod_cuimc2/build_cuimc2_leg.py, validated against his ROSMAP x Multiome leg). NIMH_HBCC x
# Catching HBCC (89 donors, overfits) is deliberately NOT included, pending the PI.
OOF_2610 = Path("/scratch/shreejoy/composition/results_2026-10-21c/deconv/oof_vs_insample.csv")
RHO_EXTRA = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/celmod_cuimc2/data/rho_extra_legs.csv")
OOF_CUIMC2 = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/celmod_cuimc2/data/oof_vs_insample__CUIMC2.csv")
EXTRA_LEGS = [("NABEC x Catching", "NABEC", "Catching_2026_NABEC", "NABEC", OOF_2610),
              ("ROSMAP x CUIMC2", "ROSMAP", "CUIMC2", "ROSMAP", OOF_CUIMC2)]
OUT = Path(__file__).resolve().parent / "data"
# Core Supertype = lower bound of the donor-bootstrap 95% CI of the median symmetric ICC
# (8 pairs since 2026-10-29, common/reliability_ci.py) > 0 (user decision 2026-10-27; was median ICC > 0)
ICC_THRESH = 0.0
SET_HI, SET_LO, SET_SUB = "Core Supertype", "Non-core Supertype", "Subclasses"


def add_extra_legs(S, B, tax):
    """Append the extra legs, with the registry build's underscore cell-type names mapped back
    to the taxonomy labels used everywhere else ('L2_3_IT_1' -> 'L2/3 IT_1', 'VLMC_Perivascular' -> 'VLMC & Perivascular')."""
    norm = lambda x: re.sub(r"[^0-9A-Za-z-]+", "_", x).strip("_")
    out = {"supertype": [S], "subclass": [B]}
    rho_x = pd.read_csv(RHO_EXTRA)
    for pair, bulk, sn, cohort, path in EXTRA_LEGS:
        o = pd.read_csv(path)
        for level, labels in (("supertype", tax.supertype_label), ("subclass", tax.subclass_label)):
            lab = {norm(x): x for x in labels.unique()}
            assert len(lab) == labels.nunique(), "name normalisation collides"
            x = o[(o.cohort == cohort) & (o.trainer == sn) & (o.level == level)].copy()
            if not len(x):
                assert level == "subclass", (pair, level)   # not drawn; no subclass leg was built
                print(f"extra leg {pair} subclass: none available, skipped")
                continue
            miss = set(x.cell_type) - set(lab)
            assert not miss, (pair, level, sorted(miss)[:5])
            rho = rho_x[(rho_x.cohort == cohort) & (rho_x.trainer == sn) & (rho_x.level == level)].set_index("cell_type")
            assert set(rho.index) == set(x.cell_type) and (rho.r_oof.reindex(x.cell_type).values == rho.r_oof.reindex(x.cell_type).values).all()
            assert np.allclose(rho.r_oof.reindex(x.cell_type).values, x.r_oof.values, atol=1e-6), "rho table is not on the saved r_oof"
            x["rho"] = rho.rho_oof.reindex(x.cell_type).values
            x["cell_type"] = x.cell_type.map(lab)
            out[level].append(pd.DataFrame({"pair": pair, "bulk": bulk, "snRNA": sn, "level": level,
                                            "cell_type": x.cell_type, "n_donors": x.n_train, "r": x.r_oof, "rho": x.rho}))
            print(f"extra leg {pair} {level}: {len(x)} cell types, n_train {int(x.n_train.iloc[0])}, "
                  f"median r_oof {x.r_oof.median():.3f}")
    return pd.concat(out["supertype"], ignore_index=True), pd.concat(out["subclass"], ignore_index=True)


def main():
    OUT.mkdir(exist_ok=True)
    tax = pd.read_csv(TAX, sep="\t")
    comp = dict(zip(tax.supertype_label, tax.compartment))
    react = dict(zip(tax.supertype_label, tax.reactive_state))
    class_of = dict(zip(tax.supertype_label, tax.class_id))
    core = set(pd.read_csv(R / "185_portable_consensus.csv").query("CORE == True").cell_type)

    rel = pd.read_csv(RELIABILITY, sep="\t")
    sup_rel = rel[rel.level == "supertype"][["cell_type", "median_icc_sym", "n_pairs", "n_pairs_ge_030"]]
    ci = pd.read_csv(RELIABILITY_CI, sep="\t").query("level == 'supertype'")[["cell_type", "ci_lo", "ci_hi", "n_pairs_ci_gt0"]]
    sup_rel = sup_rel.merge(ci, on="cell_type", how="inner", validate="1:1")

    # out-of-fold CelMod r, both levels (shared by panels a, b, c)
    S = pd.read_csv(R / "115_celmod_dfc_supertype.csv")
    B = pd.read_csv(R / "115_celmod_dfc_subclass.csv")
    S, B = add_extra_legs(S, B, tax)

    # a: per-supertype median out-of-fold r (across bulk x snRNA pairs) vs snRNA reliability
    bulk_r = S.groupby("cell_type").agg(bulk_rho=("rho", "median"), bulk_r=("r", "median"), bulk_n_pairs=("rho", "count")).reset_index()
    a = sup_rel.merge(bulk_r, on="cell_type", how="inner")
    a["compartment"] = a.cell_type.map(comp)
    a["reactive"] = a.cell_type.map(react)
    a["class_id"] = a.cell_type.map(class_of)
    a["core"] = a.cell_type.isin(core)
    a["concordant"] = a.ci_lo > ICC_THRESH
    assert a.compartment.notna().all()
    assert a.class_id.notna().all()
    a.to_csv(OUT / "fig4a_icc_vs_deconv.tsv", sep="\t", index=False)
    print(f"a: {len(a)} supertypes; {a.concordant.sum()} Core (CI lower bound > {ICC_THRESH}), "
          f"{sum(c in core for c in a[a.concordant].cell_type)} of those in CORE ({len(core)} total); "
          f"{(a.median_icc_sym >= 0.3).sum()} with median ICC >= 0.3")

    # b: out-of-fold r by group (ICC > 0 supertypes / ICC <= 0 supertypes / subclasses), per pair
    concordant = set(a[a.concordant].cell_type)
    S2 = S[S.cell_type.isin(a.cell_type)].copy()   # same supertypes as panel a (a type with no ICC is not called non-core)
    S2["set"] = [SET_HI if c in concordant else SET_LO for c in S2.cell_type]
    B2 = B.copy()
    B2["set"] = SET_SUB
    b = pd.concat([S2, B2])[["pair", "bulk", "snRNA", "level", "cell_type", "set", "n_donors", "r", "rho"]]
    b.to_csv(OUT / "fig4b_oof_by_concordance.tsv", sep="\t", index=False)
    print("b: median out-of-fold rho by pair x group:")
    print(b.groupby(["pair", "set"]).rho.median().unstack().round(3))
    pval_rows = []
    for pair, g in b.groupby("pair"):
        x = g[g.set == SET_HI].rho.dropna()
        for comp_name in (SET_LO, SET_SUB):
            y = g[g.set == comp_name].rho.dropna()
            if not len(y):   # NABEC x Catching has no subclass leg
                continue
            u = stats.mannwhitneyu(x, y, alternative="greater")
            pval_rows.append({"pair": pair, "comparison": comp_name,
                              "n_concordant": len(x), "n_comparison": len(y), "p": u.pvalue})
            print(f"  {pair}: concordant (n={len(x)}) > {comp_name.lower()} (n={len(y)}), "
                  f"Mann-Whitney p={u.pvalue:.2e}")
    pvals = pd.DataFrame(pval_rows)
    pvals.to_csv(OUT / "fig4b_pvals.tsv", sep="\t", index=False)

    # supplement: does the panel-a relationship hold per pair? Spearman rho(ICC, out-of-fold r)
    m = S.merge(sup_rel, on="cell_type")
    rows = []
    for pair, g in m.groupby("pair"):
        g2 = g.dropna(subset=["median_icc_sym", "rho"])
        rho, p = stats.spearmanr(g2.median_icc_sym, g2.rho)
        rows.append({"pair": pair, "n": len(g2), "rho": rho, "p": p})
    c = pd.DataFrame(rows).sort_values("rho")
    c.to_csv(OUT / "fig4_supp_rho_by_pair.tsv", sep="\t", index=False)
    print("supp: per-pair Spearman rho (ICC vs out-of-fold r):")
    print(c.to_string(index=False))
    print(f"wrote 4 tables to {OUT}")


if __name__ == "__main__":
    main()
