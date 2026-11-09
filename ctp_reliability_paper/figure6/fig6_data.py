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

190 and 196 are saved only as medians over the PI's CORE set (185, partly selected on
cross-cohort bulk portability) and over all supertypes. The paper groups by the Fig 4
Core Supertype set instead (reliability-only; user decision 2026-10-28), so the
per-supertype values behind both are recomputed here by calling the PI's own harness
functions (187_harness.disattenuated / triangulate) on the PI's saved deconvolution
outputs, read-only, with the PI's arms, partners and denominators. Each recomputation
is checked against the PI's published CORE/all medians before it is written. The
triangulation bootstrap only feeds `portability`, which is not used, so N_BOOT is set
small; `reliability` is the full-sample estimate either way.

Writes data/fig6_*.tsv.
"""
import glob
import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats

R = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
D = Path("/scratch/shreejoy/cell_type_bias/wp3/deconv")
OUT = Path(__file__).resolve().parent / "data"
WP3 = Path("/scratch/shreejoy/cell_type_bias/wp3")
# Fig 4 Core Supertype: lower bound of the donor-bootstrap 95% CI of the median symmetric ICC > 0
RELIABILITY_CI = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/common/data/reliability_ci_summary.tsv")
ICC_THRESH = 0.0
# 190_score_arms.py: family -> held-out snRNA dataset; 196_validate_production.py: MSBB harness
FAMILIES = {"exclude_Green": "Green_DFC", "exclude_Mathys": "Mathys_DFC", "exclude_Multiome": "Multiome_DFC"}
REL_COL = {"Green_DFC": "rel_green", "Mathys_DFC": "rel_mathys"}
REGIONS = ["MSBB_FP", "MSBB_IFG", "MSBB_STG"]
PARTNER_TRAINER = "PsychAD_HBCC_DFC"
HARNESS_ARMS = (("production consensus", "prod__consensus"),
                ("solo Green (reference)", "xfer__Green_DFC"),
                ("solo Mathys (reference)", "xfer__Mathys_DFC"))
GROUPS = {  # gene-space groups of 195_production_model.py
    "rnaseq": ["ROSMAP", "ROSMAP_PCC", "MSBB_FP", "MSBB_IFG", "MSBB_STG", "NIMH_HBCC", "CommonMind", "GVEX"],
    "microarray": ["HBTRC"],
    "gtex_regions": ["GTEx_BA9_frz", "GTEx_BA24", "GTEx_CTX", "GTEx_CBL"],
    "wp1_bulk": ["Diverse_MSSM_DLPFC", "Diverse_MSSM_STG", "Diverse_MayoEmory_DLPFC",
                 "Diverse_MayoEmory_STG", "Diverse_Rush2_DLPFC", "Diverse_Rush2_STG",
                 "Diverse_Columbia_DLPFC", "Mayo_TCX", "NABEC"],
}


def harness():
    """The PI's 187_harness, imported from (and run in) its own directory."""
    os.environ["N_BOOT"] = "20"   # bootstrap feeds only `portability`, unused here
    os.chdir(WP3)
    sys.path.insert(0, str(WP3))
    return importlib.import_module("187_harness")


def per_type_arms(H):
    """190 per supertype: every arm's disattenuated r against its family's held-out dataset."""
    core = H.core_set()
    rel_mult = pd.read_csv(WP3 / "results" / "190_rel_multiome.csv", index_col=0).iloc[:, 0].to_dict()
    out = []
    for fam, held in FAMILIES.items():
        partner = H.snrna(held)
        rel = core[REL_COL[held]].to_dict() if held in REL_COL else rel_mult
        for f in sorted(glob.glob(str(WP3 / f"deconv/ROSMAP_DFC_supertype_arm__{fam}__*.parquet"))):
            arm = os.path.basename(f).split("arm__")[1].replace(".parquet", "")
            if "consTau" in arm:   # production uses the plain mean
                continue
            X = pd.read_parquet(f); X.index = X.index.astype(str)
            d = H.disattenuated(X, partner, rel)
            d["family"], d["arm"] = fam, arm
            out.append(d)
    return pd.concat(out, ignore_index=True)


def per_type_harness(H):
    """196 per supertype: reliability in the MSBB triads, median over region pairs."""
    snm = H.snrna("PsychAD_MSBB_DFC")
    out = []
    for name, pat in HARNESS_ARMS:
        vals = []
        for r in REGIONS:
            X = pd.read_parquet(WP3 / f"deconv/{r}_DFC_supertype_{pat}.parquet"); X.index = X.index.astype(str)
            for r2 in REGIONS:
                if r2 == r:
                    continue
                P2 = pd.read_parquet(WP3 / f"deconv/{r2}_DFC_supertype_xfer__{PARTNER_TRAINER}.parquet")
                P2.index = P2.index.astype(str)
                D = H.triangulate(X, snm, P2)
                if D is not None:
                    vals.append(D.set_index("cell_type"))
        M = pd.concat(vals).groupby(level=0).agg(reliability=("reliability", "median"),
                                                 sign_ok=("sign_ok", "mean"), n=("n", "median"))
        M["arm"] = name
        out.append(M.reset_index())
    return pd.concat(out, ignore_index=True)


def check(per, pub, key, value, pub_value, is_core):
    """Recomputed medians over the PI's CORE / all must equal the PI's published ones."""
    for (k, scope), want in pub.set_index([key, "scope"])[pub_value].items():
        d = per[per[key] == k]
        if scope == "CORE":
            d = d[d.cell_type.map(lambda c: is_core.get(c, False))]
        got = d[value].median()
        assert abs(got - want) < 1e-6, f"{k} {scope}: recomputed {got:.6f} vs published {want:.6f}"
    print(f"  {key}: recomputed medians match the PI's for {len(pub)} rows")


def paired(per, key, value, ref, others, keep):
    """Paired Wilcoxon per supertype, ref arm vs each other arm, on the supertypes in `keep`."""
    rows = []
    for o in others:
        j = per[per[key] == ref].set_index("cell_type")[[value]].join(
            per[per[key] == o].set_index("cell_type")[[value]], lsuffix="_ref", rsuffix="_o").dropna()
        j = j[j.index.isin(keep)]
        w = stats.wilcoxon(j[f"{value}_ref"], j[f"{value}_o"])
        rows.append(dict(ref=ref, other=o, n_types=len(j),
                         paired_diff=float((j[f"{value}_ref"] - j[f"{value}_o"]).median()), p=w.pvalue))
    return rows


def main():
    OUT.mkdir(exist_ok=True)
    ci = pd.read_csv(RELIABILITY_CI, sep="\t").query("level == 'supertype'")
    core_sup = dict(zip(ci.cell_type, ci.ci_lo > ICC_THRESH))
    keep = {c for c, v in core_sup.items() if v}
    print(f"Core Supertypes (Fig 4): {len(keep)} of {len(core_sup)}")
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
    x["core_supertype"] = x.cell_type.map(core_sup)
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

    # ---- per-supertype recomputation of 190 and 196 (PI harness, read-only) ----
    here = os.getcwd()
    H = harness()
    is_core = H.core_set().CORE.to_dict()
    pb = per_type_arms(H)
    check(pb, g[g.vs.isna()], "arm", "r_disatt", "median_r_disatt", is_core)
    pc = per_type_harness(H)
    check(pc, c, "arm", "reliability", "median_reliability", is_core)
    os.chdir(here)
    for d in (pb, pc):
        d["core_supertype"] = d.cell_type.map(core_sup)
    pb.to_csv(OUT / "fig6b_per_type.tsv", sep="\t", index=False)
    pc.to_csv(OUT / "fig6c_per_type.tsv", sep="\t", index=False)

    # paired tests on the Core Supertype set: consensus vs each single label; production vs each reference
    tests = []
    for fam in FAMILIES:
        arms = sorted(pb[pb.family == fam].arm.unique())
        cons = [a for a in arms if "__cons_" in a]
        assert len(cons) == 1, cons
        for r_ in paired(pb, "arm", "r_disatt", cons[0], [a for a in arms if "__solo_" in a], keep):
            tests.append(dict(panel="b", family=fam, **r_))
    for r_ in paired(pc, "arm", "reliability", HARNESS_ARMS[0][0], [a for a, _ in HARNESS_ARMS[1:]], keep):
        tests.append(dict(panel="c", family="MSBB harness", **r_))
    tests = pd.DataFrame(tests)
    tests.to_csv(OUT / "fig6_core_supertype_tests.tsv", sep="\t", index=False)
    print(tests.round(4).to_string())
    print(f"wrote 8 tables to {OUT}")


if __name__ == "__main__":
    main()
