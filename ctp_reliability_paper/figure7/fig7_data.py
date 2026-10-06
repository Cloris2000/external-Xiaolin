"""Data for Figure 7: cross-region reliability of cell-type proportions.

Companion to Figure 1 (same-region, cross-dataset reliability). Figure 7 asks
the complementary question: for the SAME donors and the SAME dataset, how
reliable is a supertype's proportion when measured in a second cortical
region instead of a second dataset? Same three lineages as Figure 1 (Immune,
Sst, L2/3 IT), same metric (symmetric ICC on donor-level CLR), same pseudocount.

Two cohorts, each a single matched-donor region pair (not five independent
pairs as in Figure 1), so there is no across-pair IQR to report. Instead each
point's interval is a donor-bootstrap 25th-75th percentile (2026-10 decision):

  (1) Multiome2025 (AMP-AD Diverse multiome, 2025 release): DLPFC vs STG,
      same donors. wp3's own cache only has DLPFC (46_multiome2025_convert.R
      deliberately drops STG/caudate). STG is only available, already
      BRISC-labelled with the same DFC-2026 taxonomy (verified identical
      supertype/subclass/class sets against wp3's own Multiome2025_DLPFC_DFC
      cache), in a separate atlas pipeline's cell tables:
        /scratch/shreejoy/composition/results_2026-09-27b/atlas_umap/full/cells/
          Multiome2025_DLPFC.parquet, Multiome2025_STG.parquet
      Those tables key donors by a "BDR..." person_id that is NOT a wp3/
      registry id. Bridged here via: atlas cell_id (batch_barcode) -> raw
      metadata cell key -> individualID -> wp3's own person crosswalk
      (21_person_crosswalk.csv, "AMP-AD Diverse multiome" is REGISTRY_DIRECT)
      -> TLP person_id. Verified: 0 unmatched/unmapped cells either region,
      412/419 TLP donors (DLPFC/STG), 398 shared.

  (2) SEA-AD Multiregion 2026: DFC vs MTG, same donors. Already cached by
      174_seaad_multiregion.py as donor x region x supertype counts, native
      SEA-AD supertype labels that are already DFC-2026-taxonomy-compatible
      (verified: all 35 supertypes under Immune/Sst/L2/3 IT present). 80
      shared donors.

Outputs (data/ next to this script):
  fig7a_pairs.tsv     per cohort x level x cell_type, the three target
                      lineages only (Immune, Sst, L2/3 IT) plus their
                      supertypes: n_donors, icc_sym, and bootstrap
                      q25/q75/median (B=2000 resamples of donors)
  fig7_all_supertypes.tsv   same icc_sym (no bootstrap, to keep runtime
                      down), supertype level, ALL ~150 DFC-2026 supertypes
                      measurable in both cohorts -- for the all-cell-type
                      panel (b)
  fig7bc_clr.tsv    donor-level CLR for the subclasses + their supertypes,
                    both regions, both cohorts (for scatter panels)

Run (login node is fine for the original lineage-only pass; reads cached
parquet/csv only, no raw processing):
  PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig7_data.py

The all-supertypes pass (fig7_all_supertypes.tsv, ~150 types x 2 cohorts)
exceeds the login node's CPU-time limit, and /project is read-only on
Trillium compute nodes (see nextflow/site_env.sh), so under Slurm this must
write to scratch and be copied back from a login node afterward:
  sbatch fig7_data.sbatch
  cp /scratch/zhoux156/ctp_reliability_paper/figure7_data/*.tsv data/
FIG7_OUT overrides the output directory (default: data/ next to this script).
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

WP3_CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
WP3_RESULTS = Path("/scratch/shreejoy/cell_type_bias/wp3/results")
ATLAS_CELLS = Path("/scratch/shreejoy/composition/results_2026-09-27b/atlas_umap/full/cells")
RAW_MULTIOME_META = Path("/project/rrg-shreejoy/AMP_AD_Diverse_snRNA/Metadata/"
                         "rosmap_diversity_snmultiome_clusterdata_20241107.csv")
SEAAD_MULTIREGION_CACHE = Path("/scratch/shreejoy/seaad_multiregion/"
                               "counts_donor_region_supertype.parquet")
TAX = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/"
           "refs/taxonomy_DFC_2026.tsv")
OUT = Path(os.environ.get("FIG7_OUT", Path(__file__).resolve().parent / "data"))

SUBCLASSES = ["Immune", "Sst", "L2/3 IT"]
PSEUDO = 0.5
NEURONAL_PREFIX = "Neuronal"
MIN_DONORS = 40    # same floor as common/reliability_stats.py
N_BOOT = 2000
SEED = 0

COHORTS = [
    ("Multiome2025", "DLPFC", "STG"),
    ("SEA-AD Multiregion", "DFC", "MTG"),
]


# ---- shared CLR / ICC machinery (identical to common/reliability_stats.py) --
def count_grid(cells, level):
    grp = np.where(cells.cls.str.startswith(NEURONAL_PREFIX), "neu", "glia")
    t = pd.crosstab(cells.donor, [grp, cells[level]])
    return {g: t[g] for g in ("neu", "glia") if g in t.columns.get_level_values(0)}


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def icc_sym_from_clr(a, b):
    ac, bc = a - a.mean(), b - b.mean()
    cov = float(np.cov(ac, bc, ddof=1)[0, 1])
    va, vb = float(ac.var(ddof=1)), float(bc.var(ddof=1))
    return 2 * cov / (va + vb)


def region_stats(ca, cb, level, cohort, rng, keep_types=None, bootstrap=True):
    """Same-donor, cross-region icc_sym per cell type, with a donor-bootstrap
    25th/50th/75th percentile in place of Figure 1's across-pair IQR.

    ca/cb must carry the FULL cell population (all subclasses/supertypes),
    not just the reported ones: CLR is compositional, so restricting the
    input to a handful of types first collapses a compartment to one column
    and degenerates its CLR to zero variance (e.g. "Immune" is the only
    reported glia-group subclass, so filtering first makes its CLR constant).
    keep_types subsets which cell types are reported, after CLR is computed
    over the full compartment. bootstrap=False skips the per-type resample
    (point estimate only), for when every supertype is reported (panel b's
    ~150 types x 2 cohorts) and the per-pair whisker is not plotted anyway.
    """
    shared = sorted(set(ca.donor) & set(cb.donor))
    gA, gB = count_grid(ca, level), count_grid(cb, level)
    rows = []
    for grp in sorted(set(gA) & set(gB)):
        types = [c for c in gA[grp].columns if c in gB[grp].columns]
        dA = gA[grp].reindex(index=shared, columns=types).fillna(0).values
        dB = gB[grp].reindex(index=shared, columns=types).fillna(0).values
        keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
        dA, dB = dA[keep], dB[keep]
        donors = np.array(shared)[keep]
        yA, yB = clr(dA), clr(dB)
        n = len(donors)
        if n < MIN_DONORS:
            continue
        for j, ct in enumerate(types):
            if keep_types is not None and ct not in keep_types:
                continue
            a, b = yA[:, j], yB[:, j]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < MIN_DONORS:
                continue
            a_ok, b_ok = a[ok], b[ok]
            icc = icc_sym_from_clr(a_ok, b_ok)
            boot_median = q25 = q75 = np.nan
            if bootstrap:
                ok_idx = np.flatnonzero(ok)
                resample_pos = rng.integers(0, len(ok_idx), size=(N_BOOT, len(ok_idx)))
                boot = np.array([icc_sym_from_clr(a_ok[p], b_ok[p]) for p in resample_pos])
                boot = boot[np.isfinite(boot)]
                if len(boot):
                    boot_median, q25, q75 = (float(np.median(boot)),
                                             float(np.quantile(boot, 0.25)),
                                             float(np.quantile(boot, 0.75)))
            rows.append(dict(
                cohort=cohort, level=level, cell_type=ct, group=grp,
                n_donors=int(ok.sum()), icc_sym=icc,
                boot_median=boot_median, q25=q25, q75=q75,
            ))
    return pd.DataFrame(rows)


# ---- cohort 1: Multiome2025 DLPFC vs STG, bridged to wp3 TLP person ids -----
def load_multiome_region(region_label, raw_region_name, crosswalk):
    cells = pd.read_parquet(ATLAS_CELLS / f"Multiome2025_{region_label}.parquet",
                            columns=["cell_id", "supertype", "subclass", "class"])
    md = pd.read_csv(RAW_MULTIOME_META, usecols=["individualID", "region", "batch", "barcode"])
    md = md[md.region == raw_region_name].copy()
    md["cell_id"] = md.batch + "_" + md.barcode
    m = cells.merge(md[["cell_id", "individualID"]], on="cell_id", how="left")
    if m.individualID.isna().any():
        raise SystemExit(f"Multiome2025 {region_label}: {m.individualID.isna().sum()} "
                         f"cells did not match the raw metadata")
    m = m.merge(crosswalk, on="individualID", how="left")
    if m.TLP.isna().any():
        raise SystemExit(f"Multiome2025 {region_label}: {m.TLP.isna().sum()} "
                         f"cells had no registry crosswalk")
    m = m.rename(columns={"TLP": "donor", "class": "cls"})
    m["donor"] = m["donor"].astype(str)
    print(f"  Multiome2025 {region_label}: {m.donor.nunique()} TLP donors, {len(m):,} nuclei")
    return m[["donor", "supertype", "subclass", "cls"]]


def build_multiome_cohort():
    cw = pd.read_csv(WP3_RESULTS / "21_person_crosswalk.csv")
    cw = cw[cw.dataset == "AMP-AD Diverse multiome"][["native_donor", "person_id"]].rename(
        columns={"native_donor": "individualID", "person_id": "TLP"})
    dlpfc = load_multiome_region("DLPFC", "dorsolateral prefrontal cortex", cw)
    stg = load_multiome_region("STG", "superior temporal gyrus", cw)
    return dlpfc, stg


# ---- cohort 2: SEA-AD Multiregion DFC vs MTG, already cached --------------
def build_seaad_cohort(st2sc_full):
    """taxonomy_DFC_2026.tsv is DFC-specific and does not list every supertype
    that appears across SEA-AD Multiregion's ten regions (hippocampal
    CA/DG/Sub, entorhinal EC IT, and a handful of other-region-only states,
    e.g. L5 ET_1 -- checked 2026-10-05: none belong to Immune/Sst/L2/3 IT).
    For the subclass-level CLR these still need to occupy a column in their
    neu/glia compartment (dropping them instead of bucketing them would
    silently shrink the compositional denominator and inflate the remaining
    columns -- see region_stats docstring), so anything outside the target
    taxonomy is pooled into one "Other" column per compartment. The
    supertype-level CLR does not need a subclass label at all, so this only
    affects the subclass grouping.
    """
    cnt = pd.read_parquet(SEAAD_MULTIREGION_CACHE)
    cnt["donor"] = cnt["donor"].astype(str)

    def to_cells(region):
        sub = cnt[cnt.region == region]
        rows = sub.loc[sub.index.repeat(sub.n)]
        print(f"  SEA-AD Multiregion {region}: {sub.donor.nunique()} donors, "
              f"{int(sub.n.sum()):,} nuclei")
        subclass = rows.Supertype.map(st2sc_full).fillna("Other")
        return pd.DataFrame(dict(
            donor=rows.donor.values, supertype=rows.Supertype.values,
            subclass=subclass.values, cls=rows.cls.values))

    dfc = to_cells("DFC")
    mtg = to_cells("MTG")
    return dfc, mtg


def main():
    OUT.mkdir(exist_ok=True)
    global st2sc
    tax_full = pd.read_csv(TAX, sep="\t")
    st2sc_full = dict(zip(tax_full.supertype_label, tax_full.subclass_label))
    tax = tax_full[tax_full.subclass_label.isin(SUBCLASSES)]
    supertypes = set(tax.supertype_label)
    st2sc = dict(zip(tax.supertype_label, tax.subclass_label))   # target-only, for reporting
    keep_types = set(SUBCLASSES) | supertypes

    print("Cohort 1: Multiome2025 DLPFC vs STG")
    mome_a, mome_b = build_multiome_cohort()
    print("Cohort 2: SEA-AD Multiregion DFC vs MTG")
    seaad_a, seaad_b = build_seaad_cohort(st2sc_full)

    cohort_data = {
        "Multiome2025": (mome_a, mome_b),
        "SEA-AD Multiregion": (seaad_a, seaad_b),
    }

    rng = np.random.default_rng(SEED)
    stats_out, clr_out, all_st_out = [], [], []
    for cohort, region_a, region_b in COHORTS:
        ca_full, cb_full = cohort_data[cohort]
        ca = ca_full[ca_full.subclass.isin(SUBCLASSES)]
        cb = cb_full[cb_full.subclass.isin(SUBCLASSES)]
        shared = sorted(set(ca.donor) & set(cb.donor))
        print(f"{cohort} ({region_a} vs {region_b}): shared donors {len(shared)}")
        # CLR must be computed over the FULL population (all subclasses), not
        # just the reported three -- see region_stats docstring. keep_types
        # restricts which rows are reported.
        for level, lvl_keep in (("subclass", set(SUBCLASSES)), ("supertype", supertypes)):
            st = region_stats(ca_full, cb_full, level, cohort, rng, keep_types=lvl_keep)
            stats_out.append(st)

        # panel b: every supertype measurable in this cohort, point estimate
        # only (no bootstrap -- not plotted, and ~150 types x 2 cohorts would
        # be slow for no benefit)
        all_st = region_stats(ca_full, cb_full, "supertype", cohort, rng,
                              keep_types=None, bootstrap=False)
        all_st_out.append(all_st)

        # donor-level CLR for the subclasses + all their supertypes, for
        # scatter panels (analogous to fig1bc_clr.tsv). Same full-population
        # requirement as region_stats: CLR over the full compartment, keep_types
        # only subsets which columns are written out.
        gA, gB = count_grid(ca_full, "supertype"), count_grid(cb_full, "supertype")
        for grp in sorted(set(gA) & set(gB)):
            types = [c for c in gA[grp].columns if c in gB[grp].columns and c in keep_types]
            if not types:
                continue
            dA = gA[grp].reindex(index=shared, columns=[c for c in gA[grp].columns
                                                          if c in gB[grp].columns]).fillna(0).values
            dB = gB[grp].reindex(index=shared, columns=[c for c in gA[grp].columns
                                                          if c in gB[grp].columns]).fillna(0).values
            keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
            donors = np.array(shared)[keep]
            all_types = [c for c in gA[grp].columns if c in gB[grp].columns]
            yA, yB = clr(dA[keep]), clr(dB[keep])
            for j, ct in enumerate(all_types):
                if ct not in types:
                    continue
                a, b = yA[:, j], yB[:, j]
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() < MIN_DONORS:
                    continue
                clr_out.append(pd.DataFrame(dict(
                    cohort=cohort, region_a=region_a, region_b=region_b,
                    cell_type=ct, group=grp, donor=donors[ok],
                    clr_a=a[ok], clr_b=b[ok])))

    stats = pd.concat(stats_out, ignore_index=True)
    stats["subclass"] = np.where(stats.level == "subclass", stats.cell_type,
                                 stats.cell_type.map(st2sc))
    stats = stats[["cohort", "level", "subclass", "cell_type", "group", "n_donors",
                   "icc_sym", "boot_median", "q25", "q75"]]
    stats.to_csv(OUT / "fig7a_pairs.tsv", sep="\t", index=False)

    clr_df = pd.concat(clr_out, ignore_index=True)
    clr_df.to_csv(OUT / "fig7bc_clr.tsv", sep="\t", index=False)

    all_st = pd.concat(all_st_out, ignore_index=True)
    all_st["subclass"] = all_st.cell_type.map(st2sc_full).fillna("Other")
    all_st = all_st[["cohort", "subclass", "cell_type", "group", "n_donors", "icc_sym"]]
    all_st.to_csv(OUT / "fig7_all_supertypes.tsv", sep="\t", index=False)

    print(f"\nwrote {OUT}/fig7a_pairs.tsv ({len(stats)} rows), "
          f"fig7bc_clr.tsv ({len(clr_df)} rows), "
          f"fig7_all_supertypes.tsv ({len(all_st)} rows)")
    print(stats.groupby(["cohort", "level"]).cell_type.nunique().to_string())
    print(all_st.groupby("cohort").cell_type.nunique().to_string())


if __name__ == "__main__":
    main()
