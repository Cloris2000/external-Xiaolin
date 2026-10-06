"""Figure 8 inputs: does sequencing depth explain per-donor discordance?

Question (user, 2026-10-05): for the five independent matched-donor pairs used
in Figure 1 (common/reliability_stats.py), is a donor's cross-dataset CLR
discordance larger when the two datasets sequenced that donor to very
different depths? This is SEQUENCING depth (UMIs / genes detected per
nucleus), not the cortical SAMPLING depth already covered by Figure 3
(fig3_donor_depth.tsv, a laminar-position proxy built from the SCZ Xenium
prior -- a different QC axis that happens to share the word "depth").

Per-donor sequencing depth is not in the wp3 _cells.parquet caches (they only
carry donor/supertype/confidence/subclass/cls). Two sources:

  Green, Mathys: already extracted per donor in
    /scratch/shreejoy/cell_type_bias/wp3/results/165_sample_qc.csv
    (columns umi, genes = per-donor median; measurement Green_DFC/Mathys_DFC).

  Multiome2025_DLPFC, PsychAD_full, BrainSCOPE_CMC: not previously extracted.
  Pulled here directly from each dataset's raw h5ad obs (cell_id/index plus
  nCount_RNA/nFeature_RNA or n_counts/n_genes -- no count matrix is read, obs
  only), joined onto the BRISC label CSV's cell_id, donor recovered and
  re-keyed to registry person_id EXACTLY as 48_brisc_to_cache.py does (so the
  resulting donor ids match pair_stats_percelltype.tsv):
    - PsychAD has no donor column in its label file or h5ad; cell_id encodes
      it as "<donor>-<lane>-<barcode>-<suffix>" (48_brisc_to_cache.py's own
      recovery regex, checked there against the h5ad's own donor levels).
    - Native donor ids (individualID / Donor / recovered PsychAD donor) are
      re-keyed to registry person_id via 21_person_crosswalk.csv, same as the
      cache build, so donor ids here line up with common/reliability_stats.py.

Per-donor discordance (ALL cell types, not just inhibitory, per user
decision 2026-10-05): for each of the five pairs' shared donors, the mean
absolute difference between that donor's supertype-level CLR vectors in
dataset A and B, centred by removing each cell type's across-donor mean in
each dataset first (so a systematic offset between datasets, e.g. driven by
reference/annotation differences, is not counted as discordance). This uses
the SAME supertypes/donor filters as common/reliability_stats.py's pair_stats
(neu/glia compartments, MIN_DONORS floor skipped here since we report the
raw per-donor residual, not a summary stat that needs it).

Outputs (data/ next to this script):
  fig8_donor_depth.tsv   one row per (dataset, donor): median_umi, median_genes
  fig8_pair_donor.tsv    one row per (pair, donor): depth in A and B, depth
                          gap (log2 ratio), and mean |residual CLR| across all
                          shared supertypes

Run: must go through Slurm (PsychAD h5ad is 113 GB on disk; obs-only read is
light but the inventory docs warn against anything large on a login node --
see fig8_data.sbatch). /project is read-only on Trillium compute nodes, so
this writes to $FIG8_OUT (set by fig8_data.sbatch to a scratch dir) rather
than to data/ directly; the sbatch script copies the outputs back from a
login node after the job finishes.
  PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig8_data.py
"""
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

WP3 = Path("/scratch/shreejoy/cell_type_bias/wp3")
CACHE = WP3 / "cache"
LABELS = WP3 / "labels"
RESULTS = WP3 / "results"
OUT = Path(os.environ.get("FIG8_OUT", Path(__file__).resolve().parent / "data"))

# the five independent pairs (common/reliability_stats.py, 2026-10-05 set)
PAIRS = [
    ("Green_DFC_person", "Mathys_DFC_person", "Green x Mathys"),
    ("Green_DFC_person", "Multiome2025_DLPFC_DFC_person", "Green x Multiome"),
    ("Mathys_DFC_person", "Multiome2025_DLPFC_DFC_person", "Mathys x Multiome"),
    ("Multiome2025_DLPFC_DFC_person", "PsychAD_full_DFC_person", "Multiome x PsychAD"),
    ("BrainSCOPE_CMC_DFC_person", "PsychAD_full_DFC_person", "BrainSCOPE CMC x PsychAD"),
]
# dataset name (as used in PAIRS/cache) -> depth extraction recipe
RAW_QC = {
    "Multiome2025_DLPFC_DFC_person": dict(
        h5=Path("/scratch/shreejoy/AMPAD_multiome2025/AMPAD_multiome2025_DLPFC.h5ad"),
        idx="_index", umi="nCount_RNA", genes="nFeature_RNA",
        labels=LABELS / "Multiome2025_DLPFC_DFC_brisc_labels.csv",
        donor_col="individualID", crosswalk_dataset="AMP-AD Diverse multiome",
    ),
    "PsychAD_full_DFC_person": dict(
        h5=Path("/scratch/medney/Datasets/PsychAD_Data/merged_final_clean.h5ad"),
        idx="barcodekey", umi="n_counts", genes="n_genes",
        labels=LABELS / "PsychAD_full_DFC_brisc_labels.csv",
        donor_col=None,  # recovered from cell_id, as 48_brisc_to_cache.py
        crosswalk_dataset=None,  # spans PsychAD MSBB/HBCC/ROSMAP; no single name
    ),
    "BrainSCOPE_CMC_DFC_person": dict(
        h5=Path("/scratch/medney/Datasets/Brain_scope/CMC.h5ad"),
        idx="cell", umi="nCount_RNA", genes="nFeature_RNA",
        labels=LABELS / "BrainSCOPE_CMC_DFC_brisc_labels.csv",
        donor_col="Donor", crosswalk_dataset="BrainSCOPE CMC",
    ),
}
PSEUDO = 0.5
MIN_CELLS_DONOR = 1   # report every shared donor; no composition floor here


def short(name):
    return name.replace("_DFC_person", "").replace("_full", "")


def read_obs_col(h, name):
    """obs column, numeric or string, handling anndata's categorical and
    nullable-string (values/mask) group encodings."""
    g = h["obs"][name]
    if isinstance(g, h5py.Group):
        if "categories" in g:
            cats = [x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]]
            return pd.Series(g["codes"][:]).map(dict(enumerate(cats)))
        if "values" in g:
            arr = g["values"][:]
            return pd.Series([x.decode() if isinstance(x, bytes) else str(x) for x in arr])
        raise TypeError(f"unrecognized obs group layout for {name}: {list(g.keys())}")
    arr = g[:]
    if arr.dtype.kind in "SO":
        return pd.Series([x.decode() if isinstance(x, bytes) else str(x) for x in arr])
    return pd.Series(arr)


def extract_raw_depth(name, spec, crosswalk):
    """Per-donor median UMI/genes for one dataset not already in 165_sample_qc.csv."""
    print(f"  {name}: reading obs from {spec['h5']}", flush=True)
    with h5py.File(spec["h5"], "r") as h:
        cell_id = read_obs_col(h, spec["idx"]).astype(str)
        umi = pd.to_numeric(read_obs_col(h, spec["umi"]), errors="coerce")
        genes = pd.to_numeric(read_obs_col(h, spec["genes"]), errors="coerce")
    obs = pd.DataFrame({"cell_id": cell_id.values, "umi": umi.values, "genes": genes.values})

    lab = pd.read_csv(spec["labels"], usecols=lambda c: c in
                       {"cell_id", "passed_QC", "brisc_Supertype", spec["donor_col"]})
    lab = lab[(lab.passed_QC == True) & lab.brisc_Supertype.notna()]
    if spec["donor_col"] is None:
        donor = lab.cell_id.astype(str).str.replace(r"-\d+-[ACGTN]+-\d+$", "", regex=True)
    else:
        donor = lab[spec["donor_col"]].astype(str)
    lab = pd.DataFrame({"cell_id": lab.cell_id.astype(str), "native_donor": donor})

    m = lab.merge(obs, on="cell_id", how="left")
    miss = m.umi.isna().mean()
    if miss > 0.01:
        raise SystemExit(f"{name}: {miss:.1%} of QC-pass cells had no raw QC match "
                         f"by cell_id -- check the id convention")

    if spec["crosswalk_dataset"] is not None:
        xw = crosswalk[crosswalk.dataset == spec["crosswalk_dataset"]]
        n2p = dict(zip(xw.native_donor.astype(str), xw.person_id))
        donor_id = m.native_donor.map(n2p)
        unmapped = donor_id.isna().mean()
        if unmapped > 0.01:
            raise SystemExit(f"{name}: {unmapped:.1%} of cells had no registry "
                             f"crosswalk for dataset '{spec['crosswalk_dataset']}'")
        m["donor"] = donor_id
    else:
        # PsychAD: re-key through whichever of its three crosswalk rows matches
        xw = crosswalk[crosswalk.dataset.isin(
            ["PsychAD MSBB", "PsychAD HBCC", "PsychAD ROSMAP"])]
        n2p = dict(zip(xw.native_donor.astype(str), xw.person_id))
        donor_id = m.native_donor.map(n2p)
        unmapped = donor_id.isna().mean()
        if unmapped > 0.01:
            raise SystemExit(f"{name}: {unmapped:.1%} of PsychAD cells had no "
                             f"registry crosswalk")
        m["donor"] = donor_id

    g = m.dropna(subset=["donor"]).groupby("donor")
    d = pd.DataFrame({"dataset": short(name), "donor": g.size().index,
                      "n_nuclei_raw": g.size().values,
                      "median_umi": g.umi.median().values,
                      "median_genes": g.genes.median().values})
    print(f"    {len(d)} donors, median UMI {d.median_umi.median():.0f}, "
          f"median genes {d.median_genes.median():.0f}")
    return d


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def load_cells(name):
    c = pd.read_parquet(CACHE / f"{name}_cells.parquet",
                        columns=["donor", "supertype", "cls"])
    c["donor"] = c["donor"].astype(str)
    return c


def pair_donor_discordance(ca, cb):
    """Per-donor mean |residual CLR| across ALL shared supertypes (both neu and
    glia compartments pooled), residual = CLR minus that cell type's own
    across-donor mean in its own dataset, so a systematic between-dataset
    offset in a cell type is not counted as a donor's discordance."""
    grpA = np.where(ca.cls.str.startswith("Neuronal"), "neu", "glia")
    grpB = np.where(cb.cls.str.startswith("Neuronal"), "neu", "glia")
    tA = pd.crosstab(ca.donor, [grpA, ca.supertype])
    tB = pd.crosstab(cb.donor, [grpB, cb.supertype])
    shared_donors = sorted(set(ca.donor) & set(cb.donor))
    resid_sum = pd.Series(0.0, index=shared_donors)
    n_types = 0
    for grp in sorted(set(tA.columns.get_level_values(0)) & set(tB.columns.get_level_values(0))):
        types = [t for t in tA[grp].columns if t in tB[grp].columns]
        if not types:
            continue
        dA = tA[grp].reindex(index=shared_donors, columns=types).fillna(0).values
        dB = tB[grp].reindex(index=shared_donors, columns=types).fillna(0).values
        keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
        if keep.sum() < 2:
            continue
        yA, yB = clr(dA[keep]), clr(dB[keep])
        resA = yA - yA.mean(axis=0, keepdims=True)
        resB = yB - yB.mean(axis=0, keepdims=True)
        absdiff = np.abs(resA - resB)
        ok = np.isfinite(absdiff)
        donors_keep = np.array(shared_donors)[keep]
        for j in range(len(types)):
            col_ok = ok[:, j]
            resid_sum.loc[donors_keep[col_ok]] += absdiff[col_ok, j]
        n_types += len(types)
    return pd.DataFrame({"donor": shared_donors,
                         "mean_abs_resid_clr": (resid_sum / n_types).values,
                         "n_types": n_types})


def main():
    OUT.mkdir(exist_ok=True)
    crosswalk = pd.read_csv(RESULTS / "21_person_crosswalk.csv").dropna(subset=["person_id"])

    qc = pd.read_csv(RESULTS / "165_sample_qc.csv")
    have = qc[qc.measurement.isin(["Green_DFC", "Mathys_DFC"]) & qc.umi.notna()]
    depth_rows = [pd.DataFrame({
        "dataset": have.measurement.str.replace("_DFC", ""), "donor": have.person_id,
        "n_nuclei_raw": have.n_nuclei, "median_umi": have.umi, "median_genes": have.genes,
    })]
    print(f"Green/Mathys: {len(depth_rows[0])} donors from 165_sample_qc.csv")

    for name, spec in RAW_QC.items():
        depth_rows.append(extract_raw_depth(name, spec, crosswalk))

    depth = pd.concat(depth_rows, ignore_index=True)
    depth.to_csv(OUT / "fig8_donor_depth.tsv", sep="\t", index=False)

    Dk = depth.set_index(["dataset", "donor"])
    pair_rows = []
    for a, b, label in PAIRS:
        ca, cb = load_cells(a), load_cells(b)
        disc = pair_donor_discordance(ca, cb)
        g = lambda ds, col: Dk[col].reindex(
            list(zip([short(ds)] * len(disc), disc.donor))).values
        df = pd.DataFrame({
            "pair": label, "donor": disc.donor,
            "mean_abs_resid_clr": disc.mean_abs_resid_clr, "n_types": disc.n_types,
            "umi_a": g(a, "median_umi"), "umi_b": g(b, "median_umi"),
            "genes_a": g(a, "median_genes"), "genes_b": g(b, "median_genes"),
        })
        df = df.dropna(subset=["umi_a", "umi_b"])
        df["log2_umi_ratio"] = np.log2(df.umi_b / df.umi_a)
        df["log2_genes_ratio"] = np.log2(df.genes_b / df.genes_a)
        pair_rows.append(df)
        print(f"  {label:26s} {len(df)} shared donors with depth on both sides "
              f"(of {len(disc)} shared donors total)")
    pd.concat(pair_rows, ignore_index=True).to_csv(OUT / "fig8_pair_donor.tsv", sep="\t", index=False)
    print(f"wrote {OUT}/fig8_donor_depth.tsv, fig8_pair_donor.tsv")


if __name__ == "__main__":
    main()
