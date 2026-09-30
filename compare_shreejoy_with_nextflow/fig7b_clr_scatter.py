"""Donor-level within-compartment CLR for matched snRNA-seq pairs.

Writes a long table used by fig7b_clr_scatter.R. Does not edit the wp3 tree.
Uses the same caches and CLR definition as wp3lib.decompose (pseudocount 0.5,
CLR closed within neuronal vs non-neuronal).
"""
from pathlib import Path
import numpy as np
import pandas as pd

CACHE = Path("/scratch/shreejoy/cell_type_bias/wp3/cache")
OUT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow/data")
OUT.mkdir(parents=True, exist_ok=True)

PSEUDO = 0.5
NEURONAL_PREFIX = "Neuronal"

SHORTLIST = [
    "Astro_6-SEAAD", "OPC_2_2-SEAAD", "OPC_2",
    "Sst_25", "Sst_3", "Micro-PVM_2_1-SEAAD",
    "Sst_23", "Micro-PVM_2", "Sst_2",
]

PAIRS = [
    ("Green_person", "Mathys_person", "Green", "Mathys", "Green x Mathys"),
    ("Green_person", "Mome25_person", "Green", "Multiome 2025", "Green x Multiome 2025"),
    ("Mathys_person", "Mome25_person", "Mathys", "Multiome 2025", "Mathys x Multiome 2025"),
    ("Mome25_person", "PsychAD_BRISC_person", "Multiome 2025", "PsychAD", "Multiome 2025 x PsychAD"),
    ("BrainSCOPE_CMC", "PsychAD_BRISC_person", "BrainSCOPE CMC", "PsychAD", "BrainSCOPE CMC x PsychAD"),
    ("PsychAD_BRISC_person", "Ruzicka_MtSinai", "PsychAD", "Ruzicka", "PsychAD x Ruzicka"),
]


def load_cells(name):
    p = CACHE / f"{name}_cells.parquet"
    df = pd.read_parquet(p, columns=["donor", "supertype", "cls"])
    df["donor"] = df["donor"].astype(str)
    return df


def count_grid(cells):
    df = cells.copy()
    df["grp"] = np.where(df.cls.str.startswith(NEURONAL_PREFIX), "neu", "glia")
    t = pd.crosstab(df.donor, [df.grp, df.supertype])
    groups = {}
    for g in ("neu", "glia"):
        if g in t.columns.get_level_values(0):
            groups[g] = t[g]
    return groups


def clr(counts):
    p = (counts + PSEUDO)
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def pair_clr(gA, gB, types):
    shared = None
    rows = []
    for grp in ("neu", "glia"):
        if grp not in gA or grp not in gB:
            continue
        tA, tB = gA[grp], gB[grp]
        all_types = sorted(set(tA.columns) & set(tB.columns))
        keep_types = [c for c in types if c in all_types]
        if not keep_types:
            continue
        idx = sorted(set(tA.index.astype(str)) & set(tB.index.astype(str)))
        dA = tA.reindex(index=idx, columns=all_types).fillna(0).to_numpy()
        dB = tB.reindex(index=idx, columns=all_types).fillna(0).to_numpy()
        ok_donor = (dA.sum(axis=1) > 0) & (dB.sum(axis=1) > 0)
        dA, dB = dA[ok_donor], dB[ok_donor]
        donors = np.array(idx)[ok_donor]
        if shared is None:
            shared = len(donors)
        yA, yB = clr(dA), clr(dB)
        col = {ct: i for i, ct in enumerate(all_types)}
        for ct in keep_types:
            j = col[ct]
            a, b = yA[:, j], yB[:, j]
            finite = np.isfinite(a) & np.isfinite(b)
            for donor, va, vb in zip(donors[finite], a[finite], b[finite]):
                rows.append(dict(
                    cell_type=ct, group=grp, donor=donor,
                    clr_a=float(va), clr_b=float(vb),
                ))
    return rows, shared


def main():
    grids = {}
    out = []
    for cache_a, cache_b, lab_a, lab_b, pair in PAIRS:
        print(f"{pair} ...", flush=True)
        for name in (cache_a, cache_b):
            if name not in grids:
                grids[name] = count_grid(load_cells(name))
        rows, n = pair_clr(grids[cache_a], grids[cache_b], SHORTLIST)
        for r in rows:
            r.update(pair=pair, dataset_a=lab_a, dataset_b=lab_b, n_shared=n)
        out.extend(rows)
        have = sorted({r["cell_type"] for r in rows})
        print(f"  n_shared={n}  types={have}", flush=True)

    df = pd.DataFrame(out)
    path = OUT / "fig7b_clr_scatter.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(f"wrote {path}  rows={len(df)}", flush=True)


if __name__ == "__main__":
    main()
