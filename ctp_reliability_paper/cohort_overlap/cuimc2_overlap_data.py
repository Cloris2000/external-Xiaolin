"""Which other ROSMAP snRNA-seq datasets hold the CUIMC2 donors? (UpSet input)

CUIMC2 (Columbia 10x multiome of ROSMAP DLPFC, the follow-on to Green) enters
the paper through the 2026-10-28 overlap handoff (pairs CUIMC2 x Multiome2025,
CUIMC2 x Mathys). This script tabulates, for every CUIMC2 person, which of the
other ROSMAP-donor snRNA-seq datasets also measured that person.

Inputs (read-only): cleaned snRNA objects, vintage 2026-10-13
  /project/rrg-shreejoy/brain-data-registry/snrnaseq/2026-10-13/<dataset>.obs.parquet
A person counts as present in a dataset if it has >= 1 nucleus with
passed_QC == True (the standard keep rule; CUIMC2 demuxlet doublets fail it).
obs.person_id is registry release 2026-10-21c compatible (handoff, checked for
every donor), so datasets are joined on it directly. PsychAD is one object;
its ROSMAP arm is obs.dataset_id == PsychAD_ROSMAP.

Outputs (data/ next to this script):
  cuimc2_membership.tsv   one row per CUIMC2 person: n CUIMC2 nuclei, then 0/1
                          per other dataset
  cuimc2_overlap_sizes.tsv  per other dataset: persons shared with CUIMC2
Run (login node is fine, ~1 min, PsychAD obs is the large read):
  PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python -I cuimc2_overlap_data.py
"""
from pathlib import Path

import pandas as pd

OBJ = Path("/project/rrg-shreejoy/brain-data-registry/snrnaseq/2026-10-13")
OUT = Path(__file__).resolve().parent / "data"

# (dataset_id, obs file); the ROSMAP-donor snRNA-seq datasets in the vintage
OTHERS = [
    ("Mathys", "Mathys"),
    ("Multiome2025_DLPFC", "Multiome2025_DLPFC"),
    ("Multiome2025_STG", "Multiome2025_STG"),
    ("Green", "Green"),
    ("Cain_2023", "Cain_2023"),
    ("Zhou_2020", "Zhou_2020"),
    ("Mathys_2019", "Mathys_2019"),
    ("PsychAD_ROSMAP", "PsychAD"),
]


def qc_persons(obs_file, dataset_id):
    o = pd.read_parquet(OBJ / f"{obs_file}.obs.parquet",
                        columns=["dataset_id", "person_id", "passed_QC"])
    o = o[o.passed_QC & o.person_id.notna()]
    o = o[o.dataset_id.astype(str) == dataset_id]
    # person_id is categorical: count observed values only, never empty categories
    return o.person_id.astype(str).value_counts()


def main():
    OUT.mkdir(exist_ok=True)
    cu = qc_persons("CUIMC2", "CUIMC2")
    m = pd.DataFrame({"person_id": cu.index, "n_nuclei_CUIMC2": cu.values})
    for ds, f in OTHERS:
        m[ds] = m.person_id.isin(qc_persons(f, ds).index).astype(int)
    assert m.person_id.is_unique
    m.to_csv(OUT / "cuimc2_membership.tsv", sep="\t", index=False)

    sizes = pd.DataFrame({"dataset": [d for d, _ in OTHERS],
                          "shared_with_CUIMC2": [int(m[d].sum()) for d, _ in OTHERS]})
    sizes.to_csv(OUT / "cuimc2_overlap_sizes.tsv", sep="\t", index=False)

    others = [d for d, _ in OTHERS]
    n_any = int((m[others].sum(axis=1) > 0).sum())
    print(f"CUIMC2 persons with QC-passing nuclei: {len(m)}")
    print(sizes.to_string(index=False))
    print(f"in >= 1 other dataset: {n_any}; CUIMC2 only: {len(m) - n_any}")
    same = (m.Multiome2025_DLPFC == m.Multiome2025_STG).all()
    print(f"Multiome2025 DLPFC and STG hold the same CUIMC2 persons: {same}")


if __name__ == "__main__":
    main()
