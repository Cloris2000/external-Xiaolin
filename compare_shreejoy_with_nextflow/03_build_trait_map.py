#!/usr/bin/env python3
"""Stage A / step 3

Map Shreejoy's SEA-AD supertypes onto Xiaolin's 19 broad MGP classes, via the
subclass level of refs/taxonomy_DFC_2026.tsv.

The mapping is many-to-one by construction: several supertypes roll up into one
subclass, and several subclasses roll up into one MGP class (the three IT
subclasses all become `IT`, Sst and Sst Chodl both become `SST`). That collapse
is the thing the comparison is about, so it is made explicit here rather than
buried in a join.

Writes three files:

  data/trait_map.tsv        supertype -> subclass -> MGP class, one row each
  data/subclass_map.tsv     subclass -> MGP class, with core-supertype counts
  data/unmapped.tsv         everything with no counterpart on the other side

Read-only. Stdlib only.
"""

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TAXONOMY = Path("/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/"
                "celltype-composition/refs/taxonomy_DFC_2026.tsv")
DATA = HERE / "data"

# Xiaolin's 19 MGP traits, exactly as they appear in her result filenames.
MY_TRAITS = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
    "L5.ET", "L6b", "L6.CT", "LAMP5", "Microglia", "Oligodendrocyte", "OPC",
    "PAX6", "Pericyte", "PVALB", "SST", "VIP", "VLMC",
]

# SEA-AD subclass -> MGP class. Every one of the 24 subclasses is listed, so a
# new taxonomy version that adds one will fail the completeness check below
# rather than silently drop it.
SUBCLASS_TO_MGP = {
    "Astrocyte": "Astrocyte",
    "Chandelier": "PVALB",            # chandelier cells are a Pvalb type
    "Endothelial": "Endothelial",
    "Immune": "Microglia",            # SEA-AD Immune is microglia plus PVM
    "L2/3 IT": "IT",
    "L4 IT": "L4.IT",
    "L5/6 NP": "L5.6.NP",
    "L5 ET": "L5.ET",
    "L5 IT": "IT",
    "L6b": "L6b",
    "L6 CT": "L6.CT",
    "L6 IT": "IT",
    "L6 IT Car3": "L5.6.IT.Car3",
    "Lamp5": "LAMP5",
    "Lamp5 Lhx6": "LAMP5",
    "Oligodendrocyte": "Oligodendrocyte",
    "OPC": "OPC",
    "Pax6": "PAX6",
    "Pvalb": "PVALB",
    "Sncg": None,                     # no MGP counterpart
    "Sst": "SST",
    "Sst Chodl": "SST",
    "Vip": "VIP",
    # SEA-AD pools VLMC and pericytes into one subclass while MGP splits them.
    # Assigned to VLMC; Pericyte is therefore left without a counterpart.
    "VLMC & Perivascular": "VLMC",
}


def norm(label: str) -> str:
    """Supertype label as it appears in REGENIE filenames.

    `Sst Chodl_3-SEAAD` on disk in the taxonomy and arm_agreement tables becomes
    `Sst_Chodl_3_SEAAD` as a REGENIE trait name.
    """
    return label.replace(" ", "_").replace("-", "_").replace("/", "_")


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)

    with TAXONOMY.open() as fh:
        tax = list(csv.DictReader(fh, delimiter="\t"))

    subclasses = {r["subclass_label"] for r in tax}
    unknown = subclasses - set(SUBCLASS_TO_MGP)
    if unknown:
        sys.exit(f"FATAL: taxonomy has subclasses this script does not map: "
                 f"{sorted(unknown)}")

    with (DATA / "trait_map.tsv").open("w") as out:
        out.write("supertype_label\tregenie_trait\tsubclass_label\t"
                  "class_label\tcompartment\tcore\tmgp_trait\n")
        for r in tax:
            mgp = SUBCLASS_TO_MGP[r["subclass_label"]] or "NA"
            out.write("\t".join([
                r["supertype_label"], norm(r["supertype_label"]),
                r["subclass_label"], r["class_label"], r["compartment"],
                r["core"], mgp,
            ]) + "\n")

    # Per subclass: how many supertypes exist, and how many survived the
    # portability filter into the CORE set Shreejoy actually scanned.
    counts: dict[str, list[int]] = {}
    for r in tax:
        c = counts.setdefault(r["subclass_label"], [0, 0])
        c[0] += 1
        if r["core"] == "True":
            c[1] += 1

    with (DATA / "subclass_map.tsv").open("w") as out:
        out.write("subclass_label\tmgp_trait\tn_supertypes\tn_core\n")
        for sub in sorted(counts):
            out.write(f"{sub}\t{SUBCLASS_TO_MGP[sub] or 'NA'}\t"
                      f"{counts[sub][0]}\t{counts[sub][1]}\n")

    # An MGP trait has a counterpart only if some subclass mapping to it kept at
    # least one core supertype.
    core_by_mgp: dict[str, int] = {t: 0 for t in MY_TRAITS}
    for sub, (_, n_core) in counts.items():
        mgp = SUBCLASS_TO_MGP[sub]
        if mgp:
            core_by_mgp[mgp] += n_core

    with (DATA / "unmapped.tsv").open("w") as out:
        out.write("side\tname\treason\n")
        for t in MY_TRAITS:
            if core_by_mgp[t] == 0:
                out.write(f"mine\t{t}\tno core supertype survived the "
                          f"portability filter\n")
        for sub in sorted(counts):
            if SUBCLASS_TO_MGP[sub] is None:
                out.write(f"his\t{sub}\tno MGP class corresponds\n")

    no_counterpart = [t for t in MY_TRAITS if core_by_mgp[t] == 0]
    print(f"wrote {DATA}/trait_map.tsv, subclass_map.tsv, unmapped.tsv",
          file=sys.stderr)
    print(f"  {len(MY_TRAITS) - len(no_counterpart)} of {len(MY_TRAITS)} MGP "
          f"traits have a core counterpart", file=sys.stderr)
    print(f"  without one: {', '.join(no_counterpart)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
