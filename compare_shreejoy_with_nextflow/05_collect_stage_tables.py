#!/usr/bin/env python3
"""Stage A / step 5

Assemble the tidy per-stage tables. Stages 1, 3 and 4 need no new compute: both
pipelines already wrote the comparable quantity, so this is joining and naming,
not calculation.

  stage1_deconvolution.tsv  bulk-vs-snRNA agreement per cell type, both sides
  stage3_discovery.tsv      hits per trait with heterogeneity, both sides
  stage3_meta_vs_mega.tsv   his pooled mega against his own ancestry meta
  stage4_coloc.tsv          PP.H4 per (locus, disease, cell type), both sides

Read-only. Stdlib only.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

MINE_ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
MINE_DECONV = MINE_ROOT / "manuscript_figure" / "figure2_celltype_accuracy.tsv"
MINE_HET = (MINE_ROOT / "results" / "meta_analysis_15cohorts_hg19_v2" /
            "heterogeneity" / "meta_heterogeneity_summary.tsv")
MINE_COLOC = Path("/scratch/zhoux156/results/downstream_v2/coloc/"
                  "coloc_results_full")

HIS_ROOT = Path("/scratch/shreejoy/ctpgwas/results")

# figure2_celltype_accuracy.tsv uses lowercase slugs; results files use these.
SLUG_TO_TRAIT = {
    "astrocyte": "Astrocyte", "endothelial": "Endothelial", "it": "IT",
    "l4_it": "L4.IT", "l5_6_it_car3": "L5.6.IT.Car3", "l5_6_np": "L5.6.NP",
    "l5_et": "L5.ET", "l6b": "L6b", "l6_ct": "L6.CT", "lamp5": "LAMP5",
    "microglia": "Microglia", "oligodendrocyte": "Oligodendrocyte",
    "opc": "OPC", "pax6": "PAX6", "pericyte": "Pericyte", "pvalb": "PVALB",
    "sst": "SST", "vip": "VIP", "vlmc": "VLMC",
}

# The two loci, in hg19, so my coloc results can be assigned to them. Offsets
# come from 04_align_coords.py.
MY_LOCUS_WINDOWS = {
    "TMEM106B": (7, 11_900_000, 12_600_000),
    "GRN/FAM171A2": (17, 42_100_000, 42_800_000),
}

# Disease labels differ between the two panels; only these are shared.
DISEASE_ALIAS = {
    "AD_Bellenguez2022": "AD", "PD_Nalls2019": "PD", "LBD_Chia2021": "LBD",
    "SCZ_Trubetskoy2022": "SCZ", "SCZ_PGC3": "SCZ", "SCZ_Bigdeli2026": "SCZ",
    "BD_bip2024": "BD", "BD_2024": "BD",
    "MDD_MDD2025": "MDD", "MDD_2025_Clin": "MDD",
    "ADHD_2022": "ADHD",
}


def read_tsv(path: Path, delim: str = "\t") -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def norm(label: str) -> str:
    """Same normalisation as 03_build_trait_map.py.

    Needed because the three sources spell supertypes three ways: the taxonomy
    has `L2/3 IT_10` and `Sst Chodl_3-SEAAD`, arm_agreement.csv has `L2_3_IT_10`
    and `Sst_Chodl_3-SEAAD`, and REGENIE filenames have `Sst_Chodl_3_SEAAD`.
    """
    return label.replace(" ", "_").replace("-", "_").replace("/", "_")


def trait_map() -> dict[str, dict]:
    """Supertype lookup, keyed on every spelling that appears on disk."""
    index: dict[str, dict] = {}
    for r in read_tsv(DATA / "trait_map.tsv"):
        index[r["supertype_label"]] = r
        index[r["regenie_trait"]] = r
        index[norm(r["supertype_label"])] = r
    return index


def lookup(tmap: dict[str, dict], name: str) -> dict | None:
    return tmap.get(name) or tmap.get(norm(name))


def his_scanned_traits() -> list[str]:
    """The traits Shreejoy actually scanned, from the extracted pooled arm."""
    with (DATA / "his_loci.tsv").open() as fh:
        return sorted({r["trait"] for r in csv.DictReader(fh, delimiter="\t")
                       if r["arm"] == "pooled"})


def stage1() -> None:
    """Bulk-vs-snRNA agreement, both pipelines, on one table.

    Not like for like, and deliberately so: mine is 19 broad classes, his is 33
    supertypes, and a finer type is intrinsically harder to deconvolve. The
    comparison is therefore conservative in his favour.
    """
    tmap = trait_map()
    rows = []

    for r in read_tsv(MINE_DECONV):
        slug = r["cell_type"]
        if slug not in SLUG_TO_TRAIT:
            continue                      # skips a repeated header row
        rows.append({
            "pipeline": "mine", "resolution": "broad class",
            "label": SLUG_TO_TRAIT[slug], "mgp_trait": SLUG_TO_TRAIT[slug],
            "r": r["mean_r"], "r_disattenuated": "NA",
            "n": r["total_donors"],
        })

    for r in read_tsv(HIS_ROOT / "person_pheno" / "arm_agreement.csv", ","):
        t = lookup(tmap, r["cell_type"])
        if t is None:
            sys.exit(f"FATAL: {r['cell_type']} in arm_agreement.csv has no "
                     f"taxonomy entry")
        rows.append({
            "pipeline": "his", "resolution": "supertype",
            "label": r["cell_type"],
            "mgp_trait": t["mgp_trait"] if t else "NA",
            "r": r["r"], "r_disattenuated": r["r_disatt"], "n": r["n"],
        })

    write(DATA / "stage1_deconvolution.tsv", rows)


def stage3_discovery() -> None:
    """Hits per trait, with the heterogeneity that only a meta-analysis exposes."""
    tmap = trait_map()
    rows = []

    for r in read_tsv(MINE_HET):
        rows.append({
            "pipeline": "mine", "trait": r["cell_type"],
            "mgp_trait": r["cell_type"],
            "n_gw": r["n_gw"], "n_gw_high_het": r["n_gw_high_het"],
            "n_variants": r["n_variants"],
        })

    # His scan is pooled, so there is no per-variant I2 to report; hits are
    # counted from his own locus table, over the traits he actually scanned.
    per_trait: dict[str, int] = defaultdict(int)
    for r in read_tsv(HIS_ROOT / "loci" / "annotated_hits.csv", ","):
        if r["arm"] == "all":
            per_trait[r["trait"]] += 1
    for trait in his_scanned_traits():
        t = lookup(tmap, trait)
        rows.append({
            "pipeline": "his", "trait": trait,
            "mgp_trait": t["mgp_trait"] if t else "NA",
            "n_gw": per_trait.get(trait, 0), "n_gw_high_het": "NA",
            "n_variants": "NA",
        })

    write(DATA / "stage3_discovery.tsv", rows)


def stage3_meta_vs_mega() -> None:
    """His pooled mega against his own inverse-variance ancestry meta.

    Same people, same variants, two models. This is the only contrast available
    that isolates meta-versus-mega from everything else that differs.
    """
    rows = []
    for r in read_tsv(HIS_ROOT / "loci" / "annotated_hits.csv", ","):
        if r["arm"] != "all":
            continue
        rows.append({
            "model": "mega (pooled)", "locus": r["locus"], "trait": r["trait"],
            "CHROM": r["CHROM"], "POS": r["GENPOS"], "P": r["P"],
            "BETA": r["BETA"], "I2": "NA",
        })
    for r in read_tsv(HIS_ROOT / "meta_ancestry" / "hits.csv", ","):
        rows.append({
            "model": "meta (ancestry IVW)", "locus": "NA", "trait": r["trait"],
            "CHROM": r["CHROM"], "POS": r["GENPOS"], "P": r["P"],
            "BETA": r["BETA"], "I2": r["I2"],
        })
    write(DATA / "stage3_meta_vs_mega.tsv", rows)


def assign_locus(chrom: int, pos: int) -> str | None:
    for name, (c, lo, hi) in MY_LOCUS_WINDOWS.items():
        if chrom == c and lo <= pos <= hi:
            return name
    return None


def stage4() -> None:
    """Colocalisation posteriors, both pipelines, restricted to shared diseases.

    Nearly controlled: his refs/disease_gwas/panel.tsv points AD, PD and LBD at
    Xiaolin's own standardised hg19 files, so the disease side is the same data.
    """
    tmap = trait_map()
    rows = []

    for path in sorted(MINE_COLOC.glob("*_coloc_results.tsv")):
        for r in read_tsv(path):
            # locus_id is celltype_chrN_pos in hg19
            parts = r["locus_id"].rsplit("_", 2)
            if len(parts) != 3 or not parts[1].startswith("chr"):
                continue
            chrom, pos = int(parts[1][3:]), int(parts[2])
            rows.append({
                "pipeline": "mine", "locus": assign_locus(chrom, pos) or "other",
                "locus_id": r["locus_id"], "CHROM": chrom, "POS_hg19": pos,
                "disease": DISEASE_ALIAS.get(r["disease"], r["disease"]),
                "disease_raw": r["disease"], "cell_type": r["cell_type"],
                "mgp_trait": r["cell_type"], "PP.H4": r["PP.H4"],
                "n_snps": r["n_snps"],
            })

    for r in read_tsv(HIS_ROOT / "coloc" / "coloc_results.tsv"):
        locus = r["locus"].replace("GRN_FAM171A2", "GRN/FAM171A2")
        t = lookup(tmap, r["trait"])
        rows.append({
            "pipeline": "his", "locus": locus, "locus_id": r["locus"],
            "CHROM": "NA", "POS_hg19": "NA",
            "disease": DISEASE_ALIAS.get(r["disease"], r["disease"]),
            "disease_raw": r["disease"], "cell_type": r["trait"],
            "mgp_trait": t["mgp_trait"] if t else "NA",
            "PP.H4": r["PP.H4"], "n_snps": r["n"],
        })

    write(DATA / "stage4_coloc.tsv", rows)


def write(path: Path, rows: list[dict]) -> None:
    if not rows:
        sys.exit(f"FATAL: no rows for {path.name}")
    with path.open("w") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path.name}: {len(rows)} rows", file=sys.stderr)


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    stage1()
    stage3_discovery()
    stage3_meta_vs_mega()
    stage4()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
