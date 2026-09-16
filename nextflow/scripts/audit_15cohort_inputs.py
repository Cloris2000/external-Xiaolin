#!/usr/bin/env python3
"""Inventory RNA, VCF, SCC phenotypes and .raw_p files for the 15 bulk cohorts."""

from __future__ import annotations

import csv
from pathlib import Path

NF = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
SCC = Path("/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow/results")
MAP = NF / "validation/cohort_trillium_paths.tsv"
OUT = Path("/scratch/zhoux156/results/validation/cohort_input_audit.tsv")


def status(path: Path | None, min_bytes: int = 1) -> str:
    if path is None or str(path) == "":
        return "NA"
    if path.is_file() and path.stat().st_size >= min_bytes:
        return f"OK:{path.stat().st_size}"
    if path.is_dir():
        return "DIR"
    return "MISSING"


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(MAP.open(), delimiter="\t"))
    print(f"{'cohort':<18} {'rna':<10} {'meta':<10} {'bio':<10} {'vcf_n':>5} {'pheno':<8} {'prop':<8} {'raw_p':>6} {'pgen':<6}")
    out_rows = []
    for r in rows:
        rna = Path(r["count_matrix"])
        meta = Path(r["metadata"])
        bio = Path(r["biospecimen"]) if r["biospecimen"] else None
        vcf_dir = Path(r["vcf_dir"])
        n_vcf = len(list(vcf_dir.glob(r["vcf_glob"]))) if vcf_dir.is_dir() else 0
        scc = SCC / r["scc_results"]
        pheno = scc / "phenotypes_RINT.txt"
        prop = scc / "cell_proportions.csv"
        n_raw = len(list((scc / "regenie_step2").glob("*.regenie.raw_p"))) if (scc / "regenie_step2").is_dir() else 0
        n_pgen = len(list(scc.glob("*.QC.final.pgen")))
        print(
            f"{r['cohort']:<18} {status(rna):<10} {status(meta):<10} {status(bio):<10} "
            f"{n_vcf:>5} {status(pheno):<8} {status(prop):<8} {n_raw:>5}/19 {n_pgen:>5}"
        )
        out_rows.append({
            "cohort": r["cohort"],
            "rna": status(rna),
            "metadata": status(meta),
            "biospecimen": status(bio),
            "vcf_dir": str(vcf_dir),
            "n_vcf": n_vcf,
            "scc_pheno": status(pheno),
            "scc_prop": status(prop),
            "n_raw_p": n_raw,
            "n_pgen": n_pgen,
            "config": r["config"],
        })
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    n_ok = sum(1 for x in out_rows if x["rna"].startswith("OK") and x["metadata"].startswith("OK") and int(x["n_vcf"]) >= 22 and int(x["n_raw_p"]) == 19)
    print(f"\ncohorts with RNA+meta+>=22 VCF+19 raw_p: {n_ok}/{len(out_rows)}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
