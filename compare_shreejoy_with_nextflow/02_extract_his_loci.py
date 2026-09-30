#!/usr/bin/env python3
"""Stage A / step 2

Pull the TMEM106B and GRN windows out of Shreejoy's phase-2 REGENIE step 2, for
the three arms that matter to the comparison: pooled, snRNA-seq, deconvolved bulk.

Two guards, both of which have bitten this tree before:

  * The step2 directory also holds `*.stale_pre_repin_20260922T0917` copies left
    by the caching defect described in ctp-gwas/src/ctpgwas/traits.py, where
    REGENIE reused step 1 and step 2 never swept dropped trait files. Only the
    current directories are read.
  * The trait list on disk is checked against loci/annotated_hits.csv, so a
    stale or partial arm is caught here rather than silently plotted.

REGENIE reports LOG10P, not P, so P is recovered as 10**-LOG10P. BETA is oriented
to ALLELE1, which is ALT in the chr:pos:REF:ALT ID scheme.

Read-only. Stdlib only.
"""

import os
import sys
from pathlib import Path

HIS_ROOT = Path("/scratch/shreejoy/ctpgwas/results")
OUT = Path(__file__).resolve().parent / "data" / "his_loci.tsv"

ARMS = ["all_shrunk", "sn_shrunk", "bulk_shrunk"]
ARM_LABEL = {"all_shrunk": "pooled", "sn_shrunk": "snRNA", "bulk_shrunk": "bulk"}

# hg38 windows, matching the hg19 windows in 01_extract_my_loci.sh
WINDOWS = {
    7: (11_860_000, 12_560_000),   # TMEM106B hg38 chr7:12,211,240-12,243,367
    17: (44_020_000, 44_720_000),  # GRN      hg38 chr17:44,345,086-44,353,106
}


def expected_traits() -> set[str]:
    """Trait names Shreejoy's own locus table refers to."""
    path = HIS_ROOT / "loci" / "annotated_hits.csv"
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split(",")
        i = header.index("trait")
        return {line.split(",")[i] for line in fh if line.strip()}


def arm_traits(arm_dir: Path) -> set[str]:
    return {p.name.split("_", 1)[1].removesuffix(".regenie")
            for p in arm_dir.glob("chr7_*.regenie")}


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    want = expected_traits()

    rows = 0
    with OUT.open("w") as out:
        out.write("arm\ttrait\tvariant\tCHROM\tPOS_hg38\tREF\tALT\t"
                  "BETA\tSE\tP\tA1FREQ\tN\n")

        for arm in ARMS:
            arm_dir = HIS_ROOT / "step2" / arm
            if not arm_dir.is_dir():
                sys.exit(f"FATAL: missing arm directory {arm_dir}")

            have = arm_traits(arm_dir)
            missing = want - have
            if missing:
                sys.exit(f"FATAL: {arm} is missing traits named in "
                         f"annotated_hits.csv: {sorted(missing)}")
            print(f"  {arm}: {len(have)} traits", file=sys.stderr)

            for chrom, (lo, hi) in WINDOWS.items():
                for path in sorted(arm_dir.glob(f"chr{chrom}_*.regenie")):
                    trait = path.name.split("_", 1)[1].removesuffix(".regenie")
                    with path.open() as fh:
                        head = fh.readline().split()
                        ix = {c: i for i, c in enumerate(head)}
                        for line in fh:
                            f = line.split()
                            pos = int(f[ix["GENPOS"]])
                            if not (lo <= pos <= hi):
                                continue
                            vid = f[ix["ID"]]
                            log10p = f[ix["LOG10P"]]
                            p = "NA" if log10p in ("NA", "") else f"{10 ** -float(log10p):.6g}"
                            out.write("\t".join([
                                ARM_LABEL[arm], trait, vid, str(chrom), str(pos),
                                f[ix["ALLELE0"]], f[ix["ALLELE1"]],
                                f[ix["BETA"]], f[ix["SE"]], p,
                                f[ix["A1FREQ"]], f[ix["N"]],
                            ]) + "\n")
                            rows += 1

    print(f"wrote {OUT}: {rows} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
