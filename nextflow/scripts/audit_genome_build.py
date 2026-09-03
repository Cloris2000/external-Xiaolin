#!/usr/bin/env python3
"""Decide whether REGENIE summary stats are on hg19/b37 or hg38.

Two independent signals, so a cohort with sparse coverage still gets a call:

  1. Landmark variants whose hg19 and hg38 coordinates differ by far more than
     any plausible annotation slop.  A hit at one coordinate and not the other
     is decisive on its own.
  2. Maximum observed position per chromosome.  hg19 chromosomes are longer than
     their hg38 counterparts, so sumstats exceeding the hg38 length cannot be
     hg38.  This is one-sided: it can rule hg38 out but cannot rule hg19 out.

Usage:
    scripts/audit_genome_build.py FILE [FILE ...]
    scripts/audit_genome_build.py --label-from-path 'results/{cohort}/...' FILE...

Expects whitespace-delimited REGENIE output with CHROM and GENPOS columns.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys

# (rsid, chrom, hg19_pos, hg38_pos)
LANDMARKS = [
    ("rs429358 APOE",  "19", 45411941, 44908684),
    ("rs6265 BDNF",    "11", 27679916, 27658369),
    ("rs9939609 FTO",  "16", 53820527, 53786615),
    ("rs1801133 MTHFR", "1",  11856378, 11796321),
    ("rs4988235 LCT",   "2", 136608646, 135851076),
]

# Chromosome lengths; hg19 is longer than hg38 for every autosome here.
CHROM_LEN = {
    "1":  (249250621, 248956422), "2":  (243199373, 242193529),
    "3":  (198022430, 198295559), "4":  (191154276, 190214555),
    "6":  (171115067, 170805979), "11": (135006516, 135086622),
    "16": (90354753,  90338345),  "19": (59128983,  58617616),
}


def opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def audit(path: str):
    wanted = {(c, p19) for _, c, p19, _ in LANDMARKS} | {(c, p38) for _, c, _, p38 in LANDMARKS}
    hits: set[tuple[str, int]] = set()
    max_pos: dict[str, int] = {}
    n = 0

    with opener(path) as fh:
        header = fh.readline().split()
        try:
            ci, pi = header.index("CHROM"), header.index("GENPOS")
        except ValueError:
            return None, f"no CHROM/GENPOS columns in header: {header[:6]}"
        for line in fh:
            f = line.split()
            if len(f) <= pi:
                continue
            n += 1
            c = f[ci]
            try:
                p = int(f[pi])
            except ValueError:
                continue
            if p > max_pos.get(c, 0):
                max_pos[c] = p
            if (c, p) in wanted:
                hits.add((c, p))

    lm19 = sum((c, p) in hits for _, c, p, _ in LANDMARKS)
    lm38 = sum((c, p) in hits for _, c, _, p in LANDMARKS)

    # A position beyond the hg38 chromosome length rules hg38 out.
    over_hg38 = sum(1 for c, m in max_pos.items()
                    if c in CHROM_LEN and m > CHROM_LEN[c][1])

    if lm19 and not lm38:
        call = "hg19"
    elif lm38 and not lm19:
        call = "hg38"
    elif lm19 and lm38:
        call = "MIXED?"
    elif over_hg38:
        call = "hg19"          # exceeds hg38 lengths, no landmarks seen
    else:
        call = "UNKNOWN"

    detail = (f"n={n:>9,}  landmarks hg19={lm19}/{len(LANDMARKS)} "
              f"hg38={lm38}/{len(LANDMARKS)}  chroms_over_hg38_len={over_hg38}")
    return call, detail


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    width = max(len(os.path.basename(f)) for f in args.files)
    rc = 0
    for path in args.files:
        if not os.path.exists(path):
            print(f"{os.path.basename(path):<{width}}  MISSING")
            rc = 1
            continue
        call, detail = audit(path)
        print(f"{os.path.basename(path):<{width}}  {str(call):<8}  {detail}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
