#!/usr/bin/env python3
"""
Rewrite a PLINK2 .pvar's variant IDs into the canonical chr{N}:{pos}:{REF}:{ALT}
form used across this repo (see harmonize_meta_sumstats.py / liftover_sumstats.py).

Needed because the 1000G reference panel stores IDs WITHOUT the "chr" prefix
("1:10390:A:G") while every cohort VCF in this project stores them WITH it
("chr1:13143:G:C").  Intersecting the two on raw IDs would match nothing.

IDs are rebuilt from the CHROM/POS/REF/ALT columns rather than string-patched,
so an ID that disagrees with its own coordinate columns is corrected rather than
carried forward.

Writes a .pvar with normalized IDs (rows and order unchanged), so it can be
copied directly over the input .pvar.
"""

import argparse
import sys


def norm_id(chrom: str, pos: str, ref: str, alt: str) -> str:
    c = str(chrom).lstrip("chr").lstrip("0") or "0"
    return f"chr{c}:{pos}:{ref.upper()}:{alt.upper()}"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pvar", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    n = 0
    changed = 0
    with open(args.pvar) as fh_in, open(args.out, "w") as fh_out:
        header_written = False
        for line in fh_in:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                fh_out.write("#CHROM\tPOS\tID\tREF\tALT\n")
                header_written = True
                continue
            if not header_written:
                fh_out.write("#CHROM\tPOS\tID\tREF\tALT\n")
                header_written = True

            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            new = norm_id(chrom, pos, ref, alt)
            if new != vid:
                changed += 1
            fh_out.write(f"{chrom}\t{pos}\t{new}\t{ref.upper()}\t{alt.upper()}\n")
            n += 1

    if n == 0:
        sys.exit(f"ERROR: no variant rows read from {args.pvar}")
    print(f"[normalize_ids] {n:,} variants, {changed:,} IDs rewritten", flush=True)


if __name__ == "__main__":
    main()
