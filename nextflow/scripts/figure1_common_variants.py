#!/usr/bin/env python3
"""
Intersect the 1000G backbone with all 15 cohort .pvar files for the Figure 1 PCA.

A variant is kept only if, in EVERY cohort and in the reference, it is present
with the same REF/ALT orientation.  Because variant IDs already encode
chr:pos:REF:ALT, an exact ID match implies identical orientation; this script
additionally re-checks REF/ALT from the columns so a malformed ID cannot slip
through.

Strand-ambiguous SNPs (A/T, C/G) are dropped: with only allele codes to go on,
a flipped-strand site is indistinguishable from a correctly oriented one, and a
mis-oriented variant would distort the projected PCs.

Writes the shared ID list (--out) and a per-cohort accounting table (--report).
"""

import argparse
import sys
from pathlib import Path

STRAND_AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}


def is_strand_ambiguous(ref: str, alt: str) -> bool:
    return frozenset((ref.upper(), alt.upper())) in STRAND_AMBIGUOUS


def read_pvar(path):
    """Return {variant_id: (REF, ALT)} for a PLINK2 .pvar."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            out[f[2]] = (f[3].upper(), f[4].upper())
    return out


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref-pvar", required=True, help="1000G backbone .pvar (hg19)")
    p.add_argument("--cohort-pvars", nargs="+", required=True,
                   help="per-cohort .pvar files on the backbone")
    p.add_argument("--out", required=True, help="shared variant ID list")
    p.add_argument("--report", required=True, help="per-cohort TSV report")
    return p.parse_args()


def main():
    args = parse_args()

    ref = read_pvar(args.ref_pvar)
    print(f"[common_variants] reference: {len(ref):,} variants", flush=True)

    # Drop strand-ambiguous sites from the reference up front.
    ambig = {v for v, (r, a) in ref.items() if is_strand_ambiguous(r, a)}
    shared = set(ref) - ambig
    print(f"[common_variants] dropped {len(ambig):,} strand-ambiguous "
          f"(A/T, C/G) reference variants", flush=True)

    rows = []
    for pv in sorted(args.cohort_pvars):
        cohort = Path(pv).stem.replace("coh_", "")
        cv = read_pvar(pv)
        # Same ID and same REF/ALT as the reference.
        ok = {v for v, (r, a) in cv.items()
              if v in ref and ref[v] == (r, a) and v not in ambig}
        before = len(shared)
        shared &= ok
        rows.append((cohort, len(cv), len(ok), before, len(shared)))
        print(f"[common_variants] {cohort:18s} pvar={len(cv):8,d} "
              f"match_ref={len(ok):8,d} running_shared={len(shared):8,d}",
              flush=True)

    if not shared:
        sys.exit("ERROR: no variants shared across all cohorts")

    # Write in reference coordinate order for reproducibility.
    order = [v for v in ref if v in shared]
    with open(args.out, "w") as fh:
        fh.write("\n".join(order) + "\n")

    with open(args.report, "w") as fh:
        fh.write("cohort\tn_cohort_variants\tn_matching_reference\t"
                 "shared_before\tshared_after\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
        fh.write(f"FINAL\t\t\t\t{len(shared)}\n")

    print(f"[common_variants] final shared set: {len(shared):,} variants",
          flush=True)
    if len(shared) < 20000:
        print(f"[common_variants] WARNING: only {len(shared):,} shared variants; "
              "PC estimates may be unstable", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
