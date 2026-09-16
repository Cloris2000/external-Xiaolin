#!/usr/bin/env python3
"""
Liftover a PLINK2 .pvar from hg38 to hg19 for the Figure 1 genotype PCA.

Same conventions as scripts/liftover_sumstats.py (which lifts REGENIE output):
  - pyliftover queries are 0-based, so query POS-1 and add 1 to the result
  - minus-strand chain hits have their alleles reverse-complemented
  - strand-ambiguous SNPs (A/T, C/G) on a minus-strand hit are dropped
  - variant IDs are rewritten as chr{N}:{pos}:{REF}:{ALT}

Writes:
  --out-pvar  a .pvar with hg19 CHROM/POS/ID (rows aligned with --out-keep)
  --out-keep  the ORIGINAL ids of the kept rows, for `plink2 --extract`

The two files are row-aligned: --extract keeps exactly those variants, in the
original file order, so the new .pvar can be copied over the extracted one.

Optionally --restrict-to limits output to variants whose lifted ID appears in a
given ID list (used to keep only the 1000G backbone variants).
"""

import argparse
import sys
from pathlib import Path

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
STRAND_AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}

AUTOSOMES = {str(i) for i in range(1, 23)}


def complement_allele(a: str) -> str:
    return a.translate(_COMP).upper()


def is_strand_ambiguous(ref: str, alt: str) -> bool:
    return frozenset((ref.upper(), alt.upper())) in STRAND_AMBIGUOUS


def norm_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    c = str(chrom).lstrip("chr").lstrip("0") or "0"
    return f"chr{c}:{pos}:{ref.upper()}:{alt.upper()}"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pvar", required=True, help="input .pvar (hg38)")
    p.add_argument("--chain", required=True, help="hg38ToHg19.over.chain.gz")
    p.add_argument("--out-pvar", required=True, help="output .pvar (hg19)")
    p.add_argument("--out-keep", required=True,
                   help="original IDs of kept variants, for --extract")
    p.add_argument("--unmapped-log", required=True)
    p.add_argument("--restrict-to", default=None,
                   help="optional file of hg19 IDs to keep")
    return p.parse_args()


def main():
    args = parse_args()
    try:
        import pyliftover
    except ImportError:
        sys.exit("ERROR: pyliftover is not installed in this interpreter.")

    chain = Path(args.chain)
    if not chain.exists():
        sys.exit(f"ERROR: chain file not found: {chain}")
    print(f"[liftover_pvar] loading chain {chain}", flush=True)
    lo = pyliftover.LiftOver(str(chain))

    restrict = None
    if args.restrict_to:
        with open(args.restrict_to) as fh:
            restrict = {line.strip() for line in fh if line.strip()}
        print(f"[liftover_pvar] restricting to {len(restrict):,} IDs", flush=True)

    n_total = n_kept = n_unmapped = n_ambig = n_multi = n_filtered = 0

    with open(args.pvar) as fh_in, \
         open(args.out_pvar, "w") as fh_pvar, \
         open(args.out_keep, "w") as fh_keep, \
         open(args.unmapped_log, "w") as fh_un:

        fh_un.write("orig_chrom\torig_pos\torig_id\treason\n")
        header_written = False

        for line in fh_in:
            # Drop ## metadata: its contigs and INFO definitions no longer
            # describe the lifted file.
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                fh_pvar.write("#CHROM\tPOS\tID\tREF\tALT\n")
                header_written = True
                continue

            if not header_written:
                fh_pvar.write("#CHROM\tPOS\tID\tREF\tALT\n")
                header_written = True

            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            n_total += 1
            chrom, pos_s, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            pos = int(pos_s)

            query_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            hits = lo.convert_coordinate(query_chrom, pos - 1)

            def log_un(reason):
                fh_un.write(f"{chrom}\t{pos}\t{vid}\t{reason}\n")

            if not hits:
                n_unmapped += 1
                log_un("no_chain_hit")
                continue
            if len(hits) > 1:
                n_multi += 1
                log_un("multiple_chain_hits")
                continue

            new_chrom_full, new_pos0, strand, _score = hits[0]
            new_pos = new_pos0 + 1
            new_chrom = new_chrom_full.lstrip("chr")
            try:
                new_chrom = str(int(new_chrom))
            except ValueError:
                pass

            # The PCA backbone is autosome-only.
            if new_chrom not in AUTOSOMES:
                n_filtered += 1
                log_un("non_autosome_after_lift")
                continue

            new_ref, new_alt = ref, alt
            if strand == "-":
                if is_strand_ambiguous(ref, alt):
                    n_ambig += 1
                    log_un("strand_ambiguous_minus")
                    continue
                new_ref = complement_allele(ref)
                new_alt = complement_allele(alt)

            new_id = norm_id(new_chrom, new_pos, new_ref, new_alt)

            if restrict is not None and new_id not in restrict:
                n_filtered += 1
                continue

            fh_pvar.write(
                f"{new_chrom}\t{new_pos}\t{new_id}\t"
                f"{new_ref.upper()}\t{new_alt.upper()}\n"
            )
            fh_keep.write(f"{vid}\n")
            n_kept += 1

    pct = 100 * n_kept / n_total if n_total else 0
    print(
        f"[liftover_pvar] {n_total:,} in -> {n_kept:,} kept ({pct:.1f}%); "
        f"unmapped={n_unmapped:,} multi={n_multi:,} "
        f"ambig_minus={n_ambig:,} filtered={n_filtered:,}",
        flush=True,
    )
    if n_kept == 0:
        sys.exit("ERROR: no variants survived liftover")


if __name__ == "__main__":
    main()
