#!/usr/bin/env python3
"""
Liftover REGENIE summary statistics from hg38 to hg19.

Reads a .regenie.raw_p file, converts CHROM:GENPOS from hg38 to hg19
using a chain file (via pyliftover), and writes a new file with updated
CHROM, GENPOS, and ID columns.  Variants that cannot be lifted (unmapped
in the chain) are excluded and recorded in an unmapped log.

When the chain indicates a minus-strand alignment, alleles are reverse-
complemented so that positions remain on the + strand (VCF convention).
Strand-ambiguous SNPs (A/T or C/G) that are on a minus-strand chain entry
cannot be unambiguously complemented and are excluded with reason
"strand_ambiguous_minus".

Position conventions
--------------------
  - Input GENPOS   : 1-based (VCF / REGENIE)
  - pyliftover API : 0-based (UCSC BED) → query with GENPOS - 1
  - Output GENPOS  : 1-based → result from pyliftover + 1

Usage
-----
    python liftover_sumstats.py \\
        --input-file  GTEx_v10_Exc_L2_3_IT_step2.regenie.raw_p \\
        --cohort      GTEx_v10 \\
        --chain-file  /path/to/hg38ToHg19.over.chain.gz \\
        --output-file GTEx_v10_Exc_L2_3_IT_step2.hg19_lifted.regenie.raw_p \\
        --unmapped-log GTEx_v10_unmapped.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Nucleotide complement (used when chain maps a region on the minus strand)
# ---------------------------------------------------------------------------
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def complement_allele(allele: str) -> str:
    return allele.translate(_COMP).upper()


STRAND_AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}


def is_strand_ambiguous(ref: str, alt: str) -> bool:
    return frozenset((ref.upper(), alt.upper())) in STRAND_AMBIGUOUS


def norm_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Regenerate variant ID in the same format used by harmonize_meta_sumstats.py."""
    c = str(chrom).lstrip("chr").lstrip("0") or "0"
    return f"chr{c}:{pos}:{ref.upper()}:{alt.upper()}"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-file",   required=True,
                   help=".regenie.raw_p file (hg38 coordinates)")
    p.add_argument("--cohort",       required=True,
                   help="Cohort name (written to unmapped log)")
    p.add_argument("--chain-file",   required=True,
                   help="Path to hg38ToHg19.over.chain.gz")
    p.add_argument("--output-file",  required=True,
                   help="Output file path (hg19 coordinates)")
    p.add_argument("--unmapped-log", required=True,
                   help="TSV log of variants that could not be lifted")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main liftover logic
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    try:
        import pyliftover
    except ImportError:
        sys.exit(
            "ERROR: pyliftover is not installed.  "
            "Run: pip install pyliftover"
        )

    chain_path = Path(args.chain_file)
    if not chain_path.exists():
        sys.exit(f"ERROR: chain file not found: {chain_path}")

    print(f"Loading chain file: {chain_path}", flush=True)
    lo = pyliftover.LiftOver(str(chain_path))

    input_path  = Path(args.input_file)
    output_path = Path(args.output_file)
    unmapped_path = Path(args.unmapped_log)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    counts = dict(total=0, lifted=0, unmapped=0,
                  multi_hit=0, minus_strand_ambiguous=0)

    with (
        open(input_path)   as fh_in,
        open(output_path, "w", newline="") as fh_out,
        open(unmapped_path, "w", newline="") as fh_unmap,
    ):
        reader  = csv.DictReader(fh_in, delimiter=" ")
        fieldnames = reader.fieldnames
        if fieldnames is None:
            sys.exit(f"ERROR: empty or unreadable input file: {input_path}")

        writer  = csv.DictWriter(fh_out,   fieldnames=fieldnames,
                                 delimiter=" ", lineterminator="\n",
                                 extrasaction="ignore")
        unmap_w = csv.DictWriter(fh_unmap,
                                 fieldnames=["cohort", "orig_chrom",
                                             "orig_pos", "id", "reason"],
                                 delimiter="\t", lineterminator="\n")

        writer.writeheader()
        unmap_w.writeheader()

        for row in reader:
            counts["total"] += 1
            chrom_raw = row["CHROM"]
            genpos    = int(row["GENPOS"])
            ref       = row["ALLELE0"]
            alt       = row["ALLELE1"]

            # Ensure "chr" prefix for pyliftover
            query_chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"

            # pyliftover uses 0-based coordinates (UCSC/BED convention)
            hits = lo.convert_coordinate(query_chrom, genpos - 1)

            def _log_unmapped(reason):
                counts["unmapped"] += 1
                unmap_w.writerow({
                    "cohort":      args.cohort,
                    "orig_chrom":  chrom_raw,
                    "orig_pos":    genpos,
                    "id":          row.get("ID", "."),
                    "reason":      reason,
                })

            if not hits:
                _log_unmapped("no_chain_hit")
                continue

            if len(hits) > 1:
                counts["multi_hit"] += 1
                _log_unmapped("multiple_chain_hits")
                continue

            new_chrom_full, new_pos_0based, strand, _score = hits[0]

            # Convert 0-based result back to 1-based VCF position
            new_pos = new_pos_0based + 1

            # Strip "chr" prefix to match REGENIE output format (plain integer)
            new_chrom_stripped = new_chrom_full.lstrip("chr")
            # Preserve original zero-padding style: REGENIE uses plain integers
            # so strip leading zeros (e.g. "chr01" → "1")
            try:
                new_chrom_stripped = str(int(new_chrom_stripped))
            except ValueError:
                pass  # X, Y, MT — keep as-is

            # Handle minus-strand chain mappings: complement alleles
            new_ref, new_alt = ref, alt
            if strand == "-":
                if is_strand_ambiguous(ref, alt):
                    counts["minus_strand_ambiguous"] += 1
                    _log_unmapped("strand_ambiguous_minus")
                    continue
                new_ref = complement_allele(ref)
                new_alt = complement_allele(alt)

            new_id = norm_id(new_chrom_stripped, new_pos, new_ref, new_alt)

            new_row = dict(row)
            new_row["CHROM"]   = new_chrom_stripped
            new_row["GENPOS"]  = new_pos
            new_row["ID"]      = new_id
            new_row["ALLELE0"] = new_ref.upper()
            new_row["ALLELE1"] = new_alt.upper()

            writer.writerow(new_row)
            counts["lifted"] += 1

    pct_lifted = 100 * counts["lifted"] / counts["total"] if counts["total"] else 0
    print(
        f"Cohort {args.cohort}: "
        f"{counts['total']:,} variants in → "
        f"{counts['lifted']:,} lifted ({pct_lifted:.1f}%) | "
        f"unmapped {counts['unmapped']:,} "
        f"(no_hit: {counts['unmapped'] - counts['multi_hit'] - counts['minus_strand_ambiguous']:,}  "
        f"multi_hit: {counts['multi_hit']:,}  "
        f"strand_ambig_minus: {counts['minus_strand_ambiguous']:,})",
        flush=True,
    )


if __name__ == "__main__":
    main()
