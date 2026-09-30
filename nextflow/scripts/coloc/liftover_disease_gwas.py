#!/usr/bin/env python3
"""
Lift a standardized disease GWAS table from GRCh38 to hg19 (GRCh37).

Why this exists
---------------
02_download_disease_gwas.sh writes <Disease>_hg19.tsv, but for the GWAS-Catalog
"harmonised" sources (AD_Bellenguez2022, LBD_Chia2021, PD_Nalls2019) it copies
hm_pos / base_pair_location straight through -- and GWAS-Catalog harmonised
coordinates are GRCh38, not GRCh37.  Verified 2026-09-16: the "_hg19" position
equals the raw hm_pos for the same rsID, and rs429358 (APOE) sits at its hg38
coordinate in PD_Nalls2019_hg19.tsv.  Because the CTP meta-analysis is hg19,
coloc matched almost nothing for those three diseases (AD 15 tests vs 462 for
BD/MDD/SCZ); the PGC daner sources (BD/MDD/SCZ) are natively hg19 and unaffected.

Format (unchanged):  snp chr pos ref alt beta se p N Ncase Ncont
  snp is rebuilt as chr<chr>:<pos>:<ref>:<alt> from the lifted position.
  On minus-strand chain blocks ref/alt are reverse-complemented; strand-ambiguous
  SNPs (A/T, C/G) on such blocks are dropped as unresolvable.  Unmapped and
  multiply-mapped positions are dropped and counted.  'ref'/'alt' keep the
  writer's semantics (ref = effect allele, alt = other allele) -- 03_run_coloc.R
  only uses them for allele matching, and coloc.abf is sign-invariant.

Usage:
  liftover_disease_gwas.py --input AD_Bellenguez2022_hg19.tsv --chain hg38ToHg19.over.chain.gz \
                           --output out/AD_Bellenguez2022/AD_Bellenguez2022_hg19.tsv
"""
import argparse
import csv
import sys
from pathlib import Path

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
AMBIG = {frozenset(("A", "T")), frozenset(("C", "G"))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    try:
        from pyliftover import LiftOver
    except ImportError:
        sys.exit("pyliftover not installed in this python")
    lo = LiftOver(a.chain)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)

    n = dict(total=0, lifted=0, unmapped=0, multi=0, ambig_minus=0, other_chrom=0, bad=0)
    with open(a.input, newline="") as fi, open(a.output, "w", newline="") as fo:
        rd = csv.DictReader(fi, delimiter="\t")
        cols = ["snp", "chr", "pos", "ref", "alt", "beta", "se", "p", "N", "Ncase", "Ncont"]
        missing = [c for c in cols if c not in rd.fieldnames]
        if missing:
            sys.exit(f"{a.input}: missing columns {missing}")
        wr = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        wr.writeheader()
        for row in rd:
            n["total"] += 1
            chrom = str(row["chr"]).replace("chr", "")
            try:
                pos = int(row["pos"])
            except ValueError:
                n["bad"] += 1
                continue
            hits = lo.convert_coordinate(f"chr{chrom}", pos - 1)  # pyliftover is 0-based
            if not hits:
                n["unmapped"] += 1
                continue
            if len(hits) > 1:
                n["multi"] += 1
                continue
            new_chrom, new_pos0, strand, _ = hits[0]
            new_chrom = new_chrom.replace("chr", "")
            if new_chrom != chrom:  # keep it simple: no cross-chromosome lifts
                n["other_chrom"] += 1
                continue
            ref, alt = row["ref"].upper(), row["alt"].upper()
            if strand == "-":
                if frozenset((ref, alt)) in AMBIG:
                    n["ambig_minus"] += 1
                    continue
                ref, alt = ref.translate(COMP), alt.translate(COMP)
            new_pos = new_pos0 + 1
            row.update(chr=new_chrom, pos=str(new_pos), ref=ref, alt=alt,
                       snp=f"chr{new_chrom}:{new_pos}:{ref}:{alt}")
            wr.writerow(row)
            n["lifted"] += 1

    pct = 100 * n["lifted"] / n["total"] if n["total"] else 0
    print(f"{Path(a.input).name}: {n['total']:,} in -> {n['lifted']:,} lifted ({pct:.2f}%) | "
          f"unmapped {n['unmapped']:,}  multi {n['multi']:,}  other_chrom {n['other_chrom']:,}  "
          f"strand_ambig_minus {n['ambig_minus']:,}  bad {n['bad']:,}")
    with open(a.output + ".liftover_summary.txt", "w") as fs:
        for k, v in n.items():
            fs.write(f"{k}\t{v}\n")


if __name__ == "__main__":
    main()
