#!/usr/bin/env python3
"""
Annotate a METAL .tbl with explicit REF/ALT and effect-allele columns.

MarkerName is chr:pos:REF:ALT (REF/ALT = VCF REF/ALT at genotype QC; the
hg38 cohorts are lifted to hg19 with alleles complemented on minus-strand
chains, so the pair stays consistent).  METAL's Allele1 is always the allele
that Effect and Freq1 describe, but METAL canonically sorts allele pairs before
writing (only A/C, A/G, A/T, C/G, T/C, T/G appear), so Allele1 is REF for some
variants and ALT for others.

Output (tab-separated, all original columns first, then):
  CHROM POS REF ALT
  EFFECT_ALLELE OTHER_ALLELE     as METAL reported them (upper case)
  EFFECT_ALLELE_IS               REF | ALT | UNMATCHED
  BETA_ALT FREQ_ALT MINFREQ_ALT MAXFREQ_ALT DIRECTION_ALT
                                 effect, frequency and per-cohort direction
                                 re-expressed for the ALT allele, so every
                                 variant is in the same frame.
The original .tbl is not modified.
"""
import argparse
import csv
import sys


def flip_dir(s):
    return s.translate(str.maketrans("+-", "-+"))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tbl", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()

    counts = dict(total=0, effect_is_alt=0, effect_is_ref=0, unmatched=0, bad_marker=0)
    with open(a.tbl, newline="") as fi, open(a.out, "w", newline="") as fo:
        rd = csv.DictReader(fi, delimiter="\t")
        need = ["MarkerName", "Allele1", "Allele2", "Freq1", "MinFreq", "MaxFreq", "Effect", "Direction"]
        missing = [c for c in need if c not in rd.fieldnames]
        if missing:
            sys.exit(f"ERROR: {a.tbl} lacks columns {missing}")
        extra = ["CHROM", "POS", "REF", "ALT", "EFFECT_ALLELE", "OTHER_ALLELE", "EFFECT_ALLELE_IS",
                 "BETA_ALT", "FREQ_ALT", "MINFREQ_ALT", "MAXFREQ_ALT", "DIRECTION_ALT"]
        wr = csv.DictWriter(fo, fieldnames=rd.fieldnames + extra, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for row in rd:
            counts["total"] += 1
            parts = row["MarkerName"].split(":")
            a1 = row["Allele1"].upper()
            a2 = row["Allele2"].upper()
            row.update(EFFECT_ALLELE=a1, OTHER_ALLELE=a2)
            if len(parts) != 4:
                counts["bad_marker"] += 1
                row.update(CHROM="NA", POS="NA", REF="NA", ALT="NA", EFFECT_ALLELE_IS="UNMATCHED",
                           BETA_ALT="NA", FREQ_ALT="NA", MINFREQ_ALT="NA", MAXFREQ_ALT="NA", DIRECTION_ALT="NA")
                wr.writerow(row)
                continue
            chrom, pos, ref, alt = parts[0].replace("chr", ""), parts[1], parts[2].upper(), parts[3].upper()
            row.update(CHROM=chrom, POS=pos, REF=ref, ALT=alt)
            eff, f1, mn, mx = fnum(row["Effect"]), fnum(row["Freq1"]), fnum(row["MinFreq"]), fnum(row["MaxFreq"])
            if a1 == alt and a2 == ref:
                counts["effect_is_alt"] += 1
                row.update(EFFECT_ALLELE_IS="ALT",
                           BETA_ALT=row["Effect"], FREQ_ALT=row["Freq1"],
                           MINFREQ_ALT=row["MinFreq"], MAXFREQ_ALT=row["MaxFreq"],
                           DIRECTION_ALT=row["Direction"])
            elif a1 == ref and a2 == alt:
                counts["effect_is_ref"] += 1
                row.update(EFFECT_ALLELE_IS="REF",
                           BETA_ALT=f"{-eff:.4f}" if eff is not None else "NA",
                           FREQ_ALT=f"{1 - f1:.4f}" if f1 is not None else "NA",
                           # min/max swap under complementation
                           MINFREQ_ALT=f"{1 - mx:.4f}" if mx is not None else "NA",
                           MAXFREQ_ALT=f"{1 - mn:.4f}" if mn is not None else "NA",
                           DIRECTION_ALT=flip_dir(row["Direction"]))
            else:
                counts["unmatched"] += 1
                row.update(EFFECT_ALLELE_IS="UNMATCHED", BETA_ALT="NA", FREQ_ALT="NA",
                           MINFREQ_ALT="NA", MAXFREQ_ALT="NA", DIRECTION_ALT="NA")
            wr.writerow(row)

    with open(a.summary, "w") as fs:
        for k, v in counts.items():
            fs.write(f"{k}\t{v}\n")
    print(f"[annotate] {a.tbl}: " + "  ".join(f"{k}={v}" for k, v in counts.items()))
    if counts["total"] and counts["unmatched"] / counts["total"] > 0.001:
        sys.exit(f"ERROR: {counts['unmatched']} variants whose METAL alleles match neither REF/ALT of MarkerName")


if __name__ == "__main__":
    main()
