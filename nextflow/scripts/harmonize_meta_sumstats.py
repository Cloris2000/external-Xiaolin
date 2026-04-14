#!/usr/bin/env python3
"""
Pre-METAL harmonization of per-cohort REGENIE summary statistics.

Inputs  : one or more <cohort>_<celltype>_step2.regenie.raw_p files
Outputs : harmonized per-cohort files ready for METAL
          one drop-log TSV recording every excluded variant with reason
          one per-celltype summary TSV

Harmonization steps applied (in order):
  1. Normalise variant IDs to chr:pos:ref:alt
  2. Remove strand-ambiguous A/T and C/G SNPs
  3. Apply minimum MAF filter (any cohort failing drops the variant)
  4. Apply minimum presence-in-cohorts filter (>= min_present_cohorts)
  5. Flag (and optionally filter) large cross-cohort AF discrepancies
  6. Log all exclusions with reason

Memory strategy: two-pass.
  Pass 1: read every cohort once, store only (ref, alt, freq) per variant.
           This index is O(unique_variants) not O(total_rows).
  Pass 2: stream each cohort file again, writing only passing variants.
           Only one cohort's rows are in memory at a time.

Columns expected in .raw_p files (REGENIE output + awk-appended P column):
  CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA P
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


STRAND_AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-files", nargs="+", required=True,
                   help="Per-cohort .raw_p files for a single cell type")
    p.add_argument("--cohort-names", required=True,
                   help="JSON array of cohort names matching input files")
    p.add_argument("--cell-type", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--drop-log", required=True)
    p.add_argument("--summary-tsv", required=True)
    p.add_argument("--drop-strand-ambiguous", default="true")
    p.add_argument("--min-meta-maf", type=float, default=0.05)
    p.add_argument("--min-present-cohorts", type=int, default=2)
    p.add_argument("--af-delta-report-threshold", type=float, default=0.20)
    p.add_argument("--af-delta-filter-threshold", type=float, default=None)
    return p.parse_args()


def norm_id(chrom, pos, ref, alt):
    chrom = str(chrom).lstrip("chr").lstrip("0") or "0"
    return f"chr{chrom}:{pos}:{ref.upper()}:{alt.upper()}"


def is_strand_ambiguous(ref, alt):
    return frozenset((ref.upper(), alt.upper())) in STRAND_AMBIGUOUS


def pass1_build_index(input_files, cohort_names):
    """
    Pass 1: read each file once, build a compact index.
    Returns:
      presence  : dict  vid -> set of cohort names
      alleles   : dict  vid -> (ref, alt)  from first cohort seen
      af_by_coh : dict  vid -> dict cohort -> float freq
    """
    presence  = defaultdict(set)
    alleles   = {}
    af_by_coh = defaultdict(dict)

    for path, cohort in zip(input_files, cohort_names):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter=" ")
            for row in reader:
                raw_id = row.get("ID", "")
                if raw_id.count(":") >= 3:
                    vid = raw_id
                else:
                    vid = norm_id(row["CHROM"], row["GENPOS"],
                                  row["ALLELE0"], row["ALLELE1"])
                presence[vid].add(cohort)
                if vid not in alleles:
                    alleles[vid] = (row["ALLELE0"], row["ALLELE1"])
                try:
                    af_by_coh[vid][cohort] = float(row["A1FREQ"])
                except (KeyError, ValueError):
                    pass

    return presence, alleles, af_by_coh


def build_pass_set(presence, alleles, af_by_coh, args, drop_strand_ambiguous):
    """
    Evaluate every variant and return:
      variants_pass    : set of vids to keep
      variants_flag_af : set of vids to annotate with AF discrepancy
      drop_log_rows    : list of dicts for the drop log
      summary counts
    """
    variants_pass    = set()
    variants_flag_af = set()
    drop_log_rows    = []
    counts = dict(
        dropped_strand_ambiguous=0,
        dropped_min_presence=0,
        dropped_min_maf=0,
        dropped_af_discrepancy=0,
        flagged_af_discrepancy=0,
    )

    for vid, cohorts_present in presence.items():
        ref, alt = alleles[vid]

        # Step 1: strand-ambiguous
        if drop_strand_ambiguous and is_strand_ambiguous(ref, alt):
            for c in cohorts_present:
                drop_log_rows.append((args.cell_type, vid, c, "strand_ambiguous"))
            counts["dropped_strand_ambiguous"] += 1
            continue

        # Step 2: min presence
        if len(cohorts_present) < args.min_present_cohorts:
            for c in cohorts_present:
                drop_log_rows.append((
                    args.cell_type, vid, c,
                    f"present_in_only_{len(cohorts_present)}_cohort(s)"
                ))
            counts["dropped_min_presence"] += 1
            continue

        # Step 3: MAF
        maf_fail = False
        for c in cohorts_present:
            freq = af_by_coh[vid].get(c)
            if freq is None:
                continue
            maf = min(freq, 1.0 - freq)
            if maf < args.min_meta_maf:
                drop_log_rows.append((
                    args.cell_type, vid, c,
                    f"maf_{maf:.4f}_below_{args.min_meta_maf}"
                ))
                maf_fail = True
        if maf_fail:
            counts["dropped_min_maf"] += 1
            continue

        # Step 4: AF discrepancy
        freqs = [v for v in af_by_coh[vid].values()]
        if len(freqs) >= 2:
            af_delta = max(freqs) - min(freqs)
            if af_delta > args.af_delta_report_threshold:
                variants_flag_af.add(vid)
                counts["flagged_af_discrepancy"] += 1
                for c in cohorts_present:
                    drop_log_rows.append((
                        args.cell_type, vid, c,
                        f"af_discrepancy_flagged_delta_{af_delta:.3f}"
                    ))
            if (args.af_delta_filter_threshold is not None
                    and af_delta > args.af_delta_filter_threshold):
                for c in cohorts_present:
                    drop_log_rows.append((
                        args.cell_type, vid, c,
                        f"af_discrepancy_filtered_delta_{af_delta:.3f}"
                    ))
                counts["dropped_af_discrepancy"] += 1
                continue

        variants_pass.add(vid)

    return variants_pass, variants_flag_af, drop_log_rows, counts


def pass2_write_harmonized(input_files, cohort_names, cell_type,
                           output_dir, variants_pass, variants_flag_af):
    """
    Pass 2: stream each cohort file and write only passing variants.
    One file open at a time — memory stays low.
    """
    os.makedirs(output_dir, exist_ok=True)
    written = 0
    for path, cohort in zip(input_files, cohort_names):
        out_path = os.path.join(output_dir, f"{cohort}_{cell_type}_harmonized.raw_p")
        cohort_written = 0
        with open(path) as fh_in, open(out_path, "w", newline="") as fh_out:
            reader = csv.DictReader(fh_in, delimiter=" ")
            fieldnames = reader.fieldnames
            writer = csv.DictWriter(
                fh_out, fieldnames=fieldnames, delimiter=" ",
                extrasaction="ignore"
            )
            writer.writeheader()
            for row in reader:
                raw_id = row.get("ID", "")
                if raw_id.count(":") >= 3:
                    vid = raw_id
                else:
                    vid = norm_id(row["CHROM"], row["GENPOS"],
                                  row["ALLELE0"], row["ALLELE1"])
                if vid not in variants_pass:
                    continue
                if vid in variants_flag_af:
                    row["EXTRA"] = "AF_DISCREPANCY_FLAGGED"
                writer.writerow(row)
                cohort_written += 1
        written += cohort_written
    return written


def main():
    args = parse_args()
    drop_strand_ambiguous = args.drop_strand_ambiguous.lower() in ("true", "1", "yes")
    cohort_names = json.loads(args.cohort_names)

    assert len(args.input_files) == len(cohort_names), (
        f"Number of input files ({len(args.input_files)}) must equal "
        f"number of cohort names ({len(cohort_names)})"
    )

    print(f"[harmonize] cell_type={args.cell_type}  cohorts={len(cohort_names)}", file=sys.stderr)

    # Pass 1 — build compact index only
    print("  pass 1: building variant index...", file=sys.stderr)
    presence, alleles, af_by_coh = pass1_build_index(args.input_files, cohort_names)
    print(f"  unique variants in union: {len(presence)}", file=sys.stderr)

    # Evaluate filters
    variants_pass, variants_flag_af, drop_log_rows, counts = build_pass_set(
        presence, alleles, af_by_coh, args, drop_strand_ambiguous
    )
    print(f"  variants passed          : {len(variants_pass)}", file=sys.stderr)
    print(f"  dropped (strand ambig)   : {counts['dropped_strand_ambiguous']}", file=sys.stderr)
    print(f"  dropped (presence<{args.min_present_cohorts})     : {counts['dropped_min_presence']}", file=sys.stderr)
    print(f"  dropped (MAF<{args.min_meta_maf})      : {counts['dropped_min_maf']}", file=sys.stderr)
    print(f"  dropped (AF hard filter) : {counts['dropped_af_discrepancy']}", file=sys.stderr)
    print(f"  flagged (AF discrepancy) : {counts['flagged_af_discrepancy']}", file=sys.stderr)

    # Free the index before pass 2 to save memory
    del presence, alleles, af_by_coh

    # Pass 2 — stream each file, write only passing variants
    print("  pass 2: writing harmonized files...", file=sys.stderr)
    pass2_write_harmonized(
        args.input_files, cohort_names, args.cell_type,
        args.output_dir, variants_pass, variants_flag_af
    )

    # Write drop log
    Path(args.drop_log).parent.mkdir(parents=True, exist_ok=True)
    with open(args.drop_log, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["cell_type", "variant_id", "cohort", "reason"])
        writer.writerows(drop_log_rows)

    # Write summary
    summary = {
        "cell_type":               args.cell_type,
        "cohorts":                 len(cohort_names),
        "total_variants_union":    len(variants_pass) + sum(counts.values()),
        "dropped_strand_ambiguous": counts["dropped_strand_ambiguous"],
        "dropped_min_maf":          counts["dropped_min_maf"],
        "dropped_min_presence":     counts["dropped_min_presence"],
        "dropped_af_discrepancy":   counts["dropped_af_discrepancy"],
        "flagged_af_discrepancy":   counts["flagged_af_discrepancy"],
        "variants_passed":          len(variants_pass),
    }
    with open(args.summary_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)

    print("  done.", file=sys.stderr)


if __name__ == "__main__":
    main()
