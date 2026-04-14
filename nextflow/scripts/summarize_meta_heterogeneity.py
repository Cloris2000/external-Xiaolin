#!/usr/bin/env python3
"""
Summarize METAL heterogeneity outputs across cell types.

This script is intentionally streaming-friendly:
- Reads each METAL .tbl file row-by-row
- Computes per-cell-type heterogeneity summary metrics
- Writes a focused top-hit table for downstream MR-MEGA triage
"""

import argparse
import csv
import math
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meta-tbls", nargs="+", required=True, help="METAL .tbl files")
    parser.add_argument("--output-summary-tsv", required=True, help="Per-cell-type summary TSV")
    parser.add_argument("--output-top-hits-tsv", required=True, help="Top-hit heterogeneity TSV")
    parser.add_argument("--gw-p-thresh", type=float, default=5e-8)
    parser.add_argument("--suggestive-p-thresh", type=float, default=1e-5)
    parser.add_argument("--i2-thresh", type=float, default=50.0)
    parser.add_argument("--het-p-thresh", type=float, default=0.05)
    return parser.parse_args()


def to_float(value):
    if value is None:
        return None
    value = str(value).strip()
    if value == "" or value.upper() == "NA":
        return None
    try:
        parsed = float(value)
        if math.isnan(parsed) or math.isinf(parsed):
            return None
        return parsed
    except ValueError:
        return None


def infer_cell_type(path_str):
    name = Path(path_str).name
    suffix = "_meta_analysis_"
    if suffix in name:
        return name.split(suffix, 1)[0]
    return name.replace(".tbl", "")


def high_heterogeneity(i2, het_p, i2_thresh, het_p_thresh):
    i2_flag = i2 is not None and i2 >= i2_thresh
    het_p_flag = het_p is not None and het_p < het_p_thresh
    return i2_flag or het_p_flag


def summarize_one_file(path, args):
    cell_type = infer_cell_type(path)
    summary = {
        "cell_type": cell_type,
        "n_variants": 0,
        "n_with_het": 0,
        "n_i2_ge_thresh": 0,
        "n_hetp_lt_thresh": 0,
        "n_gw": 0,
        "n_gw_high_het": 0,
        "n_suggestive": 0,
        "n_suggestive_high_het": 0,
    }
    top_rows = []

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            summary["n_variants"] += 1

            marker = row.get("MarkerName", "")
            p = to_float(row.get("P-value"))
            i2 = to_float(row.get("HetISq"))
            het_p = to_float(row.get("HetPVal"))

            if i2 is not None or het_p is not None:
                summary["n_with_het"] += 1
            if i2 is not None and i2 >= args.i2_thresh:
                summary["n_i2_ge_thresh"] += 1
            if het_p is not None and het_p < args.het_p_thresh:
                summary["n_hetp_lt_thresh"] += 1

            if p is None:
                continue

            is_gw = p <= args.gw_p_thresh
            is_suggestive = p <= args.suggestive_p_thresh
            is_high_het = high_heterogeneity(i2, het_p, args.i2_thresh, args.het_p_thresh)

            if is_gw:
                summary["n_gw"] += 1
                if is_high_het:
                    summary["n_gw_high_het"] += 1

            if is_suggestive:
                summary["n_suggestive"] += 1
                if is_high_het:
                    summary["n_suggestive_high_het"] += 1
                    top_rows.append(
                        {
                            "cell_type": cell_type,
                            "MarkerName": marker,
                            "P_value": p,
                            "HetISq": "" if i2 is None else i2,
                            "HetPVal": "" if het_p is None else het_p,
                            "is_gw_sig": int(is_gw),
                            "high_heterogeneity_flag": 1,
                        }
                    )

    n_het = summary["n_with_het"]
    summary["pct_i2_ge_thresh"] = (
        round(100.0 * summary["n_i2_ge_thresh"] / n_het, 4) if n_het > 0 else 0.0
    )
    summary["pct_hetp_lt_thresh"] = (
        round(100.0 * summary["n_hetp_lt_thresh"] / n_het, 4) if n_het > 0 else 0.0
    )
    summary["pct_suggestive_high_het"] = (
        round(100.0 * summary["n_suggestive_high_het"] / summary["n_suggestive"], 4)
        if summary["n_suggestive"] > 0
        else 0.0
    )
    summary["mr_mega_priority"] = (
        "yes" if summary["n_gw_high_het"] > 0 or summary["n_suggestive_high_het"] >= 10 else "no"
    )

    return summary, top_rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    summary_rows = []
    top_rows = []

    cleaned_paths = [p for p in sorted(args.meta_tbls) if not Path(p).name.endswith("1.tbl")]
    for path in cleaned_paths:
        summary, top_hits = summarize_one_file(path, args)
        summary_rows.append(summary)
        top_rows.extend(top_hits)

    summary_fields = [
        "cell_type",
        "n_variants",
        "n_with_het",
        "n_i2_ge_thresh",
        "pct_i2_ge_thresh",
        "n_hetp_lt_thresh",
        "pct_hetp_lt_thresh",
        "n_gw",
        "n_gw_high_het",
        "n_suggestive",
        "n_suggestive_high_het",
        "pct_suggestive_high_het",
        "mr_mega_priority",
    ]
    top_fields = [
        "cell_type",
        "MarkerName",
        "P_value",
        "HetISq",
        "HetPVal",
        "is_gw_sig",
        "high_heterogeneity_flag",
    ]

    write_tsv(args.output_summary_tsv, summary_rows, summary_fields)
    write_tsv(args.output_top_hits_tsv, top_rows, top_fields)


if __name__ == "__main__":
    main()
