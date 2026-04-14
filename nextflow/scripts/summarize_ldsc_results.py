#!/usr/bin/env python3
"""
Summarize LDSC h2 logs into one TSV.
"""

import argparse
import csv
import re
from pathlib import Path


PATTERNS = {
    "observed_h2": re.compile(r"Total Observed scale h2:\s*([^\s]+)\s*\(([^)]+)\)"),
    "lambda_gc": re.compile(r"Lambda GC:\s*([^\s]+)"),
    "mean_chi2": re.compile(r"Mean Chi\^2:\s*([^\s]+)"),
    "intercept": re.compile(r"Intercept:\s*([^\s]+)\s*\(([^)]+)\)"),
    "ratio": re.compile(r"Ratio:\s*([^\s]+)\s*\(([^)]+)\)"),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ldsc-logs", nargs="+", required=True)
    parser.add_argument("--output-tsv", required=True)
    return parser.parse_args()


def parse_log(path):
    text = Path(path).read_text()
    row = {
        "cell_type": Path(path).name.replace(".h2.log", ""),
        "observed_h2": "",
        "observed_h2_se": "",
        "lambda_gc": "",
        "mean_chi2": "",
        "intercept": "",
        "intercept_se": "",
        "ratio": "",
        "ratio_se": "",
    }

    match = PATTERNS["observed_h2"].search(text)
    if match:
        row["observed_h2"], row["observed_h2_se"] = match.groups()

    match = PATTERNS["lambda_gc"].search(text)
    if match:
        row["lambda_gc"] = match.group(1)

    match = PATTERNS["mean_chi2"].search(text)
    if match:
        row["mean_chi2"] = match.group(1)

    match = PATTERNS["intercept"].search(text)
    if match:
        row["intercept"], row["intercept_se"] = match.groups()

    match = PATTERNS["ratio"].search(text)
    if match:
        row["ratio"], row["ratio_se"] = match.groups()

    return row


def main():
    args = parse_args()
    rows = [parse_log(path) for path in sorted(args.ldsc_logs)]

    with open(args.output_tsv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cell_type",
                "observed_h2",
                "observed_h2_se",
                "lambda_gc",
                "mean_chi2",
                "intercept",
                "intercept_se",
                "ratio",
                "ratio_se",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
