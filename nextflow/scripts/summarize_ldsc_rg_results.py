#!/usr/bin/env python3
"""Summarize LDSC genetic-correlation logs produced by run_sn_bulk_meta_ldsc_rg.sh."""

import argparse
import csv
import math
import re
from pathlib import Path


RG_HEADER_RE = re.compile(r"^p1\s+p2\s+rg\s+se\s+z\s+p", re.IGNORECASE)
H2_SECTION_RE = re.compile(
    r"Heritability of phenotype (?P<phenotype>1|2/2)\n[-]+\n"
    r"Total Observed scale h2:\s*(?P<h2>[^\s]+)\s*\((?P<se>[^)]+)\)",
    re.MULTILINE,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--rg-dir", required=True)
    parser.add_argument("--output-tsv", required=True)
    return parser.parse_args()


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def parse_rg_log(path):
    row = {
        "rg": "",
        "rg_se": "",
        "rg_z": "",
        "rg_p": "",
        "h2_sn": "",
        "h2_sn_se": "",
        "h2_bulk": "",
        "h2_bulk_se": "",
        "status": "missing_log",
        "warnings": "",
    }
    if not path.exists():
        return row

    text = path.read_text(errors="replace")
    warnings = []
    for line in text.splitlines():
        lower = line.lower()
        if "warning" in lower or "error" in lower or "failed" in lower:
            warnings.append(line.strip())

    lines = text.splitlines()
    for idx, line in enumerate(lines):
        if RG_HEADER_RE.match(line.strip()) and idx + 1 < len(lines):
            fields = lines[idx + 1].split()
            if len(fields) >= 6:
                row["rg"] = fields[2]
                row["rg_se"] = fields[3]
                row["rg_z"] = fields[4]
                row["rg_p"] = fields[5]
                row["status"] = "ok"
            break

    for match in H2_SECTION_RE.finditer(text):
        phenotype = match.group("phenotype")
        h2 = match.group("h2")
        se = match.group("se")
        if phenotype == "1":
            row["h2_sn"] = h2
            row["h2_sn_se"] = se
        else:
            row["h2_bulk"] = h2
            row["h2_bulk_se"] = se

    if "ERROR computing rg" in text:
        row["status"] = "rg_failed"
    elif row["rg"].upper() == "NA" or row["rg"].lower() == "nan":
        row["status"] = "rg_na"

    flags = []
    rg = safe_float(row["rg"])
    rg_se = safe_float(row["rg_se"])
    h2_sn = safe_float(row["h2_sn"])
    h2_bulk = safe_float(row["h2_bulk"])
    if row["status"] != "ok":
        flags.append(row["status"])
    if math.isfinite(rg) and abs(rg) > 1.25:
        flags.append("rg_outside_expected_range")
    if math.isfinite(rg_se) and rg_se > 0.5:
        flags.append("high_rg_se")
    if math.isfinite(h2_sn) and h2_sn <= 0:
        flags.append("nonpositive_sn_h2")
    if math.isfinite(h2_bulk) and h2_bulk <= 0:
        flags.append("nonpositive_bulk_h2")
    if "h2  out of bounds" in text or "h2's was out of bounds" in text:
        flags.append("h2_out_of_bounds")
    if warnings:
        flags.append("log_warning")

    row["reliability_flags"] = ";".join(dict.fromkeys(flags))
    row["warnings"] = " | ".join(warnings[:10])
    return row


def main():
    args = parse_args()
    rg_dir = Path(args.rg_dir)
    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.manifest, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        manifest_rows = list(reader)

    fieldnames = [
        "sn_cell_type",
        "bulk_cell_type",
        "rg",
        "rg_se",
        "rg_z",
        "rg_p",
        "h2_sn",
        "h2_sn_se",
        "h2_bulk",
        "h2_bulk_se",
        "status",
        "reliability_flags",
        "warnings",
        "log_file",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for manifest_row in manifest_rows:
            sn_ct = manifest_row["sn_cell_type"]
            bulk_ct = manifest_row["bulk_cell_type"]
            log_file = rg_dir / f"{sn_ct}_vs_{bulk_ct}.log"
            parsed = parse_rg_log(log_file)
            out = {
                "sn_cell_type": sn_ct,
                "bulk_cell_type": bulk_ct,
                "log_file": str(log_file),
                **parsed,
            }
            writer.writerow({key: out.get(key, "") for key in fieldnames})


if __name__ == "__main__":
    main()
