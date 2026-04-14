#!/usr/bin/env python3
"""
Prepare LDSC input from one METAL output table.

N strategy:
  For each meta-analysis variant, compute N as the sum of the REGENIE N column
  across the harmonized cohort .raw_p files that contributed that variant.

This keeps the LDSC input aligned with the actual cohort presence after the
pre-METAL harmonization filters, instead of using a single constant N.
"""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meta-tsv", required=True)
    parser.add_argument("--harmonized-files", nargs="+", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--variant-map", default="")
    return parser.parse_args()


def detect_delimiter(line):
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None


def normalize_header(name):
    return "".join(ch for ch in name.lower() if ch.isalnum())


def load_variant_map(path):
    if not path:
        return {}

    with open(path, "r", newline="") as handle:
        first_line = handle.readline().strip()
        if not first_line:
            return {}
        handle.seek(0)
        delimiter = detect_delimiter(first_line)

        if delimiter is None:
            mapping = {}
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                mapping[parts[0]] = parts[1]
            return mapping

        reader = csv.DictReader(handle, delimiter=delimiter)
        normalized = {normalize_header(name): name for name in reader.fieldnames or []}
        source_key = None
        target_key = None

        for candidate in ("variantid", "id", "snp", "markername", "oldid"):
            if candidate in normalized:
                source_key = normalized[candidate]
                break
        for candidate in ("rsid", "newid", "snp", "markername"):
            if candidate in normalized and normalized[candidate] != source_key:
                target_key = normalized[candidate]
                break

        if source_key is None or target_key is None:
            raise ValueError(
                f"Could not infer source/target columns from variant map: {path}"
            )

        mapping = {}
        for row in reader:
            source = (row.get(source_key) or "").strip()
            target = (row.get(target_key) or "").strip()
            if source and target:
                mapping[source] = target
        return mapping


def build_variant_n_map(harmonized_files):
    variant_n = {}
    for path in harmonized_files:
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=" ")
            for row in reader:
                variant_id = row.get("ID", "").strip()
                if not variant_id:
                    continue
                n_value = row.get("N", "").strip()
                if not n_value:
                    continue
                try:
                    variant_n[variant_id] = variant_n.get(variant_id, 0) + int(float(n_value))
                except ValueError:
                    continue
    return variant_n


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    args = parse_args()
    variant_map = load_variant_map(args.variant_map)
    variant_n = build_variant_n_map(args.harmonized_files)

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_tsv).parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "meta_rows_total": 0,
        "rows_written": 0,
        "rows_missing_n": 0,
        "rows_missing_variant_map": 0,
        "rows_invalid": 0,
        "rows_duplicate_snp": 0,
        "variant_map_entries": len(variant_map),
        "n_strategy": "sum_per_variant_from_harmonized_regenie_N",
    }

    seen_snps = set()

    with open(args.meta_tsv, "r", newline="") as meta_handle, open(
        args.output_tsv, "w", newline=""
    ) as out_handle:
        reader = csv.DictReader(meta_handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["SNP", "A1", "A2", "P", "BETA", "N"])

        for row in reader:
            summary["meta_rows_total"] += 1
            variant_id = (row.get("MarkerName") or "").strip()
            if not variant_id:
                summary["rows_invalid"] += 1
                continue

            mapped_id = variant_map.get(variant_id, variant_id)
            if args.variant_map and variant_id not in variant_map:
                summary["rows_missing_variant_map"] += 1
                continue
            if not mapped_id:
                summary["rows_invalid"] += 1
                continue
            if mapped_id in seen_snps:
                summary["rows_duplicate_snp"] += 1
                continue

            n_value = variant_n.get(variant_id)
            if n_value is None or n_value <= 0:
                summary["rows_missing_n"] += 1
                continue

            p_value = safe_float(row.get("P-value"))
            beta = safe_float(row.get("Effect"))
            allele1 = (row.get("Allele1") or "").upper()
            allele2 = (row.get("Allele2") or "").upper()

            if (
                p_value is None
                or beta is None
                or not allele1
                or not allele2
                or p_value <= 0.0
                or p_value > 1.0
            ):
                summary["rows_invalid"] += 1
                continue

            writer.writerow([
                mapped_id,
                allele1,
                allele2,
                f"{p_value:.12g}",
                f"{beta:.12g}",
                n_value,
            ])
            seen_snps.add(mapped_id)
            summary["rows_written"] += 1

    with open(args.summary_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(list(summary.keys()))
        writer.writerow([summary[key] for key in summary.keys()])


if __name__ == "__main__":
    main()
