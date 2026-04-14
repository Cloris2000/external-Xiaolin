#!/usr/bin/env python3

import argparse
import csv
import glob
import json
import os
from pathlib import Path


def load_metadata(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = {row["cohort"]: row for row in reader}
    return rows


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a cohort QC audit matrix for standalone meta-analysis."
    )
    parser.add_argument("--metadata", required=True, help="Metadata TSV")
    parser.add_argument("--base-dir", required=True, help="Nextflow project base dir")
    parser.add_argument("--cohorts-json", required=True, help="JSON list of cohorts")
    parser.add_argument("--output-tsv", required=True, help="Output TSV path")
    parser.add_argument("--output-md", required=True, help="Output markdown path")
    return parser.parse_args()


def first_header(path_pattern):
    matches = sorted(glob.glob(path_pattern))
    if not matches:
        return [], 0

    with open(matches[0]) as handle:
        header = handle.readline().strip().split()
    return header, len(matches)


def main():
    args = parse_args()
    metadata = load_metadata(args.metadata)
    cohorts = json.loads(args.cohorts_json)

    rows = []
    builds = set()
    for cohort in cohorts:
        meta = metadata.get(cohort, {})
        raw_p_pattern = os.path.join(
            args.base_dir, "results", cohort, "regenie_step2", "*.regenie.raw_p"
        )
        header, raw_p_count = first_header(raw_p_pattern)
        header_set = set(header)
        info_status = "present" if "INFO" in header_set else "missing"
        freq_status = "present" if "A1FREQ" in header_set else "missing"
        builds.add(meta.get("build", ""))

        rows.append(
            {
                "cohort": cohort,
                "build": meta.get("build", ""),
                "ancestry_policy": meta.get("ancestry_policy", ""),
                "mixed_ancestry": meta.get("mixed_ancestry", ""),
                "imputation_metric": meta.get("imputation_metric", ""),
                "info_threshold": meta.get("info_threshold", ""),
                "maf_threshold": meta.get("maf_threshold", ""),
                "tissue_filter": meta.get("tissue_filter", ""),
                "phenotype_normalization": meta.get("phenotype_normalization", ""),
                "deconvolution": meta.get("deconvolution", ""),
                "pc_covariates": meta.get("pc_covariates", ""),
                "clinical_covariates": meta.get("clinical_covariates", ""),
                "technical_covariates": meta.get("technical_covariates", ""),
                "low_qc_sensitivity_exclude": meta.get("low_qc_sensitivity_exclude", ""),
                "raw_p_files_found": raw_p_count,
                "raw_p_complete_19": "yes" if raw_p_count == 19 else "no",
                "info_column_in_sumstats": info_status,
                "freq_column_in_sumstats": freq_status,
                "notes": meta.get("notes", ""),
            }
        )

    fieldnames = list(rows[0].keys()) if rows else []
    with open(args.output_tsv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    build_note = "single_build" if len({b for b in builds if b}) <= 1 else "mixed_builds_detected"
    mixed_count = sum(1 for row in rows if row["mixed_ancestry"] == "yes")
    low_qc_count = sum(
        1 for row in rows if row["low_qc_sensitivity_exclude"].lower() == "yes"
    )
    incomplete_count = sum(1 for row in rows if row["raw_p_complete_19"] != "yes")

    md_lines = [
        "# Meta Cohort QC Audit",
        "",
        f"- Cohorts audited: {len(rows)}",
        f"- Build status: {build_note}",
        f"- Documented mixed-ancestry cohorts: {mixed_count}",
        f"- Low-QC sensitivity exclusions flagged: {low_qc_count}",
        f"- Incomplete raw_p cohorts: {incomplete_count}",
        "",
        "## Cohort Table",
        "",
        "| Cohort | Build | Mixed ancestry | INFO threshold | MAF | raw_p files | INFO column | Low-QC sensitivity exclude |",
        "|---|---|---|---:|---:|---:|---|---|",
    ]

    for row in rows:
        md_lines.append(
            "| {cohort} | {build} | {mixed_ancestry} | {info_threshold} | {maf_threshold} | {raw_p_files_found} | {info_column_in_sumstats} | {low_qc_sensitivity_exclude} |".format(
                **row
            )
        )

    Path(args.output_md).write_text("\n".join(md_lines) + "\n")


if __name__ == "__main__":
    main()
