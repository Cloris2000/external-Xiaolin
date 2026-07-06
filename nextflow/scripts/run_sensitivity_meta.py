#!/usr/bin/env python3
"""
Sensitivity meta-analysis runner.

Wraps standalone_meta_analysis.nf to execute a defined suite of
sensitivity analyses for multi-ancestry GWAS meta-analysis:

  1. leave_one_cohort_out  — one run per cohort, excluding that cohort
  2. exclude_mixed         — exclude all documented mixed-ancestry cohorts
  3. exclude_low_qc        — exclude cohorts flagged low_qc_sensitivity_exclude=yes
  4. hbcc_platforms_only   — meta-analysis of the three HBCC platforms only

Each sensitivity run:
  - writes to results/meta_sensitivity/<label>/
  - uses standalone_meta_analysis.nf with a generated per-run config
  - can be submitted to SLURM via sbatch or run sequentially

Usage:
  python3 scripts/run_sensitivity_meta.py \
      --base-dir /path/to/nextflow \
      --metadata docs/meta_cohort_metadata.tsv \
      --cohorts ROSMAP Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx NABEC \
                NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M GVEX \
      --metal-path /path/to/METAL/executables \
      --submit        # submit each run as sbatch job
      --dry-run       # print commands only, do not execute
"""

import argparse
import csv
import os
import subprocess
import sys
import textwrap
from pathlib import Path


ALL_CELL_TYPES = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.ET", "L5.6.IT.Car3",
    "L5.6.NP", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC",
    "Oligodendrocyte", "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-dir", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--cohorts", nargs="+", required=True)
    p.add_argument("--metal-path", required=True)
    p.add_argument("--submit", action="store_true", help="Submit as sbatch jobs")
    p.add_argument("--dry-run", action="store_true", help="Print commands, do not run")
    p.add_argument("--work-dir-base", default=None)
    return p.parse_args()


def load_metadata(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["cohort"]: row for row in reader}


def cohorts_groovy_list(cohorts):
    items = ", ".join(f'"{c}"' for c in cohorts)
    return f"[{items}]"


def write_config(config_path, base_dir, cohorts, output_dir, metal_path):
    """Write a per-sensitivity-run Nextflow config that overrides the cohort list."""
    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    groovy_list = cohorts_groovy_list(cohorts)
    config_path.write_text(
        textwrap.dedent(f"""\
        /*
         * Auto-generated sensitivity config — do not edit by hand.
         */
        includeConfig '{base_dir}/nextflow.config.standalone_meta'

        params {{
            base_dir       = "{base_dir}"
            cohorts        = {groovy_list}
            cohort_include = {groovy_list}
            output_dir     = "{output_dir}"
            metal_path     = "{metal_path}"
            harmonized_output_dir = "{output_dir}/harmonized"
            audit_output_dir      = "{output_dir}/audit"
            sensitivity_output_dir = "{output_dir}"
        }}
        """)
    )


def write_sbatch(sbatch_path, label, base_dir, config_path, work_dir, log_dir):
    sbatch_path = Path(sbatch_path)
    sbatch_path.parent.mkdir(parents=True, exist_ok=True)
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    sbatch_path.write_text(
        textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=meta_sens_{label}
        #SBATCH --partition=medium
        #SBATCH --time=8:00:00
        #SBATCH --cpus-per-task=2
        #SBATCH --mem=8G
        #SBATCH --output={log_dir}/{label}_%j.out
        #SBATCH --error={log_dir}/{label}_%j.err

        cd {base_dir}
        nextflow run {base_dir}/standalone_meta_analysis.nf \\
            -c {config_path} \\
            -work-dir {work_dir} \\
            -resume
        """)
    )


def build_sensitivity_suite(args, metadata):
    """Return list of (label, cohort_list) pairs for the sensitivity suite."""
    cohorts = args.cohorts
    meta = metadata

    suite = []

    # 1. Leave-one-cohort-out
    for left_out in cohorts:
        remaining = [c for c in cohorts if c != left_out]
        if len(remaining) >= 2:
            suite.append((f"loco_{left_out}", remaining))

    # 2. Exclude mixed-ancestry cohorts
    mixed = [c for c in cohorts if meta.get(c, {}).get("mixed_ancestry", "no") == "yes"]
    if mixed:
        excl_mixed = [c for c in cohorts if c not in mixed]
        if len(excl_mixed) >= 2:
            suite.append(("exclude_mixed_ancestry", excl_mixed))

    # 3. Exclude low-QC cohorts
    low_qc = [c for c in cohorts if meta.get(c, {}).get("low_qc_sensitivity_exclude", "no").lower() == "yes"]
    if low_qc:
        excl_low_qc = [c for c in cohorts if c not in low_qc]
        if len(excl_low_qc) >= 2:
            suite.append(("exclude_low_qc", excl_low_qc))

    # 4. HBCC platforms only
    hbcc = [c for c in cohorts if "HBCC" in c]
    if len(hbcc) >= 2:
        suite.append(("hbcc_platforms_only", hbcc))

    # 5. European-ancestry cohorts only (exclude documented mixed)
    eur_only = [c for c in cohorts if meta.get(c, {}).get("mixed_ancestry", "no") != "yes"]
    if len(eur_only) >= 2 and set(eur_only) != set(cohorts):
        suite.append(("european_only_approximate", eur_only))

    return suite


def main():
    args = parse_args()
    metadata = load_metadata(args.metadata)
    base_dir = os.path.abspath(args.base_dir)
    work_dir_base = args.work_dir_base or os.path.join(base_dir, "work", "meta_sensitivity")
    sensitivity_base = os.path.join(base_dir, "results", "meta_sensitivity")
    sens_config_dir = os.path.join(base_dir, "sensitivity_configs")
    sens_log_dir = os.path.join(base_dir, "logs", "meta_sensitivity")

    suite = build_sensitivity_suite(args, metadata)
    print(f"Sensitivity suite: {len(suite)} analyses", file=sys.stderr)

    for label, cohort_list in suite:
        output_dir = os.path.join(sensitivity_base, label)
        config_path = os.path.join(sens_config_dir, f"nextflow.config.sensitivity.{label}")
        sbatch_path = os.path.join(base_dir, f"run_sensitivity_{label}.sbatch")
        work_dir = os.path.join(work_dir_base, label)

        write_config(config_path, base_dir, cohort_list, output_dir, args.metal_path)
        write_sbatch(sbatch_path, label, base_dir, config_path, work_dir, sens_log_dir)

        cmd = ["sbatch", sbatch_path] if args.submit else [
            "nextflow", "run", os.path.join(base_dir, "standalone_meta_analysis.nf"),
            "-c", config_path, "-work-dir", work_dir, "-resume",
        ]
        print(f"\n[{label}] cohorts ({len(cohort_list)}): {', '.join(cohort_list)}")
        print(f"  config  : {config_path}")
        print(f"  output  : {output_dir}")
        print(f"  command : {' '.join(cmd)}")

        if args.dry_run:
            print("  (dry-run: not executing)")
            continue

        result = subprocess.run(cmd, cwd=base_dir)
        if result.returncode != 0:
            print(f"  WARNING: {label} exited with code {result.returncode}", file=sys.stderr)

    print("\nDone. Generated configs and sbatch scripts for all sensitivity analyses.")
    print(f"Configs: {sens_config_dir}/")
    print(f"Sbatch scripts: {base_dir}/run_sensitivity_*.sbatch")
    print(f"To submit all: ls {base_dir}/run_sensitivity_*.sbatch | xargs -n1 sbatch")


if __name__ == "__main__":
    main()
