#!/usr/bin/env python3
"""
Run METAL on pre-harmonized per-cohort summary stats (no Nextflow / liftover / harmonization).

Expects harmonized files from a prior stratified meta run:
  {harmonized_dir}/{cell_type}/harmonized/{cohort}_{cell_type}_harmonized.raw_p
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--cohorts", nargs="+", required=True)
    parser.add_argument("--harmonized-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--metal-path", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    harmonized_dir = Path(args.harmonized_dir) / args.cell_type / "harmonized"
    if not harmonized_dir.is_dir():
        sys.exit(f"Harmonized directory not found: {harmonized_dir}")

    harm_files = []
    for cohort in args.cohorts:
        path = harmonized_dir / f"{cohort}_{args.cell_type}_harmonized.raw_p"
        if not path.is_file():
            sys.exit(f"Missing harmonized file: {path}")
        harm_files.append(path)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cohort_suffix = "_".join(sorted(args.cohorts))
    prefix = output_dir / f"{args.cell_type}_meta_analysis_{cohort_suffix}"
    script_path = output_dir / f"{args.cell_type}_metal_script.txt"

    process_lines = "\n".join(f"PROCESS {path}" for path in harm_files)
    script_path.write_text(
        "\n".join(
            [
                "SCHEME STDERR",
                "AVERAGEFREQ ON",
                "MINMAXFREQ ON",
                "FLIP OFF",
                "MARKER ID",
                "ALLELE ALLELE0 ALLELE1",
                "FREQ A1FREQ",
                "EFFECT BETA",
                "STDERR SE",
                "PVAL P",
                process_lines,
                f"OUTFILE {prefix} .tbl",
                "ANALYZE HETEROGENEITY",
                "QUIT",
                "",
            ]
        )
    )

    env = os.environ.copy()
    env["PATH"] = f"{args.metal_path}:{env.get('PATH', '')}"
    print(f"Running METAL for {args.cell_type} ({len(harm_files)} cohorts)...", flush=True)
    subprocess.run(["metal", str(script_path)], check=True, env=env)

    tbl_files = sorted(output_dir.glob(f"{args.cell_type}_meta_analysis_{cohort_suffix}*.tbl"))
    if not tbl_files:
        sys.exit("METAL finished but no .tbl output was found")
    print(f"Done: {tbl_files[0]}")


if __name__ == "__main__":
    main()
