#!/usr/bin/env python3
"""
Batch-generate Manhattan and QQ plots from REGENIE *.raw_p files.

Usage:
  python scripts/plot_gwas_from_filelist.py \
      --file_list results/MSBB_sn/gwas_raw_p_file_list.txt \
      --output_dir results/MSBB_sn/plots_gwas
"""

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Manhattan + QQ for GWAS raw_p files")
    parser.add_argument(
        "--file_list",
        required=True,
        help="Text file with one absolute/relative path to *.regenie.raw_p per line",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to write PNGs and summary table",
    )
    parser.add_argument(
        "--max_points",
        type=int,
        default=0,
        help="Maximum variants per file to plot (<=0 means no downsampling; default: 0)",
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=42,
        help="Random seed for downsampling",
    )
    return parser.parse_args()


def extract_cell_type(path: Path) -> str:
    name = path.name
    prefix = "MSBB_"
    suffix = "_step2.regenie.raw_p"
    if name.startswith(prefix) and name.endswith(suffix):
        return name[len(prefix) : -len(suffix)]
    return path.stem


def load_gwas_table(path: Path) -> pd.DataFrame:
    # Regenie raw_p is space-delimited with columns including CHROM, GENPOS, P.
    df = pd.read_csv(
        path,
        sep=r"\s+",
        usecols=["CHROM", "GENPOS", "P"],
        dtype={"CHROM": str, "GENPOS": np.int64, "P": float},
        engine="c",
    )
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["CHROM", "GENPOS", "P"])
    df = df[(df["P"] > 0) & (df["P"] <= 1)]
    return df


def make_manhattan(df: pd.DataFrame, out_png: Path, title: str) -> None:
    chrom_map = {str(i): i for i in range(1, 23)}
    df = df[df["CHROM"].isin(chrom_map)].copy()
    if df.empty:
        return
    df["chrom_num"] = df["CHROM"].map(chrom_map).astype(int)
    df = df.sort_values(["chrom_num", "GENPOS"])
    df["minus_log10_p"] = -np.log10(df["P"])

    # Build cumulative x-axis by chromosome.
    x_offset = 0
    xticks = []
    xtick_labels = []
    plot_chunks = []
    for chrom in sorted(df["chrom_num"].unique()):
        c = df[df["chrom_num"] == chrom].copy()
        c["x"] = c["GENPOS"] + x_offset
        plot_chunks.append(c)
        xticks.append(c["x"].median())
        xtick_labels.append(str(chrom))
        x_offset = c["x"].max() + 1_000_000
    plot_df = pd.concat(plot_chunks, ignore_index=True)

    plt.figure(figsize=(14, 5))
    for chrom in sorted(plot_df["chrom_num"].unique()):
        c = plot_df[plot_df["chrom_num"] == chrom]
        color = "#4C78A8" if chrom % 2 else "#F58518"
        plt.scatter(c["x"], c["minus_log10_p"], s=2, c=color, alpha=0.65, linewidths=0)

    plt.axhline(-math.log10(5e-8), color="red", linestyle="--", linewidth=1, alpha=0.8)
    plt.xticks(xticks, xtick_labels, fontsize=8)
    plt.ylabel("-log10(P)")
    plt.xlabel("Chromosome")
    plt.title(f"Manhattan: {title}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def make_qq(df: pd.DataFrame, out_png: Path, title: str) -> None:
    pvals = np.sort(df["P"].values)
    n = len(pvals)
    if n == 0:
        return
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(pvals)

    m = max(expected.max(), observed.max())
    plt.figure(figsize=(5.5, 5.5))
    plt.scatter(expected, observed, s=2, alpha=0.6, color="#4C78A8", linewidths=0)
    plt.plot([0, m], [0, m], color="red", linestyle="--", linewidth=1)
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.title(f"QQ: {title}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.random_seed)

    file_list_path = Path(args.file_list)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with file_list_path.open("r") as f:
        files = [Path(line.strip()) for line in f if line.strip()]

    summary_rows = []
    for fp in files:
        if not fp.exists():
            summary_rows.append({"file": str(fp), "status": "missing"})
            continue

        cell = extract_cell_type(fp)
        try:
            df = load_gwas_table(fp)
            total_n = len(df)
            if args.max_points > 0 and total_n > args.max_points:
                keep_idx = rng.choice(total_n, size=args.max_points, replace=False)
                df = df.iloc[np.sort(keep_idx)].copy()
            used_n = len(df)

            manhattan_out = output_dir / f"{cell}.manhattan.png"
            qq_out = output_dir / f"{cell}.qq.png"

            make_manhattan(df, manhattan_out, cell)
            make_qq(df, qq_out, cell)

            summary_rows.append(
                {
                    "file": str(fp),
                    "cell_type": cell,
                    "status": "ok",
                    "variants_total": total_n,
                    "variants_plotted": used_n,
                    "manhattan_png": str(manhattan_out),
                    "qq_png": str(qq_out),
                }
            )
            print(f"[ok] {cell}: {used_n:,} plotted")
        except Exception as e:  # noqa: BLE001
            summary_rows.append({"file": str(fp), "cell_type": cell, "status": f"error: {e}"})
            print(f"[error] {cell}: {e}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / "plot_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nWrote summary: {summary_path}")


if __name__ == "__main__":
    main()
