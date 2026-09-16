#!/usr/bin/env python3
"""
Figure 1 cell-type panel: MGP score distributions across all 15 bulk cohorts.

One row per cell type; within each row the 15 cohorts are drawn as overlaid
semi-transparent ridges (KDE), so a single compact panel shows both which cell
types vary most and how consistent that spread is across cohorts.

IMPORTANT - what this panel can and cannot say
----------------------------------------------
The MGP estimates in cell_proportions{,_scaled}.csv are *mean-centred within
each cohort* (verified: every cohort has mean 0 for every cell type).  They are
relative marker-gene scores, not compositional proportions - they take negative
values and do not sum to 1 across cell types.

So this panel deliberately does NOT claim anything about absolute cell type
abundance, and does not compare cohort means: those comparisons are not
identifiable from this data.  What it shows is the *distribution* (spread and
shape) of the GWAS phenotype in each cohort, which is exactly the quantity that
matters for the GWAS.  The axis is labelled accordingly.

Reads the same per-cohort files the GWAS used, and writes only to --out-dir.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.size"] = 7
matplotlib.rcParams["axes.linewidth"] = 0.6

COHORTS = [
    ("ROSMAP", "ROSMAP"), ("ROSMAP_array", "ROSMAP_array"),
    ("Mayo", "Mayo"), ("MSBB", "MSBB"),
    ("CMC_MSSM", "CMC_MSSM"), ("CMC_PENN", "CMC_PENN"),
    ("CMC_PITT", "CMC_PITT"), ("GTEx_v10", "GTEx_v10"),
    ("NABEC", "NABEC"), ("NIMH_HBCC_1M", "NIMH_HBCC_1M"),
    ("NIMH_HBCC_h650", "NIMH_HBCC_h650"),
    ("NIMH_HBCC_Omni5M", "NIMH_HBCC_Omni5M"),
    ("GVEX", "GVEX"), ("AMP_AD_Rush", "AMP_AD_Rush"),
    ("AMP_AD_Mayo", "AMP_AD_Mayo"),
]

# Cell types grouped Inhibitory / Excitatory / Non-neuronal, with display names.
GROUPS = [
    ("Inhibitory", "#C2477B", [
        ("LAMP5", "LAMP5"), ("PAX6", "PAX6"), ("VIP", "VIP"),
        ("SST", "SST"), ("PVALB", "PVALB"),
    ]),
    ("Excitatory", "#3B6BC4", [
        ("IT", "IT"), ("L4.IT", "L4 IT"), ("L5.6.IT.Car3", "L5/6 IT Car3"),
        ("L5.ET", "L5 ET"), ("L5.6.NP", "L5/6 NP"),
        ("L6.CT", "L6 CT"), ("L6b", "L6b"),
    ]),
    ("Non-neuronal", "#2E9B6F", [
        ("Astrocyte", "Astrocyte"), ("Oligodendrocyte", "Oligodendrocyte"),
        ("OPC", "OPC"), ("Microglia", "Microglia"),
        ("Endothelial", "Endothelial"), ("Pericyte", "Pericyte"),
        ("VLMC", "VLMC"),
    ]),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geno-root",
                   default="/scratch/zhoux156/results/cohort_gwas_v2/results",
                   help="dir holding <cohort>/cell_proportions_scaled.csv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--xlim", type=float, default=3.2,
                   help="x-axis limit in SD units")
    p.add_argument("--compact", action="store_true",
                   help="narrower/shorter figure for a small Figure 1 slot")
    return p.parse_args()


def load(geno_root):
    """Per-cohort scaled MGP tables, keyed by cohort."""
    out = {}
    for cohort, _ in COHORTS:
        f = Path(geno_root) / cohort / "cell_proportions_scaled.csv"
        if not f.exists():
            print(f"[celltype] WARNING: missing {f}")
            continue
        out[cohort] = pd.read_csv(f, index_col=0)
    if not out:
        raise SystemExit("ERROR: no cell_proportions_scaled.csv found")
    return out


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    data = load(args.geno_root)
    print(f"[celltype] loaded {len(data)} cohorts")

    rows = [(ct, disp, gname, gcol)
            for gname, gcol, cts in GROUPS for ct, disp in cts]

    # A cohort contributes one faint ridge per row; the mean across cohorts is
    # drawn on top as a solid line so the row still reads at small size.
    if args.compact:
        figsize, row_h, lab_fs = (2.5, 3.4), 1.0, 5.5
    else:
        figsize, row_h, lab_fs = (3.2, 4.6), 1.0, 6.5

    fig, ax = plt.subplots(figsize=figsize)
    xs = np.linspace(-args.xlim, args.xlim, 256)

    yticks, ylabels, tickcols = [], [], []
    y = 0.0
    group_spans = {}

    for i, (ct, disp, gname, gcol) in enumerate(rows):
        if i and rows[i - 1][2] != gname:
            y += 0.9                      # gap between groups
        group_spans.setdefault(gname, [y, y, gcol])
        group_spans[gname][1] = y

        # Per-cohort ridges.
        stack = []
        for cohort, _ in COHORTS:
            df = data.get(cohort)
            if df is None or ct not in df.columns:
                continue
            v = df[ct].dropna().values
            if len(v) < 10:
                continue
            try:
                dens = gaussian_kde(v)(xs)
            except Exception:
                continue
            dens = dens / dens.max() * row_h
            # Ridges are drawn DOWNWARD from the baseline (y - dens) because the
            # rows are laid out top-to-bottom below; using invert_yaxis() instead
            # would flip the densities and make every row read as a trough.
            ax.fill_between(xs, y, y - dens, color=gcol, alpha=0.10,
                            linewidth=0, zorder=2)
            stack.append(dens)

        if stack:
            m = np.mean(stack, axis=0)
            ax.plot(xs, y - m, color=gcol, linewidth=0.9, zorder=3)

        yticks.append(y)
        ylabels.append(disp)
        tickcols.append(gcol)
        y += 1.35

    ax.axvline(0, color="#999999", linewidth=0.5, linestyle=(0, (3, 3)),
               zorder=1)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=lab_fs)
    for t, c in zip(ax.get_yticklabels(), tickcols):
        t.set_color(c)
    # Rows already run top-to-bottom (y increases downward as rows are added),
    # so the limits are set reversed rather than calling invert_yaxis().
    ax.set_ylim(y + 0.15, -row_h - 0.35)
    ax.set_xlim(-args.xlim, args.xlim)
    ax.set_xlabel("MGP estimated proportion score\n(SD, within cohort)",
                  fontsize=lab_fs + 0.5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.tick_params(axis="x", length=2.5, width=0.6, pad=1.5, labelsize=lab_fs)

    # Group labels down the right edge.
    for gname, (y0, y1, gcol) in group_spans.items():
        ax.text(args.xlim * 1.04, (y0 + y1) / 2, gname, rotation=270,
                va="center", ha="left", fontsize=lab_fs, color=gcol)

    handles = [
        Line2D([], [], color="#555555", linewidth=0.9,
               label=f"mean of {len(data)} cohorts"),
        matplotlib.patches.Patch(facecolor="#555555", alpha=0.25,
                                 edgecolor="none", label="individual cohort"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=lab_fs - 0.3,
              loc="lower center", bbox_to_anchor=(0.5, 1.0),
              handletextpad=0.4, handlelength=1.2, ncol=2, columnspacing=1.0)

    fig.tight_layout(pad=0.3)
    stem = "figure1_panelE_celltypes"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)

    # Accompanying table: per-cohort SD per cell type, from the UNSCALED file.
    # cell_proportions_scaled.csv is z-scored within cohort, so its SD is 1.0 by
    # construction; the raw file is where the spread is still informative.
    raw = {}
    for cohort, _ in COHORTS:
        f = Path(args.geno_root) / cohort / "cell_proportions.csv"
        if f.exists():
            raw[cohort] = pd.read_csv(f, index_col=0).std()
    if raw:
        sd = pd.DataFrame(raw).T
        sd.index.name = "cohort"
        sd.round(3).to_csv(out_dir / "figure1_panelE_celltype_sd.tsv", sep="\t")

    n = {c: len(d) for c, d in data.items()}
    pd.Series(n, name="n_samples").rename_axis("cohort").to_csv(
        out_dir / "figure1_panelE_sample_counts.tsv", sep="\t")

    print(f"[celltype] wrote {stem}.(svg|pdf|png) to {out_dir}")
    print(f"[celltype] samples per cohort: {n}")


if __name__ == "__main__":
    main()
