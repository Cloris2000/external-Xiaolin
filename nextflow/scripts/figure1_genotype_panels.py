#!/usr/bin/env python3
"""
Figure 1 genotype panels A, B and D for the 15-cohort meta-analysis.

Panels (each also written as a standalone SVG for assembly in Inkscape):
  A  PC1 vs PC2, cohort samples coloured by inferred genetic ancestry, on a
     grey 1000 Genomes backdrop.
  B  PC1 vs PC2, cohort samples coloured by cohort.
  D  Per-cohort ancestry composition, stacked to 100%.

Input is figure1_pca_samples.tsv from figure1_genotype_pca.sh.

Text is kept as real text (no path conversion) so labels stay editable in
Inkscape; fonts are left at matplotlib defaults for the same reason.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

# Keep text as text in the SVG so Inkscape can edit it.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.size"] = 7
matplotlib.rcParams["axes.linewidth"] = 0.6

# Superpopulation palette (colour-blind safe, distinct in greyscale order).
ANC_COLORS = {
    "EUR": "#3B4CC0",
    "AFR": "#6FA8DC",
    "AMR": "#2BA88C",
    "SAS": "#D1495B",
    "EAS": "#A4348E",
    "UNCERTAIN": "#B0B0B0",
}
ANC_ORDER = ["EUR", "AFR", "AMR", "SAS", "EAS", "UNCERTAIN"]

COHORT_ORDER = [
    "ROSMAP", "ROSMAP_array", "Mayo", "MSBB",
    "CMC_MSSM", "CMC_PENN", "CMC_PITT",
    "GTEx_v10", "NABEC",
    "NIMH_HBCC_1M", "NIMH_HBCC_h650", "NIMH_HBCC_Omni5M",
    "GVEX", "AMP_AD_Rush", "AMP_AD_Mayo",
]

# 15 distinguishable cohort colours (tab20 minus the near-greys).
COHORT_COLORS = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
    "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
    "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#17becf",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True,
                   help="figure1_pca_samples.tsv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--point-size", type=float, default=7.0)
    p.add_argument("--ref-point-size", type=float, default=2.0)
    return p.parse_args()


def pc_label(df, pc):
    col = f"{pc}_pct_var"
    if col in df.columns and df[col].notna().any():
        return f"{pc} ({df[col].dropna().iloc[0]:.1f}%)"
    return pc


def scatter_backdrop(ax, ref):
    ax.scatter(ref.PC1, ref.PC2, s=1.2, c="#D9D9D9",
               linewidths=0, rasterized=True, zorder=1)


def style_axes(ax, df):
    ax.set_xlabel(pc_label(df, "PC1"))
    ax.set_ylabel(pc_label(df, "PC2"))
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, width=0.6, pad=1.5)


def panel_a(df, out_dir, ps, rps):
    """PC1/PC2 coloured by inferred genetic ancestry."""
    ref = df[df.dataset == "reference"]
    coh = df[df.dataset == "cohort"]

    fig, ax = plt.subplots(figsize=(3.0, 2.6))
    scatter_backdrop(ax, ref)

    present = [a for a in ANC_ORDER if a in set(coh.genetic_ancestry)]
    for anc in present:
        s = coh[coh.genetic_ancestry == anc]
        ax.scatter(s.PC1, s.PC2, s=ps, c=ANC_COLORS[anc],
                   linewidths=0.2, edgecolors="white", label=f"{anc} ({len(s)})",
                   zorder=3 if anc != "UNCERTAIN" else 2)

    style_axes(ax, df)
    leg = ax.legend(frameon=False, loc="upper left", handletextpad=0.3,
                    borderaxespad=0.2, labelspacing=0.25, markerscale=1.4)
    leg.set_title(None)
    fig.tight_layout(pad=0.4)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(Path(out_dir) / f"figure1_panelA_ancestry.{ext}",
                    dpi=400, bbox_inches="tight")
    plt.close(fig)


def panel_b(df, out_dir, ps):
    """PC1/PC2 coloured by cohort."""
    ref = df[df.dataset == "reference"]
    coh = df[df.dataset == "cohort"]

    fig, ax = plt.subplots(figsize=(3.0, 2.6))
    scatter_backdrop(ax, ref)

    order = [c for c in COHORT_ORDER if c in set(coh.cohort)]
    cmap = dict(zip(order, COHORT_COLORS))
    # Colour assignment follows COHORT_ORDER (stable across runs), but the
    # DRAW order is largest-cohort-first so a big cohort cannot bury the
    # smaller ones underneath it.  Points are semi-transparent for the same
    # reason: in the dense EUR cluster every cohort overlaps every other.
    sizes = coh.cohort.value_counts()
    draw_order = sorted(order, key=lambda c: -sizes.get(c, 0))
    for c in draw_order:
        s = coh[coh.cohort == c]
        ax.scatter(s.PC1, s.PC2, s=ps, c=cmap[c], alpha=0.65,
                   linewidths=0.15, edgecolors="white", zorder=3)

    style_axes(ax, df)
    # Legend swatches stay fully opaque and in COHORT_ORDER, regardless of the
    # draw order / alpha used for the points themselves.
    handles = [Line2D([], [], marker="o", linestyle="", markersize=3,
                      markerfacecolor=cmap[c], markeredgecolor="white",
                      markeredgewidth=0.2, label=f"{c} ({sizes.get(c, 0)})")
               for c in order]
    ax.legend(handles=handles, frameon=False, fontsize=5.2,
              loc="center left", bbox_to_anchor=(1.01, 0.5),
              handletextpad=0.3, labelspacing=0.25, borderaxespad=0)
    fig.tight_layout(pad=0.4)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(Path(out_dir) / f"figure1_panelB_cohort.{ext}",
                    dpi=400, bbox_inches="tight")
    plt.close(fig)


def panel_d(df, out_dir):
    """Per-cohort ancestry composition, stacked to 100%."""
    coh = df[df.dataset == "cohort"]
    order = [c for c in COHORT_ORDER if c in set(coh.cohort)]

    comp = (coh.groupby(["cohort", "genetic_ancestry"]).size()
               .unstack(fill_value=0).reindex(order))
    n_tot = comp.sum(axis=1)
    pct = comp.div(n_tot, axis=0) * 100
    cols = [a for a in ANC_ORDER if a in pct.columns]
    pct = pct[cols]

    fig, ax = plt.subplots(figsize=(3.1, 2.9))
    ypos = range(len(order))
    left = [0.0] * len(order)
    for anc in cols:
        vals = pct[anc].values
        ax.barh(list(ypos), vals, left=left, height=0.72,
                color=ANC_COLORS[anc], label=anc,
                edgecolor="white", linewidth=0.4)
        left = [l + v for l, v in zip(left, vals)]

    ax.set_yticks(list(ypos))
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Ancestry (%)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.tick_params(axis="x", length=2.5, width=0.6, pad=1.5)

    # Sample size at the end of each bar.
    for i, c in enumerate(order):
        ax.text(101, i, f"n={int(n_tot[c])}", va="center", ha="left", fontsize=5.5)

    ax.legend(frameon=False, fontsize=5.5, ncol=len(cols),
              loc="lower center", bbox_to_anchor=(0.5, 1.0),
              handletextpad=0.3, columnspacing=0.9, handlelength=1.0)
    fig.tight_layout(pad=0.4)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(Path(out_dir) / f"figure1_panelD_composition.{ext}",
                    dpi=400, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.samples, sep="\t")
    n_ref = (df.dataset == "reference").sum()
    n_coh = (df.dataset == "cohort").sum()
    print(f"[panels] {n_coh:,} cohort samples, {n_ref:,} reference samples")

    panel_a(df, out, args.point_size, args.ref_point_size)
    panel_b(df, out, args.point_size)
    panel_d(df, out)

    # Composition table alongside the figure, for the manuscript text.
    coh = df[df.dataset == "cohort"]
    comp = (coh.groupby(["cohort", "genetic_ancestry"]).size()
               .unstack(fill_value=0).reindex(
                   [c for c in COHORT_ORDER if c in set(coh.cohort)]))
    comp["TOTAL"] = comp.sum(axis=1)
    comp.to_csv(out / "figure1_panelD_counts.tsv", sep="\t")

    print(f"[panels] wrote panels A/B/D (svg, pdf, png) to {out}")


if __name__ == "__main__":
    main()
