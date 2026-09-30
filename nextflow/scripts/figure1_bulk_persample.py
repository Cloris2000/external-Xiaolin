#!/usr/bin/env python3
"""
Per-sample cell type profile across the 15 bulk cohorts (the main datasets).

One column per donor, one row per cell type, grouped by cohort - the same
layout as a per-sample stacked bar, but colour encodes the relative MGP score
instead of a stacked share.

Why not a stacked 100% bar for the bulk cohorts
------------------------------------------------
The bulk phenotypes come from MGP (marker gene profiles), which returns the
first principal component of each cell type's marker genes.  Those scores are
*relative*: they are mean-centred within cohort, take negative values, and do
not sum to a fixed total per donor (observed row sums range from about -96 to
+15).  A stacked 100% bar needs non-negative parts that sum to a meaningful
whole, so it cannot be built from this data without inventing an absolute scale
the deconvolution never estimated.

What this panel shows instead is the quantity the GWAS actually uses: how far
each donor sits from its cohort mean, per cell type, in SD units.  That is
comparable across all 15 cohorts, and it still answers "how does cellular
makeup vary within and between cohorts" - it just answers it in relative rather
than absolute terms.

Colours follow manuscript_figure/figure2_ctp_variation_merged.R for the cell
type class labels; the heatmap itself uses a diverging blue-white-red scale
centred on each cohort's mean.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["axes.linewidth"] = 0.6

CLASS_COLORS = {
    "Excitatory":   "#D55E00",
    "Inhibitory":   "#009E73",
    "Non-neuronal": "#7570B3",
}

# Rows, grouped Excitatory -> Inhibitory -> Non-neuronal.
ROWS = [
    ("IT",              "IT",             "Excitatory"),
    ("L4.IT",           "L4 IT",          "Excitatory"),
    ("L5.6.IT.Car3",    "L5/6 IT Car3",   "Excitatory"),
    ("L5.ET",           "L5 ET",          "Excitatory"),
    ("L5.6.NP",         "L5/6 NP",        "Excitatory"),
    ("L6.CT",           "L6 CT",          "Excitatory"),
    ("L6b",             "L6b",            "Excitatory"),
    ("LAMP5",           "LAMP5",          "Inhibitory"),
    ("PAX6",            "PAX6",           "Inhibitory"),
    ("VIP",             "VIP",            "Inhibitory"),
    ("SST",             "SST",            "Inhibitory"),
    ("PVALB",           "PVALB",          "Inhibitory"),
    ("Astrocyte",       "Astrocyte",      "Non-neuronal"),
    ("Oligodendrocyte", "Oligodendrocyte","Non-neuronal"),
    ("OPC",             "OPC",            "Non-neuronal"),
    ("Microglia",       "Microglia",      "Non-neuronal"),
    ("Endothelial",     "Endothelial",    "Non-neuronal"),
    ("Pericyte",        "Pericyte",       "Non-neuronal"),
    ("VLMC",            "VLMC",           "Non-neuronal"),
]

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

# Short labels keep the x axis readable at 15 cohorts.
SHORT = {
    "ROSMAP": "ROSMAP", "ROSMAP_array": "ROSMAP array",
    "Mayo": "Mayo", "MSBB": "MSBB",
    "CMC_MSSM": "CMC MSSM", "CMC_PENN": "CMC PENN",
    "CMC_PITT": "CMC PITT", "GTEx_v10": "GTEx v10",
    "NABEC": "NABEC", "NIMH_HBCC_1M": "HBCC 1M",
    "NIMH_HBCC_h650": "HBCC h650", "NIMH_HBCC_Omni5M": "HBCC Omni5M",
    "GVEX": "GVEX", "AMP_AD_Rush": "AMP-AD Rush",
    "AMP_AD_Mayo": "AMP-AD Mayo",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geno-root",
                   default="/scratch/zhoux156/results/cohort_gwas_v2/results")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sort-by", default="IT",
                   help="cell type used to order donors within each cohort, "
                        "or 'none' (default IT, the largest excitatory class)")
    p.add_argument("--vmax", type=float, default=2.5,
                   help="colour scale limit in SD (default 2.5)")
    p.add_argument("--figsize", type=str, default="7.0x3.4")
    p.add_argument("--fontsize", type=float, default=7.0)
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # The *_scaled file is z-scored within cohort: mean 0, SD 1 per cell type.
    # That is exactly the "deviation from this cohort's mean" this panel shows.
    mats, counts = {}, {}
    for cohort, _ in COHORTS:
        f = Path(args.geno_root) / cohort / "cell_proportions_scaled.csv"
        if not f.exists():
            print(f"[bulk] WARNING: missing {f}")
            continue
        df = pd.read_csv(f, index_col=0)
        keep = [k for k, _, _ in ROWS if k in df.columns]
        mats[cohort] = df[keep]
        counts[cohort] = len(df)
    if not mats:
        raise SystemExit("ERROR: no cell_proportions_scaled.csv found")

    cohorts = [c for c, _ in COHORTS if c in mats]
    rows = [(k, lab, cls) for k, lab, cls in ROWS if k in mats[cohorts[0]].columns]

    fw, fh = (float(x) for x in args.figsize.lower().split("x"))
    fs = args.fontsize
    fig, ax = plt.subplots(figsize=(fw, fh))

    gap = max(4, int(0.004 * sum(counts.values())))
    norm = TwoSlopeNorm(vmin=-args.vmax, vcenter=0, vmax=args.vmax)

    x = 0
    centers, ticks, bounds = [], [], []
    for c in cohorts:
        m = mats[c]
        if args.sort_by != "none" and args.sort_by in m.columns:
            m = m.loc[m[args.sort_by].sort_values(ascending=False).index]
        arr = m[[k for k, _, _ in rows]].values.T      # cell types x donors
        n = arr.shape[1]
        ax.imshow(arr, aspect="auto", cmap="RdBu_r", norm=norm,
                  extent=(x, x + n, len(rows), 0), interpolation="nearest")
        centers.append(x + n / 2)
        ticks.append(f"{SHORT.get(c, c)}\nn={counts[c]}")
        x += n
        bounds.append(x)
        x += gap

    # Thin separators between cohorts.
    for b in bounds[:-1]:
        ax.axvline(b + gap / 2, color="white", linewidth=1.2)

    ax.set_xlim(0, x - gap)
    ax.set_ylim(len(rows), 0)
    ax.set_yticks(np.arange(len(rows)) + 0.5)
    ax.set_yticklabels([lab for _, lab, _ in rows], fontsize=fs - 0.5)
    for t, (_, _, cls) in zip(ax.get_yticklabels(), rows):
        t.set_color(CLASS_COLORS[cls])
    # Cohort widths are proportional to n, so the small cohorts' labels collide
    # if they all sit on one line.  Rotate them and let each label start at its
    # block, rather than centring a two-line label in a 95-donor-wide column.
    ax.set_xticks(centers)
    ax.set_xticklabels([t.replace("\n", " ") for t in ticks],
                       fontsize=fs - 1.5, rotation=45,
                       ha="right", rotation_mode="anchor")
    ax.tick_params(axis="both", length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)

    # Class brackets down the right edge.
    y = 0
    for cls in ["Excitatory", "Inhibitory", "Non-neuronal"]:
        k = sum(1 for _, _, c in rows if c == cls)
        if k:
            ax.plot([x - gap + 2, x - gap + 2], [y + 0.15, y + k - 0.15],
                    color=CLASS_COLORS[cls], linewidth=1.6,
                    clip_on=False, solid_capstyle="butt")
            ax.text(x - gap + 5, y + k / 2, cls, rotation=270, va="center",
                    ha="left", fontsize=fs - 1, color=CLASS_COLORS[cls])
            y += k

    cb = fig.colorbar(ax.images[0], ax=ax, pad=0.10, fraction=0.018,
                      ticks=[-args.vmax, 0, args.vmax])
    cb.set_label("MGP score (SD from cohort mean)", fontsize=fs - 1)
    cb.ax.tick_params(labelsize=fs - 1.5, length=2, width=0.6)
    cb.outline.set_visible(False)

    sub = (f"each column = one donor, sorted by {args.sort_by} within cohort"
           if args.sort_by != "none" else "each column = one donor")
    ax.set_xlabel(f"{sum(counts.values()):,} donors across "
                  f"{len(cohorts)} cohorts - {sub}", fontsize=fs - 0.5)

    fig.tight_layout(pad=0.3)
    stem = "figure1_bulk_persample_celltypes"
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=400)
    plt.close(fig)

    # Per-cohort spread on the RAW (un-z-scored) MGP scores, where the spread
    # is still informative - the scaled file is 1.0 by construction.
    raw = {}
    for c in cohorts:
        f = Path(args.geno_root) / c / "cell_proportions.csv"
        if f.exists():
            raw[c] = pd.read_csv(f, index_col=0).std()
    if raw:
        sd = pd.DataFrame(raw).T
        sd.index.name = "cohort"
        sd.round(3).to_csv(out_dir / f"{stem}_sd.tsv", sep="\t")

    pd.Series(counts, name="n_samples").rename_axis("cohort").to_csv(
        out_dir / f"{stem}_counts.tsv", sep="\t")

    print(f"[bulk] {sum(counts.values()):,} donors across {len(cohorts)} cohorts")
    print(f"[bulk] figure {fw} x {fh} in -> {stem}.(svg|pdf|png)")


if __name__ == "__main__":
    main()
