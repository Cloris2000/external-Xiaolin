#!/usr/bin/env python3
"""
Cell type variation across the 15 bulk cohorts - a compact, readable summary.

Replaces the per-donor heatmap, which packed 6,399 columns into a few inches
and read as noise.  Two panels, both one mark per cohort (or per cohort x cell
type), so everything is legible at figure size:

  (a) Neuron-glia axis: box + strip of each donor's neuronal minus glial mean
      MGP score, one box per cohort.  Shows how much the cellular makeup varies
      WITHIN each cohort, and lets cohorts be compared directly.

  (b) Per-cohort spread (SD of the raw MGP score) per cell type, as a cohort x
      cell type heatmap - 15 x 19 cells rather than 6,399 columns.

Why the axis is relative, not absolute
--------------------------------------
The bulk phenotypes come from MGP, which returns the first principal component
of each cell type's marker genes: relative scores, mean-centred within cohort,
negative values allowed, no fixed per-donor total.  So absolute composition
("this cohort is 40% astrocyte") is not estimable from this data, and every
cohort mean is 0 by construction.  What differs - and what these panels show -
is the SPREAD, which is exactly the phenotypic variance the GWAS is powered on.

Colours follow manuscript_figure/figure2_ctp_variation_merged.R.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["axes.linewidth"] = 0.6

CLASS_COLORS = {
    "Excitatory":   "#D55E00",
    "Inhibitory":   "#009E73",
    "Non-neuronal": "#7570B3",
}

NEURONAL = ["IT", "L4.IT", "L5.6.IT.Car3", "L5.ET", "L5.6.NP", "L6.CT", "L6b",
            "LAMP5", "PAX6", "VIP", "SST", "PVALB"]
GLIAL = ["Astrocyte", "Oligodendrocyte", "OPC", "Microglia",
         "Endothelial", "Pericyte", "VLMC"]

ROWS = [
    ("IT", "IT", "Excitatory"), ("L4.IT", "L4 IT", "Excitatory"),
    ("L5.6.IT.Car3", "L5/6 IT Car3", "Excitatory"),
    ("L5.ET", "L5 ET", "Excitatory"), ("L5.6.NP", "L5/6 NP", "Excitatory"),
    ("L6.CT", "L6 CT", "Excitatory"), ("L6b", "L6b", "Excitatory"),
    ("LAMP5", "LAMP5", "Inhibitory"), ("PAX6", "PAX6", "Inhibitory"),
    ("VIP", "VIP", "Inhibitory"), ("SST", "SST", "Inhibitory"),
    ("PVALB", "PVALB", "Inhibitory"),
    ("Astrocyte", "Astrocyte", "Non-neuronal"),
    ("Oligodendrocyte", "Oligodendrocyte", "Non-neuronal"),
    ("OPC", "OPC", "Non-neuronal"), ("Microglia", "Microglia", "Non-neuronal"),
    ("Endothelial", "Endothelial", "Non-neuronal"),
    ("Pericyte", "Pericyte", "Non-neuronal"), ("VLMC", "VLMC", "Non-neuronal"),
]

COHORTS = ["ROSMAP", "ROSMAP_array", "Mayo", "MSBB", "CMC_MSSM", "CMC_PENN",
           "CMC_PITT", "GTEx_v10", "NABEC", "NIMH_HBCC_1M", "NIMH_HBCC_h650",
           "NIMH_HBCC_Omni5M", "GVEX", "AMP_AD_Rush", "AMP_AD_Mayo"]

SHORT = {
    "ROSMAP": "ROSMAP", "ROSMAP_array": "ROSMAP array", "Mayo": "Mayo",
    "MSBB": "MSBB", "CMC_MSSM": "CMC MSSM", "CMC_PENN": "CMC PENN",
    "CMC_PITT": "CMC PITT", "GTEx_v10": "GTEx v10", "NABEC": "NABEC",
    "NIMH_HBCC_1M": "HBCC 1M", "NIMH_HBCC_h650": "HBCC h650",
    "NIMH_HBCC_Omni5M": "HBCC Omni5M", "GVEX": "GVEX",
    "AMP_AD_Rush": "AMP-AD Rush", "AMP_AD_Mayo": "AMP-AD Mayo",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geno-root",
                   default="/scratch/zhoux156/results/cohort_gwas_v2/results")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--panel", choices=["axis", "spread", "both"], default="both")
    p.add_argument("--figsize", type=str, default=None)
    p.add_argument("--fontsize", type=float, default=7.0)
    return p.parse_args()


def load(root):
    scaled, raw, n = {}, {}, {}
    for c in COHORTS:
        fs = Path(root) / c / "cell_proportions_scaled.csv"
        fr = Path(root) / c / "cell_proportions.csv"
        if not fs.exists():
            print(f"[summary] WARNING: missing {fs}")
            continue
        scaled[c] = pd.read_csv(fs, index_col=0)
        if fr.exists():
            raw[c] = pd.read_csv(fr, index_col=0)
        n[c] = len(scaled[c])
    if not scaled:
        raise SystemExit("ERROR: no cell_proportions_scaled.csv found")
    return scaled, raw, n


def panel_axis(ax, scaled, n, cohorts, fs):
    """Neuron-glia axis per donor, one box + strip per cohort."""
    rng = np.random.default_rng(0)
    data = []
    for c in cohorts:
        d = scaled[c]
        neu = [k for k in NEURONAL if k in d.columns]
        gli = [k for k in GLIAL if k in d.columns]
        data.append((d[neu].mean(axis=1) - d[gli].mean(axis=1)).values)

    for i, v in enumerate(data):
        # Strip of donors behind the box, subsampled so dense cohorts do not
        # turn into a solid block.
        show = v if len(v) <= 400 else rng.choice(v, 400, replace=False)
        ax.scatter(rng.normal(i, 0.055, len(show)), show, s=1.1,
                   color="#9AA3AF", alpha=0.45, linewidths=0, zorder=1)

    bp = ax.boxplot(data, positions=range(len(cohorts)), widths=0.6,
                    showfliers=False, patch_artist=True, zorder=3)
    for box in bp["boxes"]:
        box.set(facecolor="white", edgecolor="#333333", linewidth=0.7, alpha=0.9)
    for part in ("whiskers", "caps"):
        for ln in bp[part]:
            ln.set(color="#333333", linewidth=0.7)
    for md in bp["medians"]:
        md.set(color="#D55E00", linewidth=1.3)

    ax.axhline(0, color="#999999", linewidth=0.5, linestyle=(0, (3, 3)), zorder=0)
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels([f"{SHORT.get(c, c)}\nn={n[c]}" for c in cohorts],
                       fontsize=fs - 1.5, rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel("Neuronal - glial MGP score (SD)", fontsize=fs)
    ax.tick_params(axis="y", length=2.5, width=0.6, pad=1.5, labelsize=fs - 0.5)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.spines[["top", "right"]].set_visible(False)


def panel_spread(ax, raw, cohorts, fs, fig):
    """Cohort x cell type heatmap of the raw-score SD."""
    keys = [k for k, _, _ in ROWS if k in raw[cohorts[0]].columns]
    labs = [lab for k, lab, _ in ROWS if k in raw[cohorts[0]].columns]
    clss = [cls for k, _, cls in ROWS if k in raw[cohorts[0]].columns]

    M = np.array([[raw[c][k].std() for k in keys] for c in cohorts])

    im = ax.imshow(M, aspect="auto", cmap="viridis", interpolation="nearest")
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels(labs, fontsize=fs - 1.5, rotation=45, ha="right",
                       rotation_mode="anchor")
    for t, cls in zip(ax.get_xticklabels(), clss):
        t.set_color(CLASS_COLORS[cls])
    ax.set_yticks(range(len(cohorts)))
    ax.set_yticklabels([SHORT.get(c, c) for c in cohorts], fontsize=fs - 1.5)
    ax.tick_params(length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)

    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.030)
    cb.set_label("SD of MGP score", fontsize=fs - 1)
    cb.ax.tick_params(labelsize=fs - 1.5, length=2, width=0.6)
    cb.outline.set_visible(False)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scaled, raw, n = load(args.geno_root)
    cohorts = [c for c in COHORTS if c in scaled]
    fs = args.fontsize

    if args.panel in ("axis", "both"):
        fw, fh = ((float(x) for x in args.figsize.lower().split("x"))
                  if args.figsize and args.panel == "axis" else (5.2, 2.6))
        fig, ax = plt.subplots(figsize=(fw, fh))
        panel_axis(ax, scaled, n, cohorts, fs)
        fig.tight_layout(pad=0.3)
        for ext in ("svg", "pdf", "png"):
            fig.savefig(out_dir / f"figure1_bulk_neuron_glia_axis.{ext}", dpi=400)
        plt.close(fig)
        print("[summary] wrote figure1_bulk_neuron_glia_axis.(svg|pdf|png)")

    if args.panel in ("spread", "both") and raw:
        fw, fh = ((float(x) for x in args.figsize.lower().split("x"))
                  if args.figsize and args.panel == "spread" else (5.6, 3.0))
        fig, ax = plt.subplots(figsize=(fw, fh))
        panel_spread(ax, raw, cohorts, fs, fig)
        fig.tight_layout(pad=0.3)
        for ext in ("svg", "pdf", "png"):
            fig.savefig(out_dir / f"figure1_bulk_celltype_spread.{ext}", dpi=400)
        plt.close(fig)
        print("[summary] wrote figure1_bulk_celltype_spread.(svg|pdf|png)")

    # Tables behind both panels.
    rows = []
    for c in cohorts:
        d = scaled[c]
        neu = [k for k in NEURONAL if k in d.columns]
        gli = [k for k in GLIAL if k in d.columns]
        v = d[neu].mean(axis=1) - d[gli].mean(axis=1)
        rows.append(dict(cohort=c, n=n[c], sd=round(v.std(), 3),
                         iqr=round(v.quantile(.75) - v.quantile(.25), 3),
                         q25=round(v.quantile(.25), 3),
                         median=round(v.median(), 3),
                         q75=round(v.quantile(.75), 3)))
    pd.DataFrame(rows).to_csv(
        out_dir / "figure1_bulk_neuron_glia_axis.tsv", sep="\t", index=False)

    if raw:
        sd = pd.DataFrame({c: raw[c].std() for c in cohorts}).T
        sd.index.name = "cohort"
        sd.round(3).to_csv(out_dir / "figure1_bulk_celltype_spread.tsv", sep="\t")

    print(f"[summary] {sum(n.values()):,} donors across {len(cohorts)} cohorts")


if __name__ == "__main__":
    main()
