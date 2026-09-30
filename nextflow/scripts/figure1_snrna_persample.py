#!/usr/bin/env python3
"""
Per-sample stacked cell type composition, grouped by dataset.

One thin vertical bar per donor, stacked to 100%, with datasets side by side.
This shows what a cohort-mean bar cannot: how variable the cellular makeup is
WITHIN each dataset, and whether a dataset's mean is representative or is being
pulled by a subset of donors.

Samples are sorted within each dataset (by default on excitatory-neuron
fraction), so each dataset reads as a gradient.  That ordering is cosmetic - it
carries no donor metadata - but it makes the spread and any sub-structure
immediately visible, and keeps datasets comparable.

Data: manuscript_figure/combined_bulk_snrna_paired.tsv (snRNA-seq proportions,
genuinely compositional).  The 15 bulk cohorts CANNOT be shown this way: their
MGP scores are mean-centred within cohort, take negative values and do not sum
to 1 per donor, so there is no per-sample composition to stack.

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

CELL_ORDER = [
    ("it",             "IT",             "#8C2D04"),
    ("l4_it",          "L4 IT",          "#A63603"),
    ("l5_6_it_car3",   "L5/6 IT Car3",   "#D94801"),
    ("l5_et",          "L5 ET",          "#F16913"),
    ("l5_6_np",        "L5/6 NP",        "#FD8D3C"),
    ("l6_ct",          "L6 CT",          "#FDAE6B"),
    ("l6b",            "L6b",            "#FDD0A2"),
    ("lamp5",          "LAMP5",          "#00441B"),
    ("pax6",           "PAX6",           "#006D2C"),
    ("vip",            "VIP",            "#238B45"),
    ("sst",            "SST",            "#41AB5D"),
    ("pvalb",          "PVALB",          "#74C476"),
    ("astrocyte",      "Astrocyte",      "#08306B"),
    ("oligodendrocyte","Oligodendrocyte","#08519C"),
    ("opc",            "OPC",            "#2171B5"),
    ("microglia",      "Microglia",      "#4292C6"),
    ("endothelial",    "Endothelial",    "#6BAED6"),
    ("pericyte",       "Pericyte",       "#9ECAE1"),
    ("vlmc",           "VLMC",           "#C6DBEF"),
]

COARSE = [
    ("Excitatory", CLASS_COLORS["Excitatory"],
     ["it", "l4_it", "l5_6_it_car3", "l5_et", "l5_6_np", "l6_ct", "l6b"]),
    ("Inhibitory", CLASS_COLORS["Inhibitory"],
     ["lamp5", "pax6", "vip", "sst", "pvalb"]),
    ("Non-neuronal", CLASS_COLORS["Non-neuronal"],
     ["astrocyte", "oligodendrocyte", "opc", "microglia",
      "endothelial", "pericyte", "vlmc"]),
]

DATASET_ORDER = ["ROSMAP", "Mathys", "MSBB", "HBCC", "Ruzicka"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paired-tsv",
                   default="manuscript_figure/combined_bulk_snrna_paired.tsv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--group", action="store_true",
                   help="collapse to Excitatory / Inhibitory / Non-neuronal")
    p.add_argument("--drop-overlap", action="store_true",
                   help="remove the 213 donors Mathys shares with ROSMAP")
    p.add_argument("--sort-by", default="Excitatory",
                   choices=["Excitatory", "Inhibitory", "Non-neuronal", "none"],
                   help="sort donors within each dataset by this fraction")
    p.add_argument("--figsize", type=str, default="6.0x2.6",
                   help="WxH inches (default 6.0x2.6)")
    p.add_argument("--fontsize", type=float, default=7.0)
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(args.paired_tsv, sep="\t").dropna(subset=["snrna_proportion"])

    note = "kept"
    if args.drop_overlap:
        shared = d.groupby("sample_id").validation_cohort.nunique() > 1
        shared = set(shared[shared].index)
        d = d[~((d.validation_cohort == "Mathys") & (d.sample_id.isin(shared)))]
        note = f"dropped {len(shared)} Mathys donors shared with ROSMAP"
        print(f"[persample] {note}")

    # Wide matrix: one row per (dataset, donor), one column per cell type.
    w = d.pivot_table(index=["validation_cohort", "sample_id"],
                      columns="cell_type", values="snrna_proportion",
                      aggfunc="mean").fillna(0.0)

    # Renormalise each donor over the cell types present, so every bar is 100%.
    w = w.div(w.sum(axis=1), axis=0)

    if args.group:
        agg = {}
        for label, colour, keys in COARSE:
            present = [k for k in keys if k in w.columns]
            agg[label] = w[present].sum(axis=1)
        w = pd.DataFrame(agg)
        cells = [(lab, lab, col) for lab, col, _ in COARSE]
    else:
        cells = [(k, lab, col) for k, lab, col in CELL_ORDER if k in w.columns]

    # Fraction used for ordering donors inside each dataset.
    if args.sort_by != "none":
        keys = dict((lab, ks) for lab, _, ks in COARSE)[args.sort_by]
        if args.group:
            sort_val = w[args.sort_by]
        else:
            sort_val = w[[k for k in keys if k in w.columns]].sum(axis=1)
    else:
        sort_val = pd.Series(0.0, index=w.index)

    datasets = [c for c in DATASET_ORDER
                if c in w.index.get_level_values(0).unique()]

    fw, fh = (float(x) for x in args.figsize.lower().split("x"))
    fs = args.fontsize
    fig, ax = plt.subplots(figsize=(fw, fh))

    # Lay the datasets out along x with a gap between them.
    gap = max(6, int(0.02 * len(w)))
    x = 0
    centers, ticks = [], []
    for ds in datasets:
        sub = w.loc[ds]
        order = sort_val.loc[ds].sort_values(ascending=False).index
        sub = sub.loc[order]
        n = len(sub)
        xs = np.arange(x, x + n)

        bottom = np.zeros(n)
        for key, label, colour in cells:
            vals = (sub[key].values if key in sub.columns
                    else np.zeros(n)) * 100
            ax.bar(xs, vals, bottom=bottom, width=1.0,
                   color=colour, linewidth=0, label=label if ds == datasets[0] else None)
            bottom += vals

        centers.append(x + n / 2)
        ticks.append(f"{ds}\n(n={n})")
        x += n + gap

    ax.set_xlim(-gap / 2, x - gap / 2)
    ax.set_ylim(0, 100)
    ax.set_xticks(centers)
    ax.set_xticklabels(ticks, fontsize=fs)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Cell type composition (%)", fontsize=fs)
    ax.set_xlabel("Donors, sorted by "
                  f"{args.sort_by.lower()} fraction within dataset"
                  if args.sort_by != "none" else "Donors",
                  fontsize=fs - 0.5)
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", length=0, pad=2, labelsize=fs)
    ax.tick_params(axis="y", length=2.5, width=0.6, pad=1.5, labelsize=fs)

    ncol = 3 if args.group else 5
    ax.legend(frameon=False, fontsize=fs - 1.5, ncol=ncol,
              loc="lower center", bbox_to_anchor=(0.5, 1.0),
              handletextpad=0.4, columnspacing=1.0, handlelength=1.0)

    fig.tight_layout(pad=0.3)
    stem = ("figure1_snrna_persample_grouped" if args.group
            else "figure1_snrna_persample")
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=400)
    plt.close(fig)

    # Per-dataset spread, the quantity the panel is really showing.
    rows = []
    for ds in datasets:
        sub = w.loc[ds]
        for key, label, _ in cells:
            if key not in sub.columns:
                continue
            v = sub[key] * 100
            rows.append(dict(dataset=ds, cell_type=label, n=len(v),
                             mean=round(v.mean(), 2), sd=round(v.std(), 2),
                             q25=round(v.quantile(.25), 2),
                             median=round(v.median(), 2),
                             q75=round(v.quantile(.75), 2)))
    pd.DataFrame(rows).to_csv(out_dir / f"{stem}_stats.tsv", sep="\t", index=False)

    print(f"[persample] {len(w)} donors across {len(datasets)} datasets")
    print(f"[persample] figure {fw} x {fh} in -> {stem}.(svg|pdf|png)")
    print(f"[persample] Mathys/ROSMAP overlap: {note}")


if __name__ == "__main__":
    main()
