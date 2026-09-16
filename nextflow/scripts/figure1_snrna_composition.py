#!/usr/bin/env python3
"""
Figure 1 cell-type panel: snRNA-seq cell type composition across datasets.

Unlike the bulk MGP scores - which are mean-centred within each cohort and so
cannot show absolute abundance - the snRNA-seq proportions are genuinely
compositional (they sum to 1 per donor), so cohort differences in cell type
composition are real and directly readable.

Input is manuscript_figure/combined_bulk_snrna_paired.tsv.

Two data caveats, both handled explicitly rather than silently:

  1. Mathys and ROSMAP share 213 donors (Mathys is a snRNA study of ROSMAP
     brains).  --drop-overlap removes the shared donors from Mathys so each
     donor is counted once; the default keeps them and labels the panel, since
     they are separate snRNA datasets of the same brains.

  2. Mathys lacks l5_et / pericyte / vlmc and Ruzicka lacks vlmc, so those
     donors sum to <1 (Mathys median 0.93).  Proportions are renormalised per
     donor over the cell types actually present, so each bar sums to 100%;
     the missing types are reported in the accompanying table.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.size"] = 7
matplotlib.rcParams["axes.linewidth"] = 0.6

# Colours follow manuscript_figure/figure2_ctp_variation_merged.R so the cell
# type panels agree across figures, and are deliberately NOT the blue/pink/green
# ancestry palette used by panels A/B/D - the same colour must not mean
# "European" in one panel and "excitatory neuron" in another.
#
# Figure 2's scheme: Okabe-Ito class colours (colourblind-safe), with a
# ColorBrewer ramp per class for the individual types -
#   Excitatory   Oranges,  Inhibitory   Greens,  Non-neuronal   Blues.
CLASS_COLORS = {
    "Excitatory":   "#D55E00",
    "Inhibitory":   "#009E73",
    "Non-neuronal": "#7570B3",
}

# Display order, grouped Excitatory -> Inhibitory -> Non-neuronal.
# Ramps are ColorBrewer Oranges/Greens/Blues levels 4-9, matching figure 2.
CELL_ORDER = [
    # Excitatory - Oranges
    ("it",             "IT",             "#8C2D04"),
    ("l4_it",          "L4 IT",          "#A63603"),
    ("l5_6_it_car3",   "L5/6 IT Car3",   "#D94801"),
    ("l5_et",          "L5 ET",          "#F16913"),
    ("l5_6_np",        "L5/6 NP",        "#FD8D3C"),
    ("l6_ct",          "L6 CT",          "#FDAE6B"),
    ("l6b",            "L6b",            "#FDD0A2"),
    # Inhibitory - Greens
    ("lamp5",          "LAMP5",          "#00441B"),
    ("pax6",           "PAX6",           "#006D2C"),
    ("vip",            "VIP",            "#238B45"),
    ("sst",            "SST",            "#41AB5D"),
    ("pvalb",          "PVALB",          "#74C476"),
    # Non-neuronal - Blues
    ("astrocyte",      "Astrocyte",      "#08306B"),
    ("oligodendrocyte","Oligodendrocyte","#08519C"),
    ("opc",            "OPC",            "#2171B5"),
    ("microglia",      "Microglia",      "#4292C6"),
    ("endothelial",    "Endothelial",    "#6BAED6"),
    ("pericyte",       "Pericyte",       "#9ECAE1"),
    ("vlmc",           "VLMC",           "#C6DBEF"),
]

COHORT_ORDER = ["ROSMAP", "Mathys", "MSBB", "HBCC", "Ruzicka"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paired-tsv",
                   default="manuscript_figure/combined_bulk_snrna_paired.tsv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--drop-overlap", action="store_true",
                   help="remove the 213 donors Mathys shares with ROSMAP")
    p.add_argument("--compact", action="store_true",
                   help="smaller figure for a tight Figure 1 slot")
    p.add_argument("--group", action="store_true",
                   help="collapse the 19 types to Excitatory / Inhibitory / "
                        "Non-neuronal (far more legible at small size)")
    p.add_argument("--vertical", action="store_true",
                   help="vertical bars (datasets on x, composition on y), "
                        "giving a portrait panel taller than it is wide")
    p.add_argument("--figsize", type=str, default=None,
                   help="override figure size as WxH in inches, e.g. 2.4x4.0")
    p.add_argument("--fontsize", type=float, default=8.5,
                   help="base text size in pt for axis and category labels; "
                        "n= labels and legends are 1.5 pt smaller "
                        "(default 8.5)")
    return p.parse_args()


# Coarse grouping used by --group.
COARSE = [
    ("Excitatory", CLASS_COLORS["Excitatory"],
     ["it", "l4_it", "l5_6_it_car3", "l5_et", "l5_6_np", "l6_ct", "l6b"]),
    ("Inhibitory", CLASS_COLORS["Inhibitory"],
     ["lamp5", "pax6", "vip", "sst", "pvalb"]),
    ("Non-neuronal", CLASS_COLORS["Non-neuronal"],
     ["astrocyte", "oligodendrocyte", "opc", "microglia",
      "endothelial", "pericyte", "vlmc"]),
]


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(args.paired_tsv, sep="\t")
    d = d.dropna(subset=["snrna_proportion"])

    note = ""
    if args.drop_overlap:
        shared = (d.groupby("sample_id").validation_cohort.nunique() > 1)
        shared = set(shared[shared].index)
        before = d.sample_id.nunique()
        d = d[~((d.validation_cohort == "Mathys") & (d.sample_id.isin(shared)))]
        note = (f"dropped {before - d.sample_id.nunique()} Mathys donors "
                f"shared with ROSMAP")
        print(f"[snrna] {note}")

    # Renormalise per donor over the cell types actually present, so every bar
    # sums to 100% even where a dataset lacks some rare types.
    tot = d.groupby(["validation_cohort", "sample_id"]).snrna_proportion.transform("sum")
    d = d.assign(prop_norm=d.snrna_proportion / tot)

    cohorts = [c for c in COHORT_ORDER if c in set(d.validation_cohort)]
    mean = (d.pivot_table(index="validation_cohort", columns="cell_type",
                          values="prop_norm", aggfunc="mean")
              .reindex(cohorts))

    if args.group:
        # Sum the fine types into the three coarse classes.
        agg = {}
        for label, colour, keys in COARSE:
            present = [k for k in keys if k in mean.columns]
            agg[label] = mean[present].fillna(0.0).sum(axis=1)
        mean = pd.DataFrame(agg)
        cells = [(label, label, colour) for label, colour, _ in COARSE]
    else:
        cells = [(k, lab, col) for k, lab, col in CELL_ORDER if k in mean.columns]

    n_don = d.groupby("validation_cohort").sample_id.nunique()

    # Portrait defaults for --vertical (taller than wide); landscape otherwise.
    if args.vertical:
        # Portrait, but only moderately so (h/w ~= 1.2), and at the same width
        # as the genotype panels, which are 3.0-3.1 in wide with aspect
        # 0.87-0.94.  A much narrower/taller panel (h/w ~ 1.9) does not sit
        # next to them comfortably.
        # Slimmer than the genotype panels but the same height, so panel E can
        # sit in a narrow column beside them.  Height is unchanged from the
        # previous 3.6 in; only the width came down (3.0 -> 2.3).
        figsize = (2.3, 3.6) if args.compact else (2.7, 4.1)
    else:
        figsize = (3.0, 2.2) if args.compact else (4.0, 2.8)
    if args.figsize:
        w, h = args.figsize.lower().split("x")
        figsize = (float(w), float(h))
    # Text sizes deliberately match figure1_genotype_panels.py so panels A/B/D
    # and E look like one figure when assembled: axis labels and tick/category
    # labels at the rcParams base (7pt), annotations and legends at 5.5pt.
    # --compact changes only the figure DIMENSIONS, never the type size.
    # The panel is narrower than A/B/D, so its text is set a little LARGER than
    # their 7 pt: once the panel is scaled to a common height in Inkscape, a
    # slimmer figure's type shrinks relative to the others, and 8.5/7 pt here
    # lands close to 7/5.5 pt there.  --fontsize overrides if the assembled
    # figure says otherwise.
    fs = args.fontsize
    fs_small = args.fontsize - 1.5               # n= labels, legends
    fig, ax = plt.subplots(figsize=figsize)

    base = pd.Series(0.0, index=mean.index)
    for key, label, colour in cells:
        vals = mean[key].fillna(0.0) * 100
        if args.vertical:
            ax.bar(range(len(mean)), vals, bottom=base.values, width=0.72,
                   color=colour, edgecolor="white", linewidth=0.3, label=label)
        else:
            ax.barh(range(len(mean)), vals, left=base.values, height=0.72,
                    color=colour, edgecolor="white", linewidth=0.3, label=label)
        base += vals

    axis_label = "snRNA-seq cell type composition (%)"
    if args.vertical:
        ax.set_xticks(range(len(mean)))
        ax.set_xticklabels(list(mean.index), fontsize=fs,
                           rotation=45, ha="right", rotation_mode="anchor")
        # Headroom above 100 for the rotated n= labels, which sit outside the
        # bars.  They are drawn vertically, so the room they need grows with
        # BOTH the text size and the digit count; a fixed 112 clipped them once
        # the type got larger.
        head = 12 + 1.6 * fs_small + 2.2 * max(len(str(v)) for v in n_don)
        ax.set_ylim(0, 100 + head)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_ylabel(axis_label, fontsize=fs)
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.tick_params(axis="x", length=0, pad=1.5, labelsize=fs)
        ax.tick_params(axis="y", length=2.5, width=0.6, pad=1.5, labelsize=fs)
        # n= above each bar, clear of the 100% top.
        for i, c in enumerate(mean.index):
            ax.text(i, 102.5, f"n={n_don[c]}", ha="center", va="bottom",
                    fontsize=fs_small, color="#333333", rotation=90,
                    clip_on=False)
    else:
        ax.set_yticks(range(len(mean)))
        ax.set_yticklabels(list(mean.index), fontsize=fs)
        ax.invert_yaxis()
        # Headroom to the right of 100 for the n= labels (horizontal here, so
        # the room needed scales mainly with the digit count).
        ax.set_xlim(0, 100 + 4 + 2.6 * fs_small
                    + 2.0 * max(len(str(v)) for v in n_don))
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel(axis_label, fontsize=fs)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=1.5)
        ax.tick_params(axis="x", length=2.5, width=0.6, pad=1.5, labelsize=fs)
        # n= to the right of each bar, outside it.
        for i, c in enumerate(mean.index):
            ax.text(102, i, f"n={n_don[c]}", va="center", ha="left",
                    fontsize=fs_small, color="#333333")

    if args.group:
        # Three entries above the panel.  "Excitatory / Inhibitory /
        # Non-neuronal" on one row needs roughly 0.42 in per pt of text; below
        # that the last label is clipped, since the canvas is not grown by
        # bbox_inches.  So wrap to 2 columns on a narrow panel.
        one_row_in = 0.42 * fs_small
        ncol = 3 if figsize[0] >= one_row_in else 2
        ax.legend(frameon=False, fontsize=fs_small, ncol=ncol,
                  loc="lower center", bbox_to_anchor=(0.5, 1.0),
                  handletextpad=0.4, columnspacing=1.0, handlelength=1.0)
    elif args.vertical:
        # 19 entries beside a narrow portrait panel would dwarf it, so the
        # legend goes underneath in two columns.
        ax.legend(frameon=False, fontsize=fs_small, ncol=2,
                  loc="upper center", bbox_to_anchor=(0.5, -0.16),
                  handletextpad=0.4, labelspacing=0.22, columnspacing=0.8,
                  handlelength=0.9)
    else:
        ax.legend(frameon=False, fontsize=fs_small, ncol=1,
                  loc="center left", bbox_to_anchor=(1.02, 0.5),
                  handletextpad=0.4, labelspacing=0.22, handlelength=0.9)

    fig.tight_layout(pad=0.3)
    stem = ("figure1_panelE_snrna_composition_grouped" if args.group
            else "figure1_panelE_snrna_composition")
    if args.vertical:
        stem += "_vertical"

    # NOTE: no bbox_inches="tight" here.  That option grows the canvas to fit
    # rotated tick labels and the legend, which silently turns a portrait
    # figsize into a near-square file.  tight_layout() already fits everything
    # INSIDE the requested figure, so the saved file keeps the exact WxH asked
    # for and the portrait aspect ratio is preserved.
    for ext in ("svg", "pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=400)
    plt.close(fig)

    w, h = fig.get_size_inches()
    print(f"[snrna] figure size: {w:.2f} x {h:.2f} in "
          f"(aspect h/w = {h / w:.2f})")

    # Table behind the panel, plus what each dataset is missing.
    t = (mean * 100).round(2)
    t.insert(0, "n_donors", n_don.reindex(t.index))
    t.to_csv(out_dir / f"{stem}.tsv", sep="\t")

    raw = pd.read_csv(args.paired_tsv, sep="\t")
    miss = (raw.pivot_table(index="validation_cohort", columns="cell_type",
                            values="snrna_proportion", aggfunc="mean")
               .isna())
    with open(out_dir / "figure1_panelE_snrna_notes.txt", "w") as fh:
        fh.write("snRNA composition panel notes\n")
        fh.write("=============================\n")
        fh.write("Proportions renormalised per donor over cell types present,\n")
        fh.write("so each bar sums to 100%.\n\n")
        for c in miss.index:
            m = list(miss.columns[miss.loc[c]])
            fh.write(f"{c}: missing {m if m else 'none'}\n")
        fh.write(f"\nMathys and ROSMAP share 213 donors "
                 f"(Mathys is a snRNA study of ROSMAP brains).\n")
        fh.write(f"overlap handling: {note if note else 'kept (not dropped)'}\n")

    print(f"[snrna] wrote {stem}.(svg|pdf|png) to {out_dir}")
    print(f"[snrna] donors per dataset: {n_don.to_dict()}")


if __name__ == "__main__":
    main()
