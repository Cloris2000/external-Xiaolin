#!/usr/bin/env python3
"""Publication Figure 7a: per-supertype ICC, -SEAAD types labelled.

cnsplots Nature styling. No title, subtitle, or caption; methods live in the
legend. Writes under compare_shreejoy_with_nextflow/figures/.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D

import cnsplots as cns

OUT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow")
FIG = OUT / "figures"
WP3 = Path("/scratch/shreejoy/cell_type_bias/wp3/results")

INDEP = {
    "Green x Mathys",
    "Green x Multiome2025",
    "Mathys x Multiome2025",
    "Multiome2025 x PsychAD",
    "BrainSCOPE CMC x PsychAD",
    "PsychAD x Ruzicka (MSSM2 x MSSM1)",
}

# Nature cycle: vermillion for glia (drawn on top), navy for neurons.
NEU = "#3C5488"
GLIA = "#E64B35"
NEU_MUTED = "#8491B4"
IQR = "#B7B7B7"


def load_table() -> pd.DataFrame:
    wl = pd.read_csv(WP3 / "93_whitelist_v2.csv")
    wl["seaad"] = wl["seaad"].astype(bool)
    wl["reactive"] = wl["seaad"] | wl["cell_type"].str.endswith(("SEAAD", "SEAD"))
    wl["single_pair"] = False
    wl["independent"] = True

    pc = pd.read_csv(WP3 / "51_pair_percelltype_supertype.csv")
    seaad = pc[pc["cell_type"].str.contains("SEAAD|SEAD", regex=True)]
    missing = sorted(set(seaad["cell_type"]) - set(wl["cell_type"]))

    extra_rows = []
    for ct in missing:
        rows = seaad[seaad["cell_type"] == ct]
        indep = rows[rows["pair"].isin(INDEP)]
        use = indep if len(indep) else rows
        icc = use["icc_person_a"].to_numpy()
        extra_rows.append({
            "cell_type": ct,
            "n_pairs": len(use),
            "median_icc": float(np.median(icc)),
            "min_icc": float(np.min(icc)),
            "q25_icc": float(np.quantile(icc, 0.25)) if len(use) >= 2
                       else float(use["icc_ci_lo"].iloc[0]),
            "q75_icc": float(np.quantile(icc, 0.75)) if len(use) >= 2
                       else float(use["icc_ci_hi"].iloc[0]),
            "group": use["group"].iloc[0],
            "seaad": True,
            "reactive": True,
            "single_pair": True,
            "independent": bool(use["pair"].isin(INDEP).all()),
        })
    if extra_rows:
        wl = pd.concat([wl, pd.DataFrame(extra_rows)], ignore_index=True)

    wl["class_lab"] = np.where(wl["group"] == "glia", "Non-neuronal", "Neuronal")
    wl = wl.sort_values("median_icc", ascending=False).reset_index(drop=True)
    wl["rank"] = np.arange(1, len(wl) + 1)
    return wl


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    wl = load_table()

    neu = wl[wl["class_lab"] == "Neuronal"]
    glia = wl[(wl["class_lab"] == "Non-neuronal") & ~wl["reactive"]]
    react_multi = wl[wl["reactive"] & ~wl["single_pair"]]
    react_one = wl[wl["reactive"] & wl["single_pair"]]
    react = wl[wl["reactive"]]

    cns.figure(width=520, height=280, color_cycle="Nature")
    ax = plt.gca()

    ax.vlines(wl["rank"], wl["q25_icc"], wl["q75_icc"],
              colors=IQR, linewidth=0.6, zorder=1)
    ax.scatter(neu["rank"], neu["median_icc"], s=12, c=NEU_MUTED,
               edgecolors=NEU, linewidths=0.2, zorder=2, label="_nolegend_")
    ax.scatter(glia["rank"], glia["median_icc"], s=22, c=GLIA,
               edgecolors="#333333", linewidths=0.35, zorder=3, label="_nolegend_")
    ax.scatter(react_multi["rank"], react_multi["median_icc"], s=36,
               c=GLIA, marker="D", edgecolors="#111111", linewidths=0.5,
               zorder=4, label="_nolegend_")
    ax.scatter(react_one["rank"], react_one["median_icc"], s=36,
               facecolors="white", marker="D", edgecolors=GLIA,
               linewidths=0.9, zorder=5, label="_nolegend_")

    texts = [
        ax.text(r.rank, r.median_icc, r.cell_type, fontsize=6,
                color="#7A1F14", fontweight="medium")
        for r in react.itertuples()
    ]
    adjust_text(
        texts, ax=ax,
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        expand=(1.15, 1.3),
    )

    ax.set_xlabel("Supertypes")
    ax.set_ylabel("Intraclass correlation")
    ax.set_title("")
    ax.set_xlim(0, wl["rank"].max() + 1)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=NEU_MUTED,
               markeredgecolor=NEU, markersize=5.5, label="Neuronal"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GLIA,
               markeredgecolor="#333333", markersize=6.5, label="Non-neuronal"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=GLIA,
               markeredgecolor="#111111", markersize=6,
               label="Reactive (−SEAAD), ≥2 pairs"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
               markeredgecolor=GLIA, markersize=6,
               label="Reactive (−SEAAD), 1 pair"),
    ]
    legend_note = (
        "Points are ranked by median ICC across independent "
        "matched-donor pairs. Grey bars are the interquartile range "
        "(or the pair’s ICC confidence interval for single-pair types). "
        "Neuronal points are drawn first so they do not cover "
        "non-neuronal points. −SEAAD names are reactive states in the "
        "SEA-AD DFC 2026 taxonomy. Oligo_2_1-SEAAD is measured only in "
        "the BrainSCOPE CMC × Ruzicka identity-control pair."
    )
    leg = ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.55),
        frameon=False,
        fontsize=7,
        handletextpad=0.4,
        borderaxespad=0.0,
        title=legend_note,
        title_fontsize=6.5,
    )
    if leg.get_title() is not None:
        leg.get_title().set_multialignment("left")

    ax.figure.patch.set_facecolor("white")
    ax.set_facecolor("white")
    png = FIG / "fig7a_supertype_reliability_seaad.png"
    pdf = FIG / "fig7a_supertype_reliability_seaad.pdf"
    ax.figure.savefig(png, dpi=float(cns.settings.savefig_dpi),
                      facecolor="white", bbox_inches="tight", pad_inches=0.04)
    cns.savefig(str(pdf))
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    print(react[["rank", "cell_type", "n_pairs", "median_icc", "single_pair"]]
          .sort_values("rank").to_string(index=False))


if __name__ == "__main__":
    main()
