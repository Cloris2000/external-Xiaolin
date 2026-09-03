#!/usr/bin/env python3
"""Figure 5: bulk ↔ sn CTP GWAS concordance (cnsplots styling).

Panel A — direction concordance among loci found in sn (same vs opposite)
Panel B — bulk β vs sn β scatter for all shared suggestive loci

Supplementary — forest of top bulk lead found in sn, per cell type
(previous Panel B; kept for locus-level detail).

Panels A/B use different journal color cycles (Science vs Cell).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

import cnsplots as cns

ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
OUT_DIR = ROOT / "manuscript_figure"
HITS = ROOT / "results/sn_bulk_meta_similarity_hodge3/top_hits/bulk_suggestive_hits_sn_direction.tsv"
CONC = ROOT / "results/sn_bulk_meta_similarity_hodge3/top_hits/cell_type_concordance_summary.tsv"

# Science — Panel A (same vs opposite bars)
SCIENCE = [
    "#3B4992",  # blue
    "#EE0000",  # red
    "#008B45",
    "#631879",
    "#008280",
    "#BB0021",
]
# Cell — Panel B / supp (distinct from Science)
CELL = [
    "#C84C3A",  # brick
    "#2F7E8F",  # steel teal
    "#E1A22E",  # gold
    "#4E5A8A",  # indigo
    "#5F9862",  # green
]


def gene_label(chrom: str, pos: int) -> str:
    ch = str(chrom).removeprefix("chr")
    if ch == "7" and 12_000_000 <= pos <= 12_800_000:
        return "TMEM106B"
    if ch == "12" and 2_000_000 <= pos <= 2_800_000:
        return "CACNA1C"
    if ch == "6" and 161_000_000 <= pos <= 163_000_000:
        return "PRKN"
    if ch == "2" and 10_000_000 <= pos <= 11_500_000:
        return "LINC01954"
    if ch == "3" and 142_000_000 <= pos <= 143_000_000:
        return "PLS1"
    if ch == "8" and 57_000_000 <= pos <= 58_000_000:
        return "RPL30P10"
    if ch == "9" and 130_000_000 <= pos <= 131_000_000:
        return "HMCN2"
    return f"{ch}:{pos}"


def _style_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.tick_params(axis="both", colors="#111111", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#333333")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _opaque_save(fig: plt.Figure, out_png: Path, out_pdf: Path, out_svg: Path) -> None:
    fig.patch.set_facecolor("white")
    fig.patch.set_alpha(1.0)
    for a in fig.axes:
        a.set_facecolor("white")
        a.patch.set_alpha(1.0)
    fig.savefig(
        out_png,
        dpi=float(cns.settings.savefig_dpi),
        facecolor="white",
        edgecolor="none",
        transparent=False,
        bbox_inches="tight",
        pad_inches=0.05,
    )
    cns.savefig(str(out_pdf))
    cns.savefig(str(out_svg))


def load_panel_a() -> tuple[pd.DataFrame, list[str]]:
    conc = pd.read_csv(CONC, sep="\t")
    conc["pct_of_found"] = np.where(
        conc["n_found_sn"] > 0,
        np.round(100 * conc["n_concordant"] / conc["n_found_sn"]).astype(int),
        np.nan,
    )
    conc = conc.sort_values("pct_of_found", ascending=False).reset_index(drop=True)
    return conc, conc["bulk_cell_type"].tolist()


def load_shared_loci() -> pd.DataFrame:
    hits = pd.read_csv(HITS, sep="\t")
    found = hits[hits["sn_found"] == True].copy()  # noqa: E712
    found["same_sign"] = (found["bulk_effect"] * found["sn_effect"]) > 0
    found["neglog10_bulk_p"] = -np.log10(found["bulk_p"].clip(lower=1e-300))
    return found


def load_top_leads(ct_order: list[str]) -> pd.DataFrame:
    found = load_shared_loci()
    top = (
        found.sort_values("bulk_p")
        .groupby("bulk_cell_type", as_index=False)
        .first()
    )
    top["cell_type"] = pd.Categorical(
        top["bulk_cell_type"], categories=ct_order, ordered=True
    )
    top = top.sort_values("cell_type").dropna(subset=["cell_type"]).reset_index(drop=True)
    top["gene"] = [
        gene_label(c, int(p)) for c, p in zip(top["chrom"], top["pos"])
    ]
    top["row_label"] = top["cell_type"].astype(str) + "  (" + top["gene"] + ")"
    top["bulk_lo"] = top["bulk_effect"] - 1.96 * top["bulk_stderr"]
    top["bulk_hi"] = top["bulk_effect"] + 1.96 * top["bulk_stderr"]
    top["sn_lo"] = top["sn_effect"] - 1.96 * top["sn_stderr"]
    top["sn_hi"] = top["sn_effect"] + 1.96 * top["sn_stderr"]
    return top


def plot_panel_a(ax: plt.Axes, conc: pd.DataFrame, ct_order: list[str]) -> None:
    """Horizontal stacked bars: same (Science blue) vs opposite (Science red)."""
    same_c, opp_c = SCIENCE[0], SCIENCE[1]
    order = list(reversed(ct_order))
    y = np.arange(len(order))
    same = conc.set_index("bulk_cell_type").loc[order, "n_concordant"].to_numpy()
    opp = conc.set_index("bulk_cell_type").loc[order, "n_discordant"].to_numpy()
    pct = conc.set_index("bulk_cell_type").loc[order, "pct_of_found"].to_numpy()

    _style_axes(ax)
    ax.barh(y, same, color=same_c, height=0.7, label="Same direction", zorder=2)
    ax.barh(y, opp, left=same, color=opp_c, height=0.7, label="Opposite", zorder=2)
    for yi, s, o, p in zip(y, same, opp, pct):
        ax.text(
            s + o + 0.4, yi, f"{int(p)}%", va="center", ha="left",
            fontsize=7, color="#222222",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8, color="#111111")
    ax.set_xlabel("Number of loci", fontsize=9, color="#111111")
    ax.set_xlim(0, max(same + opp) * 1.18)
    leg = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, -0.14),
        ncol=2, frameon=False, fontsize=8,
    )
    for t in leg.get_texts():
        t.set_color("#111111")


def plot_panel_b_scatter(ax: plt.Axes, shared: pd.DataFrame) -> None:
    """Bulk β vs sn β for all shared loci (Cell palette: teal=same, gold=opposite)."""
    same_c, opp_c = CELL[1], CELL[2]
    _style_axes(ax)

    # size ~ -log10(bulk P), clipped for readability
    sizes = 18 + 6 * np.clip(shared["neglog10_bulk_p"] - 5, 0, 8)

    same = shared[shared["same_sign"]]
    opp = shared[~shared["same_sign"]]
    ax.scatter(
        same["bulk_effect"], same["sn_effect"],
        s=sizes[shared["same_sign"].to_numpy()],
        c=same_c, alpha=0.75, edgecolors="none", zorder=3, label="Same direction",
    )
    ax.scatter(
        opp["bulk_effect"], opp["sn_effect"],
        s=sizes[~shared["same_sign"].to_numpy()],
        c=opp_c, alpha=0.80, edgecolors="none", zorder=3, label="Opposite",
    )

    # axes + 1:1 guide
    lim = 0.42
    ax.axhline(0, color="#BBBBBB", linewidth=0.7, zorder=1)
    ax.axvline(0, color="#BBBBBB", linewidth=0.7, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], color="#888888", linestyle="--",
            linewidth=0.8, zorder=1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Bulk CTP GWAS $\beta$", fontsize=9, color="#111111")
    ax.set_ylabel(r"sn CTP GWAS $\beta$", fontsize=9, color="#111111")

    r, p = stats.spearmanr(shared["bulk_effect"], shared["sn_effect"])
    pct = 100 * shared["same_sign"].mean()
    ax.text(
        0.03, 0.97,
        f"n = {len(shared)}\n"
        f"Same direction = {pct:.0f}%\n"
        f"Spearman r = {r:.2f}\n"
        r"Size $\propto$ $-$log$_{10}P_{\mathrm{bulk}}$",
        transform=ax.transAxes, va="top", ha="left",
        fontsize=8, color="#111111",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2),
    )

    leg = ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, -0.18),
        ncol=2, frameon=False, fontsize=8,
        markerscale=1.2,
    )
    for t in leg.get_texts():
        t.set_color("#111111")


def plot_supp_forest(ax: plt.Axes, top: pd.DataFrame) -> None:
    """Supplementary: top bulk lead found in sn, per CT (Cell teal/gold)."""
    bulk_c, sn_c = CELL[1], CELL[2]
    labels = top["row_label"].tolist()
    y = np.arange(len(labels))[::-1]

    _style_axes(ax)
    for yi, row in zip(y, top.itertuples(index=False)):
        ax.errorbar(
            row.bulk_effect, yi + 0.15,
            xerr=[[row.bulk_effect - row.bulk_lo], [row.bulk_hi - row.bulk_effect]],
            fmt="o", color=bulk_c, markersize=4.5, elinewidth=1.0,
            capsize=2, zorder=3,
        )
        ax.errorbar(
            row.sn_effect, yi - 0.15,
            xerr=[[row.sn_effect - row.sn_lo], [row.sn_hi - row.sn_effect]],
            fmt="^", color=sn_c, markersize=4.5, elinewidth=1.0,
            capsize=2, zorder=3,
        )

    ax.axvline(0, color="#666666", linestyle="--", linewidth=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5, color="#111111")
    ax.set_xlabel(r"$\beta$ (95% CI)", fontsize=9, color="#111111")
    ax.set_xlim(-0.55, 0.45)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.6, zorder=0)

    handles = [
        Line2D([0], [0], marker="o", color=bulk_c, linestyle="None",
               markersize=6, label="Bulk CTP GWAS"),
        Line2D([0], [0], marker="^", color=sn_c, linestyle="None",
               markersize=6, label="sn CTP GWAS"),
    ]
    leg = ax.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.10),
        ncol=2, frameon=False, fontsize=8,
    )
    for t in leg.get_texts():
        t.set_color("#111111")


def main() -> None:
    conc, ct_order = load_panel_a()
    shared = load_shared_loci()
    top = load_top_leads(ct_order)

    # ---- Main Figure 5 ----
    mp = cns.multipanel(max_width=440)

    ax_a = mp.panel(
        "A", width=420, height=340,
        color_cycle="Science",
        margin_bottom=30,
    )
    plot_panel_a(ax_a, conc, ct_order)

    ax_b = mp.panel(
        "B", width=360, height=360,
        color_cycle="Cell",
        below="A",
        margin_top=24,
        margin_bottom=36,
    )
    plot_panel_b_scatter(ax_b, shared)

    out_png = OUT_DIR / "figure5_sn_concordance.png"
    out_pdf = OUT_DIR / "figure5_sn_concordance.pdf"
    out_svg = OUT_DIR / "figure5_sn_concordance.svg"
    _opaque_save(ax_a.figure, out_png, out_pdf, out_svg)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")
    print(f"Wrote {out_svg}")

    # ---- Supplementary forest (previous Panel B) ----
    plt.close("all")
    cns.figure(width=420, height=440, color_cycle="Cell")
    ax_s = plt.gca()
    plot_supp_forest(ax_s, top)

    supp_png = OUT_DIR / "figure5_sn_concordance_supp_forest.png"
    supp_pdf = OUT_DIR / "figure5_sn_concordance_supp_forest.pdf"
    supp_svg = OUT_DIR / "figure5_sn_concordance_supp_forest.svg"
    _opaque_save(ax_s.figure, supp_png, supp_pdf, supp_svg)
    print(f"Wrote {supp_png}")
    print(f"Wrote {supp_pdf}")
    print(f"Wrote {supp_svg}")

    r, _ = stats.spearmanr(shared["bulk_effect"], shared["sn_effect"])
    print(
        f"Shared loci n={len(shared)}; "
        f"same-sign={100 * shared['same_sign'].mean():.1f}%; "
        f"Spearman r={r:.3f}"
    )
    print("Panel A (Science): same=", SCIENCE[0], "opposite=", SCIENCE[1])
    print("Panel B (Cell):    same=", CELL[1], "opposite=", CELL[2])


if __name__ == "__main__":
    main()
