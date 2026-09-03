#!/usr/bin/env python3
"""
Multi-panel ancestry heterogeneity figure.

Panel A: I² bar chart per lead SNP (total pooled heterogeneity), colour-coded
         by ancestry-level effect concordance (EUR vs AFR direction).
Panel B: Effect-size forest plot — ancestry strata (EUR / AFR / AMR / Pooled)
         for each lead SNP × cell type, with I² displayed.

Follows the style of the 2023 Nature Genetics Parkinson's GWAS (Fig. 3):
  top bar = heterogeneity magnitude, bottom = signed effect direction per group.

Outputs:
  results/meta_sensitivity/ancestry_lead_effects/figures/ancestry_het_figure.svg
  results/meta_sensitivity/ancestry_lead_effects/figures/ancestry_het_figure.png
"""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

# ─── paths ──────────────────────────────────────────────────────────────────
ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
EFF  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

# ─── load data ───────────────────────────────────────────────────────────────
def _f(v):
    try: return float(v)
    except: return None

eff_rows = list(csv.DictReader(open(EFF), delimiter="\t"))
het_rows = list(csv.DictReader(open(HET), delimiter="\t"))

# Build lookup: (cell_type, marker, stratum) → dict
eff = {(r["cell_type"], r["marker"], r["stratum"]): r for r in eff_rows}
het = {(r["cell_type"], r["marker"], r["stratum"]): r for r in het_rows}

# Ordered lead list (unique cell_type+marker, pooled rank order)
seen = set()
leads = []
for r in het_rows:
    k = (r["cell_type"], r["marker"])
    if k not in seen and r["stratum"] == "Pooled":
        seen.add(k)
        leads.append({"cell_type": r["cell_type"], "marker": r["marker"],
                      "i2": _f(r["i2"]), "het_p": _f(r["het_p"]),
                      "beta": _f(r["beta"]), "p": _f(r["p"])})

n = len(leads)  # 19

# Short labels: CellType  chr:pos
def short_label(r):
    parts = r["marker"].split(":")
    loc = f"chr{parts[0].replace('chr','')}:{int(parts[1]):,}"
    return f"{r['cell_type']}\n{loc}"

# Panel A uses short cell-type-only labels; Panel B uses cell_type + locus
a_labels = [r["cell_type"] + ("\n(rank 2)" if r.get("lead_rank","1") == "2"
             else ("\n(rank 3)" if r.get("lead_rank","1") == "3" else ""))
            for r in leads]
labels = [short_label(r) for r in leads]

# ─── colours ─────────────────────────────────────────────────────────────────
# Nature palette (colour-blind safe)
CLR = {
    "EUR":    "#2166AC",  # blue
    "AFR":    "#D95F02",  # orange
    "AMR":    "#1B9E77",  # teal
    "Pooled": "#333333",  # dark grey / diamond
}
HET_LOW  = "#BDBDBD"   # I² < 30 %
HET_MED  = "#F4A460"   # 30–75 %
HET_HIGH = "#B22222"   # ≥ 75 %

def i2_colour(i2):
    if i2 is None: return HET_LOW
    if i2 >= 75:   return HET_HIGH
    if i2 >= 30:   return HET_MED
    return HET_LOW

def i2_category(i2):
    if i2 is None: return "< 30 %"
    if i2 >= 75:   return "≥ 75 %"
    if i2 >= 30:   return "30–75 %"
    return "< 30 %"

# ─── Figure layout ────────────────────────────────────────────────────────────
# Nature full-width ~ 180 mm → 7.1 in.  Two panels stacked: top = I² bars (~1.5in),
# bottom = forest (~5.5in).
fig = plt.figure(figsize=(14, 11))
fig.patch.set_facecolor("white")

gs = gridspec.GridSpec(2, 1, figure=fig,
                       height_ratios=[1.8, 4.5],
                       hspace=0.32,
                       left=0.16, right=0.79, top=0.95, bottom=0.04)

ax_top = fig.add_subplot(gs[0])
ax_bot = fig.add_subplot(gs[1])

x = np.arange(n)
bar_w = 0.65

# ─── Panel A : I² bars ───────────────────────────────────────────────────────
i2_vals  = [r["i2"] if r["i2"] is not None else 0.0 for r in leads]
bar_cols = [i2_colour(r["i2"]) for r in leads]

bars = ax_top.bar(x, i2_vals, width=bar_w, color=bar_cols, linewidth=0, zorder=2)

# 30 % and 75 % guide lines
for yref, ls in [(30, "--"), (75, "-.")]:
    ax_top.axhline(yref, color="black", linewidth=0.6, linestyle=ls, alpha=0.5, zorder=1)
    ax_top.text(n - 0.3, yref + 1.5, f"I²={yref}%", fontsize=6.5,
                va="bottom", ha="right", color="black", alpha=0.6)

# Annotate het P for high-het bars
for i, r in enumerate(leads):
    if r["het_p"] is not None and r["het_p"] < 0.05:
        stars = "***" if r["het_p"] < 0.001 else ("**" if r["het_p"] < 0.01 else "*")
        ax_top.text(i, (r["i2"] or 0) + 1.5, stars,
                    ha="center", va="bottom", fontsize=7, color="black")

ax_top.set_xlim(-0.6, n - 0.4)
ax_top.set_ylim(0, 115)
ax_top.set_xticks(x)
ax_top.set_xticklabels(a_labels, fontsize=6.5, rotation=40, ha="right")
ax_top.set_ylabel("I² (%)\n(pooled)", fontsize=8, labelpad=4)
ax_top.set_title("Ancestry heterogeneity at genome-wide significant lead SNPs",
                 fontsize=10, fontweight="bold", pad=6)
ax_top.tick_params(axis="y", labelsize=7)
ax_top.tick_params(axis="x", bottom=True, length=2)
ax_top.spines["top"].set_visible(False)
ax_top.spines["right"].set_visible(False)
ax_top.text(-0.015, 1.06, "A", transform=ax_top.transAxes,
            fontsize=12, fontweight="bold", va="top")

# Legend for I² colours
leg_patches = [
    mpatches.Patch(color=HET_LOW,  label="I² < 30%"),
    mpatches.Patch(color=HET_MED,  label="I² 30–75%"),
    mpatches.Patch(color=HET_HIGH, label="I² ≥ 75%"),
]
ax_top.legend(handles=leg_patches, fontsize=7, loc="upper left",
              frameon=True, framealpha=0.9, edgecolor="none",
              handlelength=1, handleheight=0.8)

# ─── Panel B : Ancestry forest (row = lead SNP, col = stratum) ──────────────
STRATA = ["EUR", "AFR", "Pooled"]
S_LABEL = {"EUR": "EUR (White)", "AFR": "AFR (Black)", "Pooled": "Pooled (15 cohorts)"}
row_h = 0.9 / n   # fraction of axes height per lead

# Each lead = a horizontal strip; strata plotted as overlapping points
# Y positions: leads along y (top to bottom), strata offset slightly
offsets = {"EUR": 0.18, "AFR": 0.0, "Pooled": -0.18}

y_leads = np.arange(n)[::-1]  # reverse so top lead is highest y

ax_bot.axvline(0, color="grey", linewidth=0.6, linestyle="--", zorder=1)

for i, (lead, yl) in enumerate(zip(leads, y_leads)):
    ct  = lead["cell_type"]
    mk  = lead["marker"]

    # light alternating row shading
    if i % 2 == 0:
        ax_bot.axhspan(yl - 0.5, yl + 0.5, color="#F5F5F5", zorder=0)

    # Add I² label on far right (in its own column before P-value)
    i2_str = f"I²={lead['i2']:.0f}%" if lead["i2"] is not None else "I²=NA"
    ax_bot.text(1.01, yl, i2_str, transform=ax_bot.get_yaxis_transform(),
                ha="left", va="center", fontsize=6,
                color=i2_colour(lead["i2"]))

    for s in STRATA:
        e = eff.get((ct, mk, s))
        if e is None:
            continue
        beta = _f(e["beta"])
        se   = _f(e["se"])
        p    = _f(e["p"])
        if beta is None or se is None:
            continue

        yo = yl + offsets[s]
        ci_lo, ci_hi = beta - 1.96 * se, beta + 1.96 * se
        is_pooled = (s == "Pooled")

        ax_bot.errorbar(beta, yo, xerr=1.96 * se,
                        fmt="D" if is_pooled else "o",
                        color=CLR[s],
                        markersize=5.5 if is_pooled else 3.5,
                        linewidth=1.1 if is_pooled else 0.7,
                        capsize=1.8, capthick=0.7, zorder=3)

        # P annotation for pooled (further right to avoid I² overlap)
        if is_pooled and p is not None:
            p_str = f"P={p:.1e}" if p < 0.001 else f"P={p:.3f}"
            ax_bot.text(1.12, yl, p_str, transform=ax_bot.get_yaxis_transform(),
                        ha="left", va="center", fontsize=6, color="black")

ax_bot.set_ylim(-0.8, n - 0.2)
ax_bot.set_yticks(y_leads)
ax_bot.set_yticklabels(labels, fontsize=7)
ax_bot.set_xlabel("Beta (95% CI)\nharmonized to pooled effect allele", fontsize=8)
ax_bot.tick_params(axis="x", labelsize=7)
ax_bot.spines["top"].set_visible(False)
ax_bot.spines["right"].set_visible(False)
ax_bot.text(-0.15, 1.02, "B", transform=ax_bot.transAxes,
            fontsize=12, fontweight="bold", va="top")

# Legend for strata
leg_elems = [
    Line2D([0],[0], marker="o", color=CLR["EUR"],    linestyle="none", markersize=5,
           label="EUR meta (White, N≈2,818)"),
    Line2D([0],[0], marker="o", color=CLR["AFR"],    linestyle="none", markersize=5,
           label="AFR meta (African American, N≈270)"),
    Line2D([0],[0], marker="D", color=CLR["Pooled"], linestyle="none", markersize=6,
           label="Pooled 15-cohort meta"),
]
ax_bot.legend(handles=leg_elems, fontsize=7, loc="lower right",
              frameon=True, framealpha=0.9, edgecolor="none")

# ─── Sub-title note on AMR ───────────────────────────────────────────────────
fig.text(0.97, 0.005,
         "Latino/AMR (N≈178, single cohort) omitted: lead SNPs absent from Mayo AMR sumstats for 18/19 loci.",
         ha="right", va="bottom", fontsize=6, color="grey", style="italic")

# ─── Export ──────────────────────────────────────────────────────────────────
for ext in ["png", "svg", "pdf"]:
    out = OUTD / f"ancestry_het_figure.{ext}"
    fig.savefig(out, dpi=300 if ext == "png" else 150, bbox_inches="tight",
                facecolor="white")
    print(f"Saved {out}")

plt.close(fig)
