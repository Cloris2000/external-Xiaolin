#!/usr/bin/env python3
"""
Supplementary figure: meta-analysis robustness for the 19 lead SNPs.

  A : fixed-effect (IVW) vs random-effects (DL) pooled β ± 95% CI per lead,
      with the leave-one-cohort-out (LOO) β range shaded behind.
  B : −log10 P for fixed-effect, random-effects and the worst LOO leave-out,
      relative to the genome-wide line (5e-8).

Reads lead_meta_sensitivity.tsv (produced by lead_meta_sensitivity.py).
Outputs: results/.../figures/supp_meta_sensitivity.{png,svg,pdf}
"""
import csv
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

ROOT = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
TSV  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/lead_meta_sensitivity.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

def f(v):
    try: return float(v)
    except: return None

def short_ct(ct):
    return ct.replace('Oligodendrocyte', 'Oligo')

def locus_label(marker):
    p = marker.split(':')
    return f"{p[0].replace('chr','')}:{int(p[1])/1e6:.1f}M"

rows = list(csv.DictReader(open(TSV), delimiter='\t'))
rows = rows[::-1]  # so file order reads top→bottom on the y-axis
n = len(rows)
y = np.arange(n)

FE_CLR, RE_CLR, LOO_CLR = '#1f6fb4', '#c1442e', '#C9C9C9'
GW = 5e-8
GW_LINE = -math.log10(GW)

fig = plt.figure(figsize=(12, 9))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.6, 1.0], wspace=0.06,
                       left=0.26, right=0.955, top=0.94, bottom=0.08)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1], sharey=ax_a)

# ── Panel A : FE vs RE effect sizes with LOO range ─────────────────────────────
ax_a.axvline(0, color='grey', lw=0.9, ls='--', zorder=1)
for i, r in enumerate(rows):
    fb, fse = f(r['fe_beta']), f(r['fe_se'])
    rb, rse = f(r['re_beta']), f(r['re_se'])
    lmin, lmax = f(r['loo_beta_min']), f(r['loo_beta_max'])
    # LOO range bar
    ax_a.add_patch(plt.Rectangle((lmin, y[i]-0.30), lmax-lmin, 0.60,
                                 color=LOO_CLR, alpha=0.55, lw=0, zorder=1))
    # fixed effect (upper)
    ax_a.errorbar(fb, y[i]+0.16, xerr=1.96*fse, fmt='o', color=FE_CLR,
                  markersize=5.5, elinewidth=1.6, capsize=2.4, zorder=4)
    # random effect (lower)
    ax_a.errorbar(rb, y[i]-0.16, xerr=1.96*rse, fmt='D', color=RE_CLR,
                  markersize=5.0, elinewidth=1.6, capsize=2.4, zorder=4)

labels = [f"{short_ct(r['cell_type'])}  ({locus_label(r['marker'])})" for r in rows]
ax_a.set_yticks(y)
ax_a.set_yticklabels(labels, fontsize=8.5)
ax_a.set_ylim(-0.6, n-0.4)
ax_a.set_xlabel("Pooled effect size (β, 95% CI)", fontsize=9.5)
ax_a.tick_params(axis='x', labelsize=8.5)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(-0.34, 1.015, 'A', transform=ax_a.transAxes,
          fontsize=15, fontweight='bold', va='bottom', ha='left')
ax_a.legend(handles=[
    Line2D([0],[0], marker='o', color=FE_CLR, ls='none', markersize=6, label='Fixed-effect (IVW)'),
    Line2D([0],[0], marker='D', color=RE_CLR, ls='none', markersize=6, label='Random-effects (DL)'),
    plt.Rectangle((0,0),1,1, color=LOO_CLR, alpha=0.55, label='Leave-one-out β range'),
], fontsize=8, loc='lower left', frameon=True, framealpha=0.9, edgecolor='none')

# ── Panel B : significance robustness (−log10 P) ───────────────────────────────
def nlp(p):
    p = f(p)
    if p is None or p <= 0: return 0.0
    return min(-math.log10(p), 70)

for i, r in enumerate(rows):
    fp, rp, lp = nlp(r['fe_p']), nlp(r['re_p']), nlp(r['loo_p_max'])
    ax_b.plot([lp, fp], [y[i], y[i]], color='#DDDDDD', lw=1.2, zorder=1)
    ax_b.scatter(fp, y[i], s=34, color=FE_CLR, zorder=3)
    ax_b.scatter(rp, y[i], s=30, marker='D', color=RE_CLR, zorder=3)
    ax_b.scatter(lp, y[i], s=26, marker='|', color='#555555', linewidths=1.6, zorder=3)

ax_b.axvline(GW_LINE, color='#B22222', lw=1.1, ls=':', zorder=2)
ax_b.text(GW_LINE, n-0.2, ' 5×10⁻⁸', color='#B22222', fontsize=7.5,
          va='top', ha='left', rotation=90)
ax_b.set_xlabel("−log₁₀ P", fontsize=9.5)
ax_b.tick_params(axis='x', labelsize=8.5)
ax_b.tick_params(axis='y', left=False, labelleft=False)
ax_b.set_xlim(left=0)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.02, 1.015, 'B', transform=ax_b.transAxes,
          fontsize=15, fontweight='bold', va='bottom', ha='left')
ax_b.legend(handles=[
    Line2D([0],[0], marker='o', color=FE_CLR, ls='none', markersize=6, label='Fixed-effect'),
    Line2D([0],[0], marker='D', color=RE_CLR, ls='none', markersize=6, label='Random-effects'),
    Line2D([0],[0], marker='|', color='#555555', ls='none', markersize=9, label='Worst leave-one-out'),
], fontsize=8, loc='lower right', frameon=True, framealpha=0.9, edgecolor='none')

for ext in ['png', 'svg', 'pdf']:
    out = OUTD / f"supp_meta_sensitivity.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")
plt.close(fig)
