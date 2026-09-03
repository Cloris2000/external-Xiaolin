#!/usr/bin/env python3
"""
SUPPLEMENTARY cross-ANCESTRY heterogeneity figure — the complete 19-locus
version. The compact main-paper version now lives in
plot_main_heterogeneity_figure.py (Panels A–C).

  A (top)    : across-ancestry I² bar chart per lead SNP (shared x with B).
  B (bottom) : per-lead forest — Meta / EUR / AFR / LAT-AMR shown as points with
               vertical 95% CI whiskers (y = effect size β).

Outputs: results/.../figures/supp_cross_ancestry_forest.{png,svg,pdf}
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
EFF  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
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

# ── load ──────────────────────────────────────────────────────────────────────
eff_rows = {(r['cell_type'], r['marker'], r['stratum']): r
            for r in csv.DictReader(open(EFF), delimiter='\t')}
het_all  = list(csv.DictReader(open(HET), delimiter='\t'))

seen, leads = set(), []
for r in het_all:
    k = (r['cell_type'], r['marker'])
    if k not in seen and r['stratum'] == 'Pooled':
        seen.add(k)
        leads.append({'cell_type': r['cell_type'], 'marker': r['marker']})
N = len(leads)
x = np.arange(N)

# ── colours ───────────────────────────────────────────────────────────────────
HET_LOW, HET_MED, HET_HIGH = "#BDBDBD", "#F4A460", "#B22222"
def i2_colour(v):
    if v is None: return HET_LOW
    if v >= 75:   return HET_HIGH
    if v >= 30:   return HET_MED
    return HET_LOW

STRATA      = ['Pooled', 'EUR', 'AFR', 'AMR']
STRATA_LBL  = {'Pooled': 'Meta', 'EUR': 'EUR', 'AFR': 'AFR', 'AMR': 'LAT/AMR'}
STRATA_CLR  = {'Pooled': '#333333', 'EUR': '#2166AC', 'AFR': '#D95F02', 'AMR': '#1B9E77'}
STRATA_OFF  = {'Pooled': -0.27, 'EUR': -0.09, 'AFR': 0.09, 'AMR': 0.27}

# ── across-ancestry I² ────────────────────────────────────────────────────────
def across_ancestry_i2(ct, mk):
    est = []
    for s in ['EUR', 'AFR', 'AMR']:
        e = eff_rows.get((ct, mk, s))
        if e is None: continue
        b = f(e['beta']); se = f(e['se'])
        if b is not None and se is not None and se > 0:
            est.append((b, se))
    k = len(est)
    if k < 2:
        return None
    w = [1/se**2 for _, se in est]
    bbar = sum(wi*bi for (bi, _), wi in zip(est, w)) / sum(w)
    Q = sum(wi*(bi-bbar)**2 for (bi, _), wi in zip(est, w))
    return max(0.0, (Q-(k-1))/Q*100) if Q > 0 else 0.0

anc_i2 = [across_ancestry_i2(l['cell_type'], l['marker']) for l in leads]
col_labels = [f"{short_ct(l['cell_type'])}\n{locus_label(l['marker'])}" for l in leads]

# ── figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 8))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.5, 2.4], hspace=0.10,
                       left=0.085, right=0.90, top=0.95, bottom=0.20)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1], sharex=ax_a)

# Panel A
for j, i2v in enumerate(anc_i2):
    if i2v is None:
        continue
    ax_a.bar(j, i2v, width=0.7, color=i2_colour(i2v), linewidth=0, zorder=2)
ax_a.set_xlim(-0.6, N - 0.4)
ax_a.set_ylim(0, 105)
ax_a.set_ylabel("Across-ancestry\nheterogeneity (I², %)", fontsize=9, labelpad=6)
ax_a.tick_params(axis='y', labelsize=8)
ax_a.tick_params(axis='x', bottom=False, labelbottom=False)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(-0.075, 1.02, 'A', transform=ax_a.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')
ax_a.legend(handles=[mpatches.Patch(color=HET_LOW, label='I² < 30%'),
                     mpatches.Patch(color=HET_MED, label='30–75%'),
                     mpatches.Patch(color=HET_HIGH, label='≥ 75%')],
            fontsize=8, loc='upper right', frameon=True, framealpha=0.9,
            edgecolor='none', handlelength=1.0)

# Panel B — per-lead ancestry forest
ax_b.axhline(0, color='grey', lw=0.8, ls='--', zorder=1)
all_lo, all_hi = [], []
for j, lead in enumerate(leads):
    for s in STRATA:
        e = eff_rows.get((lead['cell_type'], lead['marker'], s))
        if e is None: continue
        b = f(e['beta'])
        if b is None: continue
        lo = f(e.get('ci_lo')); hi = f(e.get('ci_hi'))
        if lo is None or hi is None:
            se = f(e['se'])
            if se is None: continue
            lo, hi = b - 1.96*se, b + 1.96*se
        ax_b.errorbar(j + STRATA_OFF[s], b, yerr=[[b-lo], [hi-b]],
                      fmt='D' if s == 'Pooled' else 'o', color=STRATA_CLR[s],
                      markersize=4.8 if s == 'Pooled' else 4.0,
                      markeredgewidth=0, elinewidth=1.0, capsize=0, zorder=3)
        all_lo.append(lo); all_hi.append(hi)
for j in range(N - 1):
    ax_b.axvline(j + 0.5, color='#EEEEEE', lw=0.5, zorder=0)
lo_b, hi_b = min(all_lo), max(all_hi)
pad = 0.06 * (hi_b - lo_b)
ax_b.set_xlim(-0.6, N - 0.4)
ax_b.set_ylim(lo_b - pad, hi_b + pad)
ax_b.set_ylabel("Effect size (β, 95% CI)", fontsize=9, labelpad=6)
ax_b.tick_params(axis='y', labelsize=8)
ax_b.set_xticks(x)
ax_b.set_xticklabels(col_labels, fontsize=8, rotation=45, ha='right')
ax_b.tick_params(axis='x', length=3)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.075, 1.03, 'B', transform=ax_b.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')
ax_b.legend(handles=[Line2D([0],[0], marker='D' if s == 'Pooled' else 'o',
                            color=STRATA_CLR[s], ls='none', markersize=6,
                            label=STRATA_LBL[s]) for s in STRATA],
            fontsize=8, loc='upper left', bbox_to_anchor=(1.005, 1.0),
            frameon=True, framealpha=0.9, edgecolor='none',
            title='Ancestry', title_fontsize=8)

for ext in ['png', 'svg', 'pdf']:
    out = OUTD / f"supp_cross_ancestry_forest.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")
plt.close(fig)
