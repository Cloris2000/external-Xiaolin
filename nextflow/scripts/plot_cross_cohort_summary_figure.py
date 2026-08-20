#!/usr/bin/env python3
"""
Cross-COHORT heterogeneity — compact MAIN-figure version, matched in layout to
the cross-ancestry figure (plot_cross_ancestry_figure.py).

  A (top)    : across-cohort I² bar chart per lead SNP (METAL pooled 15-cohort).
  B (bottom) : per-lead scatter of every contributing cohort's β (small points)
               around the pooled meta estimate (black diamond ± 95% CI).
               y = effect size β; shared x with A.

The full 19-panel per-cohort forest grid (plot_cross_cohort_figure.py) remains
the supplementary version.

Outputs: results/.../figures/figure6_cross_cohort.{png,svg,pdf}  (main)
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

ROOT = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
EFF  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
COH  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_cohort_effects.tsv"
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

# leads in file order (SAME order as the cross-ancestry figure → matched pair)
seen, leads = set(), []
for r in het_all:
    k = (r['cell_type'], r['marker'])
    if k not in seen and r['stratum'] == 'Pooled':
        seen.add(k)
        leads.append({'cell_type': r['cell_type'], 'marker': r['marker'],
                      'i2': f(r['i2']), 'het_p': f(r['het_p'])})
N = len(leads)
x = np.arange(N)

coh_full = {}
for r in csv.DictReader(open(COH), delimiter='\t'):
    if r['beta'] and r['status'] in ('ok', 'flipped'):
        coh_full.setdefault((r['cell_type'], r['marker']), []).append(
            (r['cohort'], f(r['beta'])))

# ── cohort identity (colour + label per cohort) ───────────────────────────────
EUR_ORDER = ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN',
             'CMC_PITT','GTEx_v10','NABEC','GVEX']
AFR_ORDER = ['NIMH_HBCC_1M','NIMH_HBCC_h650','NIMH_HBCC_Omni5M','AMP_AD_Rush']
AMR_ORDER = ['AMP_AD_Mayo']
COHORT_ORDER = EUR_ORDER + AFR_ORDER + AMR_ORDER
COHORT_RANK  = {c: i for i, c in enumerate(COHORT_ORDER)}
COHORT_LABEL = {
    'ROSMAP':'ROSMAP','ROSMAP_array':'ROSMAP-array','Mayo':'Mayo','MSBB':'MSBB',
    'CMC_MSSM':'CMC-MSSM','CMC_PENN':'CMC-PENN','CMC_PITT':'CMC-PITT',
    'GTEx_v10':'GTEx v10','NABEC':'NABEC','GVEX':'GVEX',
    'NIMH_HBCC_1M':'HBCC-1M','NIMH_HBCC_h650':'HBCC-h650','NIMH_HBCC_Omni5M':'HBCC-5M',
    'AMP_AD_Rush':'Rush','AMP_AD_Mayo':'Mayo-LAT',
}
_tab = plt.get_cmap('tab20')
COHORT_COLORS = {c: _tab(i % 20) for i, c in enumerate(COHORT_ORDER)}

# ── colours ───────────────────────────────────────────────────────────────────
HET_LOW, HET_MED, HET_HIGH = "#BDBDBD", "#F4A460", "#B22222"
def i2_colour(v):
    if v is None: return HET_LOW
    if v >= 75:   return HET_HIGH
    if v >= 30:   return HET_MED
    return HET_LOW
META_CLR = '#111111'

col_labels = [f"{short_ct(l['cell_type'])}\n{locus_label(l['marker'])}" for l in leads]

# ── figure (matched to plot_cross_ancestry_figure.py) ─────────────────────────
fig = plt.figure(figsize=(13.5, 8))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.5, 2.4], hspace=0.10,
                       left=0.085, right=0.87, top=0.95, bottom=0.20)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1], sharex=ax_a)

# Panel A — across-cohort I² bars
for j, lead in enumerate(leads):
    i2v = lead['i2']
    if i2v is None:
        continue
    ax_a.bar(j, i2v, width=0.7, color=i2_colour(i2v), linewidth=0, zorder=2)
ax_a.set_xlim(-0.6, N - 0.4)
ax_a.set_ylim(0, 105)
ax_a.set_ylabel("Across-cohort\nheterogeneity (I², %)", fontsize=9, labelpad=6)
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

# Panel B — per-lead cohort scatter + pooled meta diamond
ax_b.axhline(0, color='grey', lw=0.8, ls='--', zorder=1)
for j in range(N - 1):
    ax_b.axvline(j + 0.5, color='#EEEEEE', lw=0.5, zorder=0)

all_vals = []
cohorts_present = set()
for j, lead in enumerate(leads):
    ct, mk = lead['cell_type'], lead['marker']
    pairs = [(c, b) for c, b in coh_full.get((ct, mk), []) if b is not None]
    # sort by β so points spread cleanly along the column
    pairs.sort(key=lambda cb: cb[1])
    n = len(pairs)
    if n:
        offs = np.linspace(-0.24, 0.24, n) if n > 1 else np.array([0.0])
        for k, (coh, b) in enumerate(pairs):
            ax_b.scatter(j + offs[k], b, s=26,
                         color=COHORT_COLORS.get(coh, '#888888'), alpha=0.9,
                         edgecolors='white', linewidths=0.4, zorder=3)
            cohorts_present.add(coh)
        all_vals.extend(b for _, b in pairs)
    # pooled meta diamond + 95% CI
    e = eff_rows.get((ct, mk, 'Pooled'))
    if e is not None:
        b = f(e['beta']); lo = f(e.get('ci_lo')); hi = f(e.get('ci_hi'))
        if b is not None:
            if lo is None or hi is None:
                se = f(e['se']); lo, hi = (b-1.96*se, b+1.96*se) if se else (b, b)
            ax_b.errorbar(j, b, yerr=[[b-lo], [hi-b]], fmt='D', color=META_CLR,
                          markersize=5.5, markeredgecolor='white', markeredgewidth=0.6,
                          elinewidth=1.3, capsize=0, zorder=5)
            all_vals.extend([lo, hi])

lo_b, hi_b = min(all_vals), max(all_vals)
pad = 0.06 * (hi_b - lo_b)
ax_b.set_xlim(-0.6, N - 0.4)
ax_b.set_ylim(lo_b - pad, hi_b + pad)
ax_b.set_ylabel("Effect size (β)", fontsize=9, labelpad=6)
ax_b.tick_params(axis='y', labelsize=8)
ax_b.set_xticks(x)
ax_b.set_xticklabels(col_labels, fontsize=8, rotation=45, ha='right')
ax_b.tick_params(axis='x', length=3)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.075, 1.03, 'B', transform=ax_b.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')

# legend: one entry per cohort (fixed ancestry-grouped order) + meta diamond
coh_handles = [
    Line2D([0], [0], marker='o', color=COHORT_COLORS[c], ls='none', markersize=6,
           markeredgecolor='white', markeredgewidth=0.4, label=COHORT_LABEL[c])
    for c in COHORT_ORDER if c in cohorts_present
]
meta_handle = Line2D([0], [0], marker='D', color=META_CLR, ls='none', markersize=6,
                     markeredgecolor='white', label='Meta (95% CI)')
ax_b.legend(handles=coh_handles + [meta_handle], fontsize=7.2, loc='upper left',
            bbox_to_anchor=(1.005, 1.02), frameon=True, framealpha=0.9,
            edgecolor='none', handletextpad=0.4, labelspacing=0.35,
            title='Cohort', title_fontsize=8)

for ext in ['png', 'svg', 'pdf']:
    out = OUTD / f"figure6_cross_cohort.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")
plt.close(fig)
