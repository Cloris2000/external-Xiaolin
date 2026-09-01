#!/usr/bin/env python3
"""
MAIN heterogeneity figure — consolidates the old cross-ancestry (Figure 5) and
cross-cohort (Figure 6) figures into a single progressive story:

  A : global comparison of across-cohort vs across-ancestry I² for every lead
      (each point = one lead SNP; the 3 representative loci are highlighted).
  B : detailed ancestry + cohort forest for 3 representative loci
      (overall meta, EUR/AFR subtotals, and every contributing cohort).
  C : leave-one-cohort-out estimates for the same 3 loci, showing no single
      cohort drives the association.

The complete 19-locus ancestry and cohort plots move to supplementary figures
(plot_cross_ancestry_figure.py / plot_cross_cohort_figure.py).

Outputs: results/.../figures/figure5_heterogeneity_main.{png,svg,pdf}
"""
import csv
import math
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

ROOT = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
EFF  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
COH  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_cohort_effects.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

SQRT2 = math.sqrt(2.0)

def f(v):
    try: return float(v)
    except (TypeError, ValueError): return None

def short_ct(ct):
    return ct.replace('Oligodendrocyte', 'Oligo').replace('.', '-')

def locus_label(marker):
    p = marker.split(':')
    return f"{p[0].replace('chr','')}:{int(p[1])/1e6:.1f}M"

COHORT_LABEL = {
    'ROSMAP':'ROSMAP','ROSMAP_array':'ROSMAP-array','Mayo':'Mayo','MSBB':'MSBB',
    'CMC_MSSM':'CMC-MSSM','CMC_PENN':'CMC-PENN','CMC_PITT':'CMC-PITT',
    'GTEx_v10':'GTEx v10','NABEC':'NABEC','GVEX':'GVEX',
    'NIMH_HBCC_1M':'HBCC-1M','NIMH_HBCC_h650':'HBCC-h650','NIMH_HBCC_Omni5M':'HBCC-5M',
    'AMP_AD_Rush':'Rush','AMP_AD_Mayo':'Mayo-LAT',
}
COHORT_ORDER = ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN','CMC_PITT',
                'GTEx_v10','NABEC','GVEX','NIMH_HBCC_1M','NIMH_HBCC_h650',
                'NIMH_HBCC_Omni5M','AMP_AD_Rush','AMP_AD_Mayo']
COHORT_RANK = {c: i for i, c in enumerate(COHORT_ORDER)}

EUR_CLR   = '#2166AC'
AFR_CLR   = '#D95F02'
META_CLR  = '#111111'
EUR_LT    = '#7FA8D0'
AFR_LT    = '#E6A366'
LOO_CLR   = '#4C4C4C'

# ── representative loci (low → moderate → high I²) ─────────────────────────────
REPS = [
    ('Microglia', 'chr16:31298939:T:G'),
    ('L6b',       'chr3:142612091:C:T'),
    ('L5.ET',     'chr2:10862188:G:A'),
]

# ── load ──────────────────────────────────────────────────────────────────────
eff_rows = {(r['cell_type'], r['marker'], r['stratum']): r
            for r in csv.DictReader(open(EFF), delimiter='\t')}

het_all = list(csv.DictReader(open(HET), delimiter='\t'))
leads, seen = [], set()
pooled_i2 = {}
for r in het_all:
    k = (r['cell_type'], r['marker'])
    if r['stratum'] == 'Pooled':
        pooled_i2[k] = f(r['i2'])
        if k not in seen:
            seen.add(k)
            leads.append({'cell_type': r['cell_type'], 'marker': r['marker']})

coh = defaultdict(list)  # (ct, marker) -> [(cohort, ancestry, beta, se)]
for r in csv.DictReader(open(COH), delimiter='\t'):
    if r['status'] in ('ok', 'flipped'):
        b, s = f(r['beta']), f(r['se'])
        if b is not None and s is not None and s > 0:
            coh[(r['cell_type'], r['marker'])].append((r['cohort'], r['ancestry'], b, s))

def across_ancestry_i2(ct, mk):
    est = []
    for s in ['EUR', 'AFR', 'AMR']:
        e = eff_rows.get((ct, mk, s))
        if e is None: continue
        b, se = f(e['beta']), f(e['se'])
        if b is not None and se and se > 0:
            est.append((b, se))
    k = len(est)
    if k < 2: return None
    w = [1/se**2 for _, se in est]
    bbar = sum(wi*bi for (bi, _), wi in zip(est, w)) / sum(w)
    Q = sum(wi*(bi-bbar)**2 for (bi, _), wi in zip(est, w))
    return max(0.0, (Q-(k-1))/Q*100) if Q > 0 else 0.0

def fixed_effect(pairs):
    w = [1/(s*s) for _, s in pairs]
    sw = sum(w)
    b = sum(wi*bi for (bi, _), wi in zip(pairs, w)) / sw
    se = math.sqrt(1/sw)
    return b, se

def eff_bcise(ct, mk, stratum):
    e = eff_rows.get((ct, mk, stratum))
    if e is None: return None
    b = f(e['beta'])
    if b is None: return None
    lo, hi = f(e.get('ci_lo')), f(e.get('ci_hi'))
    if lo is None or hi is None:
        se = f(e['se'])
        if se is None: return None
        lo, hi = b - 1.96*se, b + 1.96*se
    return b, lo, hi

# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(13, 15))
fig.patch.set_facecolor('white')
outer = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1.05, 1.55, 1.35],
                          hspace=0.44, left=0.075, right=0.965,
                          top=0.955, bottom=0.055)

# ── Panel A : global heterogeneity comparison (scatter) ───────────────────────
axA = fig.add_subplot(outer[0])
pos = axA.get_position()
axA.set_position([pos.x0, pos.y0, pos.width*0.5, pos.height])

axA.plot([0, 100], [0, 100], ls='--', color='#BBBBBB', lw=1.0, zorder=1)
rep_set = set(REPS)
for lead in leads:
    ct, mk = lead['cell_type'], lead['marker']
    xc = pooled_i2.get((ct, mk))
    ya = across_ancestry_i2(ct, mk)
    if xc is None or ya is None:
        continue
    if (ct, mk) in rep_set:
        continue
    axA.scatter(xc, ya, s=42, color='#9AA4B0', edgecolors='white',
                linewidths=0.6, alpha=0.9, zorder=3)
rep_colors = ['#1B9E77', '#7570B3', '#B22222']
for (ct, mk), rc in zip(REPS, rep_colors):
    xc = pooled_i2.get((ct, mk)); ya = across_ancestry_i2(ct, mk)
    if xc is None or ya is None: continue
    axA.scatter(xc, ya, s=115, color=rc, edgecolors='white', linewidths=1.1,
                zorder=5)
    axA.annotate(f"{short_ct(ct)}\n{locus_label(mk)}", (xc, ya),
                 textcoords='offset points', xytext=(8, 6), fontsize=8.2,
                 fontweight='bold', color=rc)
axA.set_xlim(-3, 103); axA.set_ylim(-3, 103)
axA.set_xlabel("Across-cohort heterogeneity (I², %)", fontsize=9.5)
axA.set_ylabel("Across-ancestry\nheterogeneity (I², %)", fontsize=9.5)
axA.tick_params(labelsize=8.5)
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.text(30, 92, "line = equal cohort &\nancestry heterogeneity",
         fontsize=7.6, color='#888888', style='italic', ha='left', va='top')
axA.text(-0.14, 1.03, 'A', transform=axA.transAxes,
         fontsize=15, fontweight='bold', va='bottom', ha='left')
axA.set_title("Global comparison of heterogeneity across the 19 lead loci",
              fontsize=10, pad=8, loc='left')

# right-hand key / take-home for Panel A
axK = fig.add_axes([pos.x0 + pos.width*0.60, pos.y0, pos.width*0.40, pos.height])
axK.axis('off')
axK.text(0.0, 1.0, "Representative loci (shown in B–C)", fontsize=9.5,
         fontweight='bold', va='top', ha='left', transform=axK.transAxes)
ky = 0.88
for (ct, mk), rc in zip(REPS, rep_colors):
    xc = pooled_i2.get((ct, mk)); ya = across_ancestry_i2(ct, mk)
    axK.scatter(0.03, ky, s=90, color=rc, edgecolors='white', linewidths=1.0,
                transform=axK.transAxes, clip_on=False)
    axK.text(0.10, ky, f"{short_ct(ct)}  {locus_label(mk)}", fontsize=8.6,
             va='center', ha='left', color=rc, fontweight='bold',
             transform=axK.transAxes)
    axK.text(0.10, ky - 0.075,
             f"I²(cohort)={xc:.0f}%   I²(ancestry)={ya:.0f}%",
             fontsize=7.8, va='center', ha='left', color='#555555',
             transform=axK.transAxes)
    ky -= 0.20
axK.text(0.0, 0.20,
         "Heterogeneity is driven mainly across cohorts,\n"
         "not across ancestries; effect directions stay\n"
         "concordant even where I² is high (see B–C).",
         fontsize=8.0, va='top', ha='left', color='#333333',
         style='italic', transform=axK.transAxes)

# ── Panels B & C : per-locus detail for the 3 representatives ──────────────────
gsB = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1], wspace=0.42)
gsC = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[2], wspace=0.42)

for col, (ct, mk) in enumerate(REPS):
    entries = coh.get((ct, mk), [])
    entries.sort(key=lambda e: (e[1] != 'EUR', COHORT_RANK.get(e[0], 99)))
    eur = [e for e in entries if e[1] == 'EUR']
    afr = [e for e in entries if e[1] == 'AFR']

    # ---- Panel B : grouped ancestry + cohort forest -------------------------
    axB = fig.add_subplot(gsB[col])
    rowsB = []  # (label, beta, lo, hi, kind)
    m = eff_bcise(ct, mk, 'Pooled')
    if m: rowsB.append(('Meta', *m, 'meta'))
    e = eff_bcise(ct, mk, 'EUR')
    if e: rowsB.append(('EUR meta', *e, 'eur_sub'))
    for c, _, b, s in eur:
        rowsB.append((COHORT_LABEL.get(c, c), b, b-1.96*s, b+1.96*s, 'eur'))
    a = eff_bcise(ct, mk, 'AFR')
    if a: rowsB.append(('AFR meta', *a, 'afr_sub'))
    for c, _, b, s in afr:
        rowsB.append((COHORT_LABEL.get(c, c), b, b-1.96*s, b+1.96*s, 'afr'))

    nrows = len(rowsB)
    yb = np.arange(nrows)[::-1]
    axB.axvline(0, color='grey', lw=0.8, ls='--', zorder=1)
    style = {
        'meta':    dict(marker='D', color=META_CLR, ms=8,  ew=1.5),
        'eur_sub': dict(marker='D', color=EUR_CLR,  ms=7,  ew=1.4),
        'afr_sub': dict(marker='D', color=AFR_CLR,  ms=7,  ew=1.4),
        'eur':     dict(marker='o', color=EUR_LT,   ms=5,  ew=1.0),
        'afr':     dict(marker='o', color=AFR_LT,   ms=5,  ew=1.0),
    }
    for yy, (lab, b, lo, hi, kind) in zip(yb, rowsB):
        st = style[kind]
        axB.errorbar(b, yy, xerr=[[b-lo], [hi-b]], fmt=st['marker'],
                     color=st['color'], markersize=st['ms'], markeredgecolor='white',
                     markeredgewidth=0.5, elinewidth=st['ew'], capsize=0, zorder=4)
    axB.set_yticks(yb)
    lab_colors = {'meta':META_CLR,'eur_sub':EUR_CLR,'afr_sub':AFR_CLR,
                  'eur':'#333333','afr':'#333333'}
    axB.set_yticklabels([r[0] for r in rowsB], fontsize=7.4)
    for tick, (_, _, _, _, kind) in zip(axB.get_yticklabels(), rowsB):
        tick.set_color(lab_colors[kind])
        if kind in ('meta', 'eur_sub', 'afr_sub'):
            tick.set_fontweight('bold')
    axB.set_ylim(-0.6, nrows - 0.4)
    axB.tick_params(axis='x', labelsize=8)
    axB.set_xlabel("Effect size (β)", fontsize=8.5)
    axB.spines['top'].set_visible(False)
    axB.spines['right'].set_visible(False)
    i2c = pooled_i2.get((ct, mk))
    axB.set_title(f"{short_ct(ct)}  {locus_label(mk)}\nI²(cohort)={i2c:.0f}%",
                  fontsize=8.8, pad=5)
    if col == 0:
        axB.text(-0.55, 1.06, 'B', transform=axB.transAxes,
                 fontsize=15, fontweight='bold', va='bottom', ha='left')

    # ---- Panel C : leave-one-cohort-out -------------------------------------
    axC = fig.add_subplot(gsC[col])
    pairs = [(b, s) for _, _, b, s in entries]
    labels_c = [c for c, _, _, _ in entries]
    b_all, se_all = fixed_effect(pairs)
    loo = []
    for i in range(len(pairs)):
        sub = pairs[:i] + pairs[i+1:]
        bb, ss = fixed_effect(sub)
        loo.append((labels_c[i], bb, ss))

    rowsC = [('All cohorts', b_all, se_all, 'all')] + \
            [(COHORT_LABEL.get(c, c), bb, ss, 'loo') for c, bb, ss in loo]
    nC = len(rowsC)
    yc = np.arange(nC)[::-1]
    # shaded band = full-data 95% CI
    axC.axvspan(b_all-1.96*se_all, b_all+1.96*se_all, color='#E8EDF3', zorder=0)
    axC.axvline(b_all, color=META_CLR, lw=1.0, ls=':', zorder=1)
    axC.axvline(0, color='grey', lw=0.8, ls='--', zorder=1)
    for yy, (lab, b, s, kind) in zip(yc, rowsC):
        if kind == 'all':
            axC.errorbar(b, yy, xerr=1.96*s, fmt='D', color=META_CLR, markersize=7,
                         markeredgecolor='white', markeredgewidth=0.6,
                         elinewidth=1.5, capsize=0, zorder=4)
        else:
            axC.errorbar(b, yy, xerr=1.96*s, fmt='o', color=LOO_CLR, markersize=4.2,
                         markeredgecolor='white', markeredgewidth=0.4,
                         elinewidth=1.0, capsize=0, zorder=3)
    axC.set_yticks(yc)
    axC.set_yticklabels([r[0] for r in rowsC], fontsize=7.4)
    axC.get_yticklabels()[0].set_fontweight('bold')
    axC.set_ylim(-0.6, nC - 0.4)
    axC.tick_params(axis='x', labelsize=8)
    axC.set_xlabel("Pooled β (leave-one-out)", fontsize=8.5)
    axC.spines['top'].set_visible(False)
    axC.spines['right'].set_visible(False)
    axC.set_title(f"{short_ct(ct)}  {locus_label(mk)}", fontsize=8.8, pad=5)
    if col == 0:
        axC.text(-0.55, 1.06, 'C', transform=axC.transAxes,
                 fontsize=15, fontweight='bold', va='bottom', ha='left')

# shared legends
figB_handles = [
    Line2D([0],[0], marker='D', color=META_CLR, ls='none', markersize=8, label='Overall meta'),
    Line2D([0],[0], marker='D', color=EUR_CLR, ls='none', markersize=7, label='EUR / AFR subtotal'),
    Line2D([0],[0], marker='o', color=EUR_LT, ls='none', markersize=6, label='EUR cohort'),
    Line2D([0],[0], marker='o', color=AFR_LT, ls='none', markersize=6, label='AFR cohort'),
]
fig.legend(handles=figB_handles, loc='center', bbox_to_anchor=(0.5, 0.345),
           ncol=4, fontsize=8.3, frameon=False, columnspacing=1.6,
           handletextpad=0.4)
figC_handles = [
    Line2D([0],[0], marker='D', color=META_CLR, ls='none', markersize=7, label='All cohorts (±95% CI)'),
    Line2D([0],[0], marker='o', color=LOO_CLR, ls='none', markersize=5, label='Leave-one-cohort-out'),
    plt.Rectangle((0,0),1,1, color='#E8EDF3', label='Full-data 95% CI'),
]
fig.legend(handles=figC_handles, loc='center', bbox_to_anchor=(0.5, 0.020),
           ncol=3, fontsize=8.3, frameon=False, columnspacing=1.6,
           handletextpad=0.4)

for ext in ['png', 'svg', 'pdf']:
    out = OUTD / f"figure5_heterogeneity_main.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")
plt.close(fig)
