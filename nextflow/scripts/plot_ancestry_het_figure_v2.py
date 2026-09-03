#!/usr/bin/env python3
"""
Three-panel ancestry heterogeneity figure (Nature Genetics PD-GWAS Fig.3 style).

Panel A  (top)    : I² bar chart per lead SNP – total cohort heterogeneity,
                    colour by level (<30%, 30-75%, ≥75%), significance stars.
Panel B  (middle) : Per-cohort signed-effect direction matrix.
                    Rows = 15 cohorts grouped by ancestry (EUR / AFR).
                    Columns = 19 lead SNPs.
                    Colour: blue = effect same as pooled direction (+),
                            red  = opposite direction (-),
                            grey = variant absent (?).
                    Circle size proportional to abs(β) of the pooled estimate
                    (all cells same column same size – shows locus magnitude).
Panel C  (right)  : Ancestry-group forest plots (EUR / AFR / Pooled).

Outputs:
  results/meta_sensitivity/ancestry_lead_effects/figures/ancestry_het_figure_v2.{svg,png,pdf}
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

# ── paths ─────────────────────────────────────────────────────────────────────
ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
EFF  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
DIR  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_directions.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

# ── cohort order & ancestry ───────────────────────────────────────────────────
COHORTS = ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN','CMC_PITT',
           'GTEx_v10','NABEC','GVEX','NIMH_HBCC_1M','NIMH_HBCC_h650',
           'NIMH_HBCC_Omni5M','AMP_AD_Rush','AMP_AD_Mayo']
# Reorder: EUR first (10), then AFR (5)
EUR_COHORTS = ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN',
               'CMC_PITT','GTEx_v10','NABEC','GVEX']
AFR_COHORTS = ['NIMH_HBCC_1M','NIMH_HBCC_h650','NIMH_HBCC_Omni5M',
               'AMP_AD_Rush','AMP_AD_Mayo']
# METAL Direction string index for each cohort (from nextflow.config cohorts order)
METAL_IDX = {c: i for i, c in enumerate(
    ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN','CMC_PITT',
     'GTEx_v10','NABEC','NIMH_HBCC_1M','NIMH_HBCC_h650','NIMH_HBCC_Omni5M',
     'GVEX','AMP_AD_Rush','AMP_AD_Mayo'])}
DISPLAY_ORDER = EUR_COHORTS + AFR_COHORTS

COHORT_LABEL = {
    'ROSMAP':'ROSMAP','ROSMAP_array':'ROSMAP-array','Mayo':'Mayo',
    'MSBB':'MSBB','CMC_MSSM':'CMC-MSSM','CMC_PENN':'CMC-PENN',
    'CMC_PITT':'CMC-PITT','GTEx_v10':'GTEx v10','NABEC':'NABEC','GVEX':'GVEX',
    'NIMH_HBCC_1M':'HBCC-1M','NIMH_HBCC_h650':'HBCC-h650',
    'NIMH_HBCC_Omni5M':'HBCC-5M','AMP_AD_Rush':'Rush-AFR','AMP_AD_Mayo':'Mayo-AFR',
}

# ── colours ───────────────────────────────────────────────────────────────────
CLR_POS    = "#2166AC"   # blue  – same direction as pooled
CLR_NEG    = "#B2182B"   # red   – opposite direction
CLR_MISS   = "#D9D9D9"   # grey  – variant absent
CLR_EUR_BG = "#EEF4FB"
CLR_AFR_BG = "#FEF3EC"
HET_LOW    = "#BDBDBD"
HET_MED    = "#F4A460"
HET_HIGH   = "#B22222"

def i2_colour(v):
    if v is None: return HET_LOW
    if v >= 75:   return HET_HIGH
    if v >= 30:   return HET_MED
    return HET_LOW

def f(v):
    try: return float(v)
    except: return None

# ── load data ─────────────────────────────────────────────────────────────────
dir_rows  = {(r['cell_type'], r['marker']): r
             for r in csv.DictReader(open(DIR), delimiter='\t')}
het_rows  = {(r['cell_type'], r['marker'], r['stratum']): r
             for r in csv.DictReader(open(HET), delimiter='\t')}
eff_rows  = {(r['cell_type'], r['marker'], r['stratum']): r
             for r in csv.DictReader(open(EFF), delimiter='\t')}

# Ordered lead list
seen, leads = set(), []
for r in csv.DictReader(open(HET), delimiter='\t'):
    k = (r['cell_type'], r['marker'])
    if k not in seen and r['stratum'] == 'Pooled':
        seen.add(k)
        leads.append({'cell_type': r['cell_type'], 'marker': r['marker'],
                      'i2': f(r['i2']), 'het_p': f(r['het_p']),
                      'beta': f(r['beta']), 'p': f(r['p'])})

N  = len(leads)    # 19 columns
NC = len(DISPLAY_ORDER)  # 15 rows

def short_ct(ct):
    return ct.replace('Oligodendrocyte','Oligo').replace('Endothelial','Endoth.')

col_labels = [f"{short_ct(r['cell_type'])}\n{r['marker'].split(':')[0].replace('chr','')}:"
              f"{int(r['marker'].split(':')[1]):,}" for r in leads]
row_labels  = [COHORT_LABEL[c] for c in DISPLAY_ORDER]

# ── build direction matrix ────────────────────────────────────────────────────
# dir_mat[cohort_row, lead_col] ∈ {+1, -1, 0}   (0 = missing)
dir_mat  = np.zeros((NC, N), dtype=float)
for j, lead in enumerate(leads):
    dr = dir_rows.get((lead['cell_type'], lead['marker']))
    if dr is None: continue
    direction_str = dr['direction']
    pooled_sign   = 1 if float(dr['beta']) > 0 else -1
    for i, cohort in enumerate(DISPLAY_ORDER):
        idx = METAL_IDX.get(cohort)
        if idx is None or idx >= len(direction_str): continue
        ch = direction_str[idx]
        if ch == '+': dir_mat[i, j] = +1  # same as allele1
        elif ch == '-': dir_mat[i, j] = -1
        else: dir_mat[i, j] = 0

# Flip sign so +1 = same direction as pooled beta, -1 = opposite
for j, lead in enumerate(leads):
    if (lead['beta'] or 0) < 0:
        dir_mat[:, j] = -dir_mat[:, j]

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(20, 14))
fig.patch.set_facecolor('white')

# Left portion: panels A (I² bars) + B (direction matrix)
# Right portion: panel C (forest)
# Use explicit axes positions to guarantee no overlap
# [left, bottom, width, height] in figure fraction

# Panel A: top-left, I² bars
ax_a = fig.add_axes([0.09, 0.72, 0.55, 0.24])
# Panel B: bottom-left, direction matrix (shares x-extent with A)
ax_b = fig.add_axes([0.09, 0.06, 0.55, 0.62])
# Panel C: right column, forest
ax_c = fig.add_axes([0.68, 0.06, 0.26, 0.90])

x = np.arange(N)

# ── Panel A – I² bars ─────────────────────────────────────────────────────────
i2_vals  = [r['i2'] or 0 for r in leads]
bar_cols = [i2_colour(r['i2']) for r in leads]
ax_a.bar(x, i2_vals, width=0.7, color=bar_cols, linewidth=0, zorder=2)

for yref, ls, lbl in [(30,'--','I²=30%'),(75,'-.','I²=75%')]:
    ax_a.axhline(yref, color='black', lw=0.55, ls=ls, alpha=0.45, zorder=1)
    ax_a.text(N - 0.3, yref + 1.5, lbl, fontsize=8, va='bottom', ha='right',
              color='grey')

for i, r in enumerate(leads):
    if r['het_p'] and r['het_p'] < 0.05:
        stars = '***' if r['het_p'] < 0.001 else '**' if r['het_p'] < 0.01 else '*'
        ax_a.text(i, (r['i2'] or 0) + 2, stars,
                  ha='center', va='bottom', fontsize=8, color='black')

ax_a.set_xlim(-0.6, N - 0.4)
ax_a.set_ylim(0, 115)
ax_a.set_xticks(x); ax_a.set_xticklabels(['']*N)
ax_a.set_ylabel('I² (%)', fontsize=9)
ax_a.tick_params(axis='y', labelsize=8)
ax_a.spines['top'].set_visible(False); ax_a.spines['right'].set_visible(False)
# Panel label placed inside the axes, top-left corner — no overflow into other panels
ax_a.text(0.01, 0.97, 'A', transform=ax_a.transAxes,
          fontsize=14, fontweight='bold', va='top', ha='left')

leg_a = [mpatches.Patch(color=HET_LOW,  label='I² < 30%'),
         mpatches.Patch(color=HET_MED,  label='I² 30–75%'),
         mpatches.Patch(color=HET_HIGH, label='I² ≥ 75%')]
ax_a.legend(handles=leg_a, fontsize=8, loc='upper right',
            frameon=True, framealpha=0.9, edgecolor='none', handlelength=0.9)

# ── Panel B – direction matrix ─────────────────────────────────────────────────
# Background bands for EUR / AFR rows
n_eur = len(EUR_COHORTS); n_afr = len(AFR_COHORTS)
ax_b.axhspan(n_eur - 0.5, NC - 0.5, color=CLR_AFR_BG, zorder=0)
ax_b.axhspan(-0.5, n_eur - 0.5, color=CLR_EUR_BG, zorder=0)

# Column separators
for j in np.arange(N) + 0.5:
    ax_b.axvline(j, color='white', lw=0.8, zorder=2)

CIRCLE_BASE = 280  # pt² base scatter size
for j, lead in enumerate(leads):
    abs_beta = abs(lead['beta'] or 0)
    sz = max(30, CIRCLE_BASE * abs_beta / 0.45)   # normalise to ~max |β| = 0.45
    for i in range(NC):
        v = dir_mat[i, j]
        if v == 0:
            col = CLR_MISS; marker = 'x'; alpha = 0.4
        elif v > 0:
            col = CLR_POS; marker = 'o'; alpha = 0.85
        else:
            col = CLR_NEG; marker = 'o'; alpha = 0.85
        ax_b.scatter(j, i, s=sz if v != 0 else 25,
                     c=col, marker=marker, alpha=alpha,
                     linewidths=0.3, edgecolors='white', zorder=3)

# Ancestry side labels
ax_b.text(-1.3, (n_eur - 1) / 2, 'EUR', rotation=90, va='center', ha='center',
          fontsize=9, fontweight='bold', color='#2166AC')
ax_b.text(-1.3, n_eur + (n_afr - 1) / 2, 'AFR', rotation=90, va='center', ha='center',
          fontsize=9, fontweight='bold', color='#D95F02')
# dividing line EUR/AFR
ax_b.axhline(n_eur - 0.5, color='grey', lw=0.8, ls='--', zorder=4)

ax_b.set_xlim(-0.6, N - 0.4)
ax_b.set_ylim(-0.6, NC - 0.4)
ax_b.set_xticks(x)
ax_b.set_xticklabels(col_labels, fontsize=7.5, rotation=45, ha='right')
ax_b.set_yticks(np.arange(NC))
ax_b.set_yticklabels(row_labels, fontsize=8.5)
ax_b.tick_params(axis='x', length=2)
ax_b.tick_params(axis='y', length=0)
ax_b.spines['top'].set_visible(False); ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_visible(False); ax_b.spines['bottom'].set_visible(False)
# Panel label inside top-left of B axes
ax_b.text(0.01, 0.99, 'B', transform=ax_b.transAxes,
          fontsize=14, fontweight='bold', va='top', ha='left')

leg_b = [Line2D([0],[0], marker='o', color=CLR_POS, ls='none', markersize=7,
                label='Same direction as pooled'),
         Line2D([0],[0], marker='o', color=CLR_NEG, ls='none', markersize=7,
                label='Opposite direction'),
         Line2D([0],[0], marker='x', color=CLR_MISS, ls='none', markersize=7,
                label='Variant absent')]
ax_b.legend(handles=leg_b, fontsize=8, loc='lower right',
            frameon=True, framealpha=0.9, edgecolor='none')

# ── Panel C – ancestry forest ──────────────────────────────────────────────────
STRATA    = ['EUR','AFR','Pooled']
CLR_F     = {'EUR':'#2166AC','AFR':'#D95F02','Pooled':'#333333'}
offsets   = {'EUR': 0.18, 'AFR': 0.0, 'Pooled': -0.18}
y_leads   = np.arange(N)[::-1]

ax_c.axvline(0, color='grey', lw=0.6, ls='--', zorder=1)

for i, (lead, yl) in enumerate(zip(leads, y_leads)):
    ct, mk = lead['cell_type'], lead['marker']
    if i % 2 == 0:
        ax_c.axhspan(yl - 0.5, yl + 0.5, color='#F5F5F5', zorder=0)

    # I² label
    i2_str = f"I²={lead['i2']:.0f}%" if lead['i2'] is not None else 'I²=NA'
    ax_c.text(1.01, yl, i2_str, transform=ax_c.get_yaxis_transform(),
              ha='left', va='center', fontsize=7.5,
              color=i2_colour(lead['i2']))

    for s in STRATA:
        e = eff_rows.get((ct, mk, s))
        if e is None: continue
        beta = f(e['beta']); se = f(e['se']); p = f(e['p'])
        if beta is None or se is None: continue
        yo = yl + offsets[s]
        is_pooled = (s == 'Pooled')
        ax_c.errorbar(beta, yo, xerr=1.96*se,
                      fmt='D' if is_pooled else 'o',
                      color=CLR_F[s],
                      markersize=5.0 if is_pooled else 3.2,
                      linewidth=1.0 if is_pooled else 0.65,
                      capsize=1.5, capthick=0.6, zorder=3)
        if is_pooled and p is not None:
            p_str = f"{p:.1e}" if p < 0.001 else f"{p:.3f}"
            ax_c.text(1.13, yl, f"P={p_str}",
                      transform=ax_c.get_yaxis_transform(),
                      ha='left', va='center', fontsize=7, color='black')

ax_c.set_ylim(-0.8, N - 0.2)
ax_c.set_yticks(y_leads)
ax_c.set_yticklabels([f"{short_ct(r['cell_type'])}" for r in leads], fontsize=8.5)
ax_c.set_xlabel('Beta (95% CI)\n(pooled effect allele)', fontsize=9)
ax_c.tick_params(axis='x', labelsize=8)
ax_c.spines['top'].set_visible(False); ax_c.spines['right'].set_visible(False)
# Panel label inside top-left of C axes
ax_c.text(0.02, 0.99, 'C', transform=ax_c.transAxes,
          fontsize=14, fontweight='bold', va='top', ha='left')

leg_c = [Line2D([0],[0], marker='o', color=CLR_F['EUR'], ls='none', markersize=6,
                label="EUR (N≈2,818)"),
         Line2D([0],[0], marker='o', color=CLR_F['AFR'], ls='none', markersize=6,
                label="AFR (N≈270)"),
         Line2D([0],[0], marker='D', color=CLR_F['Pooled'], ls='none', markersize=7,
                label="Pooled (15 cohorts)")]
ax_c.legend(handles=leg_c, fontsize=8, loc='lower right',
            frameon=True, framealpha=0.9, edgecolor='none')

# ── export ────────────────────────────────────────────────────────────────────
for ext in ['png','svg','pdf']:
    out = OUTD / f"ancestry_het_figure_v2.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")
plt.close(fig)
