#!/usr/bin/env python3
"""
Cross-COHORT heterogeneity figure (split from the old combined Figure 5).

A 4×5 grid of per-cohort forest plots — one per lead variant (all 19 leads).
Each cohort's harmonized estimate (β ± 95% CI) is a neutral-coloured circle
(cohorts can be ancestry-mixed); the pooled 15-cohort meta is a black diamond.
Titles give cell type and dbSNP rsID (chr:pos fallback for the indel with no rsID).

Outputs: results/.../figures/supp_cross_cohort_forest.{png,svg,pdf}  (supplementary)
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
EFF  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
COH  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_cohort_effects.tsv"
RSID = ROOT / "results/meta_sensitivity/ancestry_lead_effects/lead_rsids.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

def f(v):
    try: return float(v)
    except: return None

def short_ct(ct):
    return ct.replace('Oligodendrocyte', 'Oligo')

RSID_MAP = {}
if RSID.exists():
    RSID_MAP = {r['marker']: r['rsid']
                for r in csv.DictReader(open(RSID), delimiter='\t')}
def rsid_of(mk):
    rs = RSID_MAP.get(mk, '.')
    return rs if rs and rs != '.' else mk.split(':')[0].replace('chr', '') + ':' + mk.split(':')[1]

eff_rows = {(r['cell_type'], r['marker'], r['stratum']): r
            for r in csv.DictReader(open(EFF), delimiter='\t')}
het_all  = list(csv.DictReader(open(HET), delimiter='\t'))

seen, leads = set(), []
for r in het_all:
    k = (r['cell_type'], r['marker'])
    if k not in seen and r['stratum'] == 'Pooled':
        seen.add(k)
        # i2 / het_p here are the across-COHORT heterogeneity (METAL pooled 15-cohort)
        leads.append({'cell_type': r['cell_type'], 'marker': r['marker'],
                      'i2': f(r['i2']), 'het_p': f(r['het_p'])})
# order by chromosome then position (keeps co-localised loci adjacent)
leads.sort(key=lambda l: (int(l['marker'].split(':')[0].replace('chr','')),
                          int(l['marker'].split(':')[1])))
N = len(leads)

coh_rows = list(csv.DictReader(open(COH), delimiter='\t'))
coh_full = {}
for r in coh_rows:
    if r['beta'] and r['status'] in ('ok', 'flipped'):
        coh_full.setdefault((r['cell_type'], r['marker']), []).append(
            {'cohort': r['cohort'], 'beta': f(r['beta']), 'se': f(r['se'])})

EUR_ORDER = ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN',
             'CMC_PITT','GTEx_v10','NABEC','GVEX']
AFR_ORDER = ['NIMH_HBCC_1M','NIMH_HBCC_h650','NIMH_HBCC_Omni5M','AMP_AD_Rush']
AMR_ORDER = ['AMP_AD_Mayo']
COHORT_RANK = {c: i for i, c in enumerate(EUR_ORDER + AFR_ORDER + AMR_ORDER)}
COHORT_LABEL = {
    'ROSMAP':'ROSMAP','ROSMAP_array':'ROSMAP-array','Mayo':'Mayo','MSBB':'MSBB',
    'CMC_MSSM':'CMC-MSSM','CMC_PENN':'CMC-PENN','CMC_PITT':'CMC-PITT',
    'GTEx_v10':'GTEx v10','NABEC':'NABEC','GVEX':'GVEX',
    'NIMH_HBCC_1M':'HBCC-1M','NIMH_HBCC_h650':'HBCC-h650','NIMH_HBCC_Omni5M':'HBCC-5M',
    'AMP_AD_Rush':'Rush','AMP_AD_Mayo':'Mayo-LAT',
}
META_CLR = '#333333'
COH_CLR  = '#5B7FA6'

# global symmetric x-range for comparability
allb = []
for k, recs in coh_full.items():
    for rec in recs:
        allb.append(abs(rec['beta']) + 1.96 * rec['se'])
XR = (np.percentile(allb, 98) if allb else 0.6) * 1.05

NCOL, NROW = 4, 5
fig = plt.figure(figsize=(17, 21))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(NROW, NCOL, figure=fig, hspace=0.42, wspace=0.55,
                       left=0.07, right=0.985, top=0.965, bottom=0.045)

for idx, lead in enumerate(leads):
    ax = fig.add_subplot(gs[idx])
    ct, mk = lead['cell_type'], lead['marker']
    recs = sorted(coh_full.get((ct, mk), []), key=lambda r: COHORT_RANK.get(r['cohort'], 99))

    rows = []   # (label, beta, se, colour, is_meta)
    e = eff_rows.get((ct, mk, 'Pooled'))
    if e is not None and f(e['beta']) is not None and f(e['se']) is not None:
        rows.append(('Meta', f(e['beta']), f(e['se']), META_CLR, True))
    for rec in reversed(recs):
        rows.append((COHORT_LABEL.get(rec['cohort'], rec['cohort']),
                     rec['beta'], rec['se'], COH_CLR, False))

    yy = np.arange(len(rows))
    ax.axvline(0, color='grey', lw=0.6, ls='--', zorder=1)
    n_meta = sum(1 for r in rows if r[4])
    if n_meta:
        ax.axhline(n_meta - 0.5, color='#CCCCCC', lw=0.8, zorder=1)
    for yi, (lbl, b, se, col, is_meta) in zip(yy, rows):
        ax.errorbar(b, yi, xerr=1.96 * se,
                    fmt='D' if is_meta else 'o', color=col,
                    markersize=6.0 if is_meta else 5.0, markeredgewidth=0,
                    linewidth=1.5 if is_meta else 1.0, capsize=0, zorder=3)

    ax.set_xlim(-XR, XR)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.set_xticks([t for t in ax.get_xticks() if -XR <= t <= XR])
    ax.tick_params(axis='x', labelsize=8.5)
    ax.set_xlabel("β (95% CI)", fontsize=9.5)
    ax.set_title(f"{short_ct(ct)}  ({rsid_of(mk)})",
                 fontsize=10.5, fontweight='bold', pad=4, loc='left')
    # across-cohort heterogeneity annotation (top-right of each panel)
    i2v, hp = lead.get('i2'), lead.get('het_p')
    if i2v is not None:
        htxt = f"I²={i2v:.0f}%"
        if hp is not None:
            htxt += f"\n$P_{{het}}$=" + (f"{hp:.0e}" if hp < 1e-3 else f"{hp:.2f}")
        ax.text(0.97, 0.97, htxt, transform=ax.transAxes, ha='right', va='top',
                fontsize=7.5, color='#555',
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#DDDDDD', lw=0.5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# hide unused grid cells
for idx in range(N, NCOL * NROW):
    ax = fig.add_subplot(gs[idx]); ax.axis('off')

fig.legend(handles=[
    Line2D([0],[0], marker='o', color=COH_CLR, ls='none', markersize=7, label='cohort'),
    Line2D([0],[0], marker='D', color=META_CLR, ls='none', markersize=7, label='meta-analysis'),
], fontsize=10, loc='lower right', frameon=True, framealpha=0.9, edgecolor='none',
   bbox_to_anchor=(0.985, 0.012), ncol=2)

for ext in ['png', 'svg', 'pdf']:
    out = OUTD / f"supp_cross_cohort_forest.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")
plt.close(fig)
