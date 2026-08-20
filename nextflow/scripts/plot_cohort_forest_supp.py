#!/usr/bin/env python3
"""
Supplementary: per-cohort forest plots for all 19 lead SNPs.

Layout: one tall figure; each lead SNP occupies a horizontal block.
Rows within a block = individual cohorts (coloured by ancestry: blue=EUR, orange=AFR, green=AMR),
  followed by ancestry meta-analysis diamonds (EUR/AFR/Pooled).
A thin grey band separates blocks.

Outputs:
  results/meta_sensitivity/ancestry_lead_effects/figures/supp_cohort_forest.{png,svg,pdf}
"""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

ROOT = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
COH  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_cohort_effects.tsv"
EFF  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

# ── cohort display order (EUR first, then AFR, then AMR) ─────────────────────
EUR_ORDER = ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN',
             'CMC_PITT','GTEx_v10','NABEC','GVEX']
AFR_ORDER = ['NIMH_HBCC_1M','NIMH_HBCC_h650','NIMH_HBCC_Omni5M','AMP_AD_Rush']
AMR_ORDER = ['AMP_AD_Mayo']
COHORT_ORDER = EUR_ORDER + AFR_ORDER + AMR_ORDER
COHORT_LABEL = {
    'ROSMAP':'ROSMAP','ROSMAP_array':'ROSMAP-array','Mayo':'Mayo (EUR)',
    'MSBB':'MSBB','CMC_MSSM':'CMC-MSSM','CMC_PENN':'CMC-PENN',
    'CMC_PITT':'CMC-PITT','GTEx_v10':'GTEx v10','NABEC':'NABEC','GVEX':'GVEX',
    'NIMH_HBCC_1M':'HBCC-1M','NIMH_HBCC_h650':'HBCC-h650',
    'NIMH_HBCC_Omni5M':'HBCC-5M','AMP_AD_Rush':'Rush (AFR)',
    'AMP_AD_Mayo':'Mayo (LAT/AMR)',
}
ANCESTRY_CLR = {'EUR':'#2166AC','AFR':'#D95F02','AMR':'#1B9E77','?':'grey'}
META_CLR     = {'EUR':'#2166AC','AFR':'#D95F02','Pooled':'#333333'}

def f(v):
    try: return float(v)
    except: return None

def short_ct(ct):
    return (ct.replace('Oligodendrocyte','Oligo').replace('Endothelial','Endoth.'))

def locus_label(mk):
    p = mk.split(':')
    return f"chr{p[0].replace('chr','')}: {int(p[1]):,}  {p[2]}>{p[3]}"

# ── load data ─────────────────────────────────────────────────────────────────
coh_rows = list(csv.DictReader(open(COH), delimiter='\t'))
eff_rows = {(r['cell_type'],r['marker'],r['stratum']): r
            for r in csv.DictReader(open(EFF), delimiter='\t')}
het_rows = {(r['cell_type'],r['marker'],r['stratum']): r
            for r in csv.DictReader(open(HET), delimiter='\t')}

# per-lead cohort data
by_lead = {}
for r in coh_rows:
    k = (r['cell_type'],r['marker'])
    by_lead.setdefault(k,[]).append(r)

# ordered leads
seen, leads = set(), []
for r in csv.DictReader(open(HET), delimiter='\t'):
    k = (r['cell_type'],r['marker'])
    if k not in seen and r['stratum']=='Pooled':
        seen.add(k)
        leads.append({'cell_type':r['cell_type'],'marker':r['marker'],
                      'i2':f(r['i2']),'het_p':f(r['het_p']),'beta':f(r['beta'])})

# ── compute total row count for figure height ─────────────────────────────────
SUMMARY_STRATA = ['EUR','AFR','Pooled']
ROW_HT_IN = 0.30   # inches per row
SEP_HT_IN = 0.40   # separator between blocks
HEADER_HT = 0.40   # lead-SNP header

blocks = []
for lead in leads:
    k = (lead['cell_type'], lead['marker'])
    cohort_data = {r['cohort']: r for r in by_lead.get(k,[]) if r['beta']}
    # display in fixed order
    cohort_list = [c for c in COHORT_ORDER if c in cohort_data]
    # summary rows that have data
    sum_list = [s for s in SUMMARY_STRATA
                if f(eff_rows.get((k[0],k[1],s),{}).get('beta','')) is not None]
    n_rows = len(cohort_list) + len(sum_list)
    blocks.append({'lead':lead,'cohort_list':cohort_list,'sum_list':sum_list,
                   'cohort_data':cohort_data,'n_rows':n_rows})

total_h = sum(b['n_rows'] * ROW_HT_IN + HEADER_HT + SEP_HT_IN for b in blocks) + 1.5
total_h = max(total_h, 12)

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, total_h))
fig.patch.set_facecolor('white')

# We'll draw everything using data coordinates y=0..total_rows (top to bottom)
# Accumulate y position
LABEL_W = 2.5    # inches for left label — will use axis fraction
ax.set_xlim(-0.6, 0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', left=False, labelleft=False)
ax.set_xlabel("Beta (95% CI)\n(harmonized to pooled effect allele)", fontsize=9)
ax.tick_params(axis='x', labelsize=8)
ax.axvline(0, color='grey', lw=0.6, ls='--', zorder=1)

# Accumulate rows top→bottom
y_cur = 0          # will increment downward (we flip y at end)
y_ticks = []
y_labels = []
y_anc_labels = []  # (y, label, color) for ancestry labels on left

for bi, block in enumerate(blocks):
    lead = block['lead']
    ct   = lead['cell_type']
    mk   = lead['marker']

    # ── block header (lead SNP title) ──
    y_header = y_cur + 0.5
    ax.text(-0.59, y_header,
            f"{short_ct(ct)}  |  {locus_label(mk)}  |  I²={lead['i2']:.0f}%",
            fontsize=9, fontweight='bold', va='center', ha='left',
            color='#333333')
    y_cur += 1.0

    # ── cohort rows ──
    cohort_rows_y = []
    for cohort in block['cohort_list']:
        r   = block['cohort_data'][cohort]
        bv  = f(r['beta']); se = f(r['se'])
        anc = r.get('ancestry','?')
        col = ANCESTRY_CLR.get(anc, 'grey')
        if bv is not None and se is not None:
            ax.errorbar(bv, y_cur, xerr=1.96*se,
                        fmt='o', color=col, markersize=3.2,
                        linewidth=0.7, capsize=1.5, capthick=0.6, zorder=3)
        y_ticks.append(y_cur)
        y_labels.append(COHORT_LABEL.get(cohort, cohort))
        cohort_rows_y.append((y_cur, r.get('ancestry','?')))
        y_cur += 1.0

    # ancestry bracket labels (EUR / AFR / AMR) on far left
    if cohort_rows_y:
        anc_groups = {}
        for yy, anc in cohort_rows_y:
            anc_groups.setdefault(anc, []).append(yy)
        for anc, ys in anc_groups.items():
            mid = (min(ys) + max(ys)) / 2
            y_anc_labels.append((mid, anc, ANCESTRY_CLR.get(anc,'grey')))

    # separator line
    ax.axhline(y_cur - 0.5, color='#DDDDDD', lw=0.8, zorder=0)

    # ── summary diamonds ──
    for s in block['sum_list']:
        e   = eff_rows.get((ct, mk, s), {})
        bv  = f(e.get('beta','')); se = f(e.get('se',''))
        pv  = f(e.get('p',''))
        col = META_CLR.get(s,'#333333')
        if bv is not None and se is not None:
            ax.errorbar(bv, y_cur, xerr=1.96*se,
                        fmt='D', color=col, markersize=5.5,
                        linewidth=1.1, capsize=2, capthick=0.8, zorder=4)
            if pv is not None:
                ps = f"{pv:.1e}" if pv < 0.001 else f"{pv:.3f}"
                ax.text(0.61, y_cur, f"P={ps}", va='center', ha='left',
                        fontsize=7, color=col,
                        transform=ax.get_yaxis_transform())
        slbl = {'EUR':'EUR meta','AFR':'AFR meta','Pooled':'Pooled (15 cohorts)'}[s]
        y_ticks.append(y_cur)
        y_labels.append(slbl)
        y_cur += 1.0

    # gap between blocks
    y_cur += 0.6

total_rows = y_cur
ax.set_ylim(total_rows, -0.5)   # flip: 0 at top, total at bottom
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels, fontsize=8.5)

# ancestry labels far left (requires axes transform)
for yy, lbl, col in y_anc_labels:
    ax.text(-0.63, yy, lbl, va='center', ha='right',
            fontsize=8.5, fontweight='bold', color=col,
            transform=ax.get_yaxis_transform())

# ── legend ───────────────────────────────────────────────────────────────────
leg = [
    Line2D([0],[0], marker='o', color=ANCESTRY_CLR['EUR'], ls='none',
           markersize=5, label='EUR cohort'),
    Line2D([0],[0], marker='o', color=ANCESTRY_CLR['AFR'], ls='none',
           markersize=5, label='AFR cohort'),
    Line2D([0],[0], marker='o', color=ANCESTRY_CLR['AMR'], ls='none',
           markersize=5, label='LAT/AMR cohort'),
    Line2D([0],[0], marker='D', color=META_CLR['EUR'],    ls='none',
           markersize=6, label='EUR meta-analysis'),
    Line2D([0],[0], marker='D', color=META_CLR['AFR'],    ls='none',
           markersize=6, label='AFR meta-analysis'),
    Line2D([0],[0], marker='D', color=META_CLR['Pooled'], ls='none',
           markersize=6, label='Pooled (15-cohort) meta'),
]
ax.legend(handles=leg, fontsize=8, loc='lower right',
          frameon=True, framealpha=0.9, edgecolor='none')

fig.subplots_adjust(left=0.25, right=0.88, top=0.98, bottom=0.05)

for ext in ['png','svg','pdf']:
    out = OUTD / f"supp_cohort_forest.{ext}"
    fig.savefig(out, dpi=300 if ext=='png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")

plt.close(fig)
