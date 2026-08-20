#!/usr/bin/env python3
"""
Cross-cohort direction-concordance matrix for genome-wide-significant lead SNPs.

Primary message: top SNPs show consistent effect direction across the 15 cohorts.

Layout:
  Top strip  : per-lead concordance bar (k/N cohorts agreeing with the pooled sign).
  Matrix     : rows = 15 cohorts (grouped EUR / AFR / LAT), columns = lead SNPs.
               Colour = effect direction relative to the pooled meta-analysis:
                 blue  = same direction (concordant),
                 red   = opposite direction (discordant),
                 grey  = variant absent / allele mismatch.
               Colour saturation encodes |beta| (stronger effect = more saturated).
               A small black dot marks cohorts with nominal significance (P < 0.05).
  Bottom     : column labels = cell type + locus; pooled I² annotated.

Columns are ordered by number of cohorts observed (descending) so the best-powered,
most-consistent loci lead the figure.

Outputs:
  results/meta_sensitivity/ancestry_lead_effects/figures/concordance_matrix.{png,svg,pdf}
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

ROOT = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
COH  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_cohort_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

# ── cohort order & ancestry ───────────────────────────────────────────────────
EUR = ['ROSMAP','ROSMAP_array','Mayo','MSBB','CMC_MSSM','CMC_PENN','CMC_PITT',
       'GTEx_v10','NABEC','GVEX']
AFR = ['NIMH_HBCC_1M','NIMH_HBCC_h650','NIMH_HBCC_Omni5M','AMP_AD_Rush']
AMR = ['AMP_AD_Mayo']
COHORT_ORDER = EUR + AFR + AMR
ANC = {**{c:'EUR' for c in EUR}, **{c:'AFR' for c in AFR}, **{c:'AMR' for c in AMR}}
COHORT_LABEL = {
    'ROSMAP':'ROSMAP','ROSMAP_array':'ROSMAP-array','Mayo':'Mayo',
    'MSBB':'MSBB','CMC_MSSM':'CMC-MSSM','CMC_PENN':'CMC-PENN','CMC_PITT':'CMC-PITT',
    'GTEx_v10':'GTEx v10','NABEC':'NABEC','GVEX':'GVEX',
    'NIMH_HBCC_1M':'HBCC-1M','NIMH_HBCC_h650':'HBCC-h650','NIMH_HBCC_Omni5M':'HBCC-5M',
    'AMP_AD_Rush':'Rush','AMP_AD_Mayo':'Mayo-LAT',
}
ANC_CLR = {'EUR':'#2166AC','AFR':'#D95F02','AMR':'#1B9E77'}

NC = len(COHORT_ORDER)

def f(v):
    try: return float(v)
    except: return None

def short_ct(ct):
    return ct.replace('Oligodendrocyte','Oligo').replace('Endothelial','Endoth.')

def locus(mk):
    p = mk.split(':')
    return f"{p[0].replace('chr','')}:{int(p[1])/1e6:.1f}M"

# ── load ──────────────────────────────────────────────────────────────────────
coh = list(csv.DictReader(open(COH), delimiter='\t'))
het = {(r['cell_type'],r['marker']): r
       for r in csv.DictReader(open(HET), delimiter='\t') if r['stratum']=='Pooled'}

# index cohort effects
cell = {}
for r in coh:
    cell[(r['cell_type'], r['marker'], r['cohort'])] = r

# ── build leads with stats, ordered by n_cohort desc ──────────────────────────
seen, leads = set(), []
for r in coh:
    k = (r['cell_type'], r['marker'])
    if k in seen: continue
    seen.add(k)
    hb = het.get(k)
    if hb is None: continue
    pooled_beta = f(hb['beta'])
    betas = []
    for c in COHORT_ORDER:
        rr = cell.get((k[0],k[1],c))
        if rr and rr['status'] in ('ok','flipped') and rr['beta']:
            betas.append(f(rr['beta']))
    n = len(betas)
    same = sum(1 for b in betas if (b>0)==(pooled_beta>0)) if n else 0
    leads.append({'cell_type':k[0],'marker':k[1],'pooled_beta':pooled_beta,
                  'i2':f(hb['i2']),'n':n,'same':same})

leads.sort(key=lambda d: (-d['n'], -(d['same']/d['n'] if d['n'] else 0)))
N = len(leads)

# max |beta| across cohorts for saturation scaling
all_abs = [abs(f(r['beta'])) for r in coh if r['beta'] and r['status'] in ('ok','flipped')]
BMAX = np.percentile(all_abs, 95) if all_abs else 0.3

# ── colours ───────────────────────────────────────────────────────────────────
C_SAME = np.array([33, 102, 172]) / 255   # blue
C_OPP  = np.array([178, 24, 43])  / 255    # red
C_MISS = "#E0E0E0"

def cell_colour(beta, pooled_beta):
    """Return RGBA. Saturation ~ |beta|."""
    same = (beta > 0) == (pooled_beta > 0)
    base = C_SAME if same else C_OPP
    sat  = 0.35 + 0.65 * min(abs(beta) / BMAX, 1.0)
    # blend toward white for lower saturation
    rgb = base * sat + np.array([1,1,1]) * (1 - sat)
    return (*rgb, 1.0)

# ── figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 8.5))
fig.patch.set_facecolor('white')

# grid: top concordance bar + matrix
gs = fig.add_gridspec(2, 1, height_ratios=[1, 4.4], hspace=0.04,
                      left=0.11, right=0.85, top=0.94, bottom=0.20)
ax_top = fig.add_subplot(gs[0])
ax_mat = fig.add_subplot(gs[1], sharex=ax_top)

x = np.arange(N)

# ── top: concordance bar ──────────────────────────────────────────────────────
conc = [100*d['same']/d['n'] if d['n'] else 0 for d in leads]
bar_c = ['#2166AC' if c>=90 else '#7FB0D6' if c>=75 else '#F4A460' if c>=50 else '#B2182B'
         for c in conc]
ax_top.bar(x, conc, width=0.72, color=bar_c, zorder=2)
for i,d in enumerate(leads):
    ax_top.text(i, conc[i]+3, f"{d['same']}/{d['n']}", ha='center', va='bottom',
                fontsize=6.8, color='#333')
ax_top.axhline(100, color='grey', lw=0.5, ls=':', alpha=0.6)
ax_top.set_ylim(0, 118)
ax_top.set_ylabel("Cohorts\nconcordant (%)", fontsize=8.5)
ax_top.tick_params(axis='y', labelsize=8)
ax_top.tick_params(axis='x', bottom=False, labelbottom=False)
ax_top.spines['top'].set_visible(False)
ax_top.spines['right'].set_visible(False)
ax_top.set_title("Cross-cohort effect-direction concordance at genome-wide significant lead SNPs",
                 fontsize=11, fontweight='bold', pad=8)

# ── matrix ────────────────────────────────────────────────────────────────────
# ancestry background bands
n_eur, n_afr = len(EUR), len(AFR)
ax_mat.axhspan(NC-0.5-len(AMR), NC-0.5, color='#EAF6F0', zorder=0)          # AMR
ax_mat.axhspan(n_eur-0.5, NC-0.5-len(AMR), color='#FDF0E6', zorder=0)       # AFR
ax_mat.axhspan(-0.5, n_eur-0.5, color='#EEF4FB', zorder=0)                  # EUR

for j, d in enumerate(leads):
    pb = d['pooled_beta']
    for i, c in enumerate(COHORT_ORDER):
        yy = NC - 1 - i   # cohort 0 at top
        rr = cell.get((d['cell_type'], d['marker'], c))
        if (rr is None) or (not rr['beta']) or (rr['status'] not in ('ok','flipped')):
            fc = C_MISS
            ax_mat.add_patch(mpatches.Rectangle((j-0.42, yy-0.42), 0.84, 0.84,
                             facecolor=fc, edgecolor='white', linewidth=0.5, zorder=2))
            continue
        bv = f(rr['beta']); pv = f(rr['p'])
        fc = cell_colour(bv, pb)
        ax_mat.add_patch(mpatches.Rectangle((j-0.42, yy-0.42), 0.84, 0.84,
                         facecolor=fc, edgecolor='white', linewidth=0.5, zorder=2))
        if pv is not None and pv < 0.05:
            ax_mat.scatter(j, yy, s=7, c='black', marker='.', zorder=3)

# ancestry divider lines
ax_mat.axhline(NC-0.5-n_eur, color='grey', lw=0.8, ls='--', zorder=4)
ax_mat.axhline(NC-0.5-n_eur-n_afr, color='grey', lw=0.8, ls='--', zorder=4)

ax_mat.set_xlim(-0.6, N-0.4)
ax_mat.set_ylim(-0.6, NC-0.4)
ax_mat.set_xticks(x)
ax_mat.set_xticklabels([f"{short_ct(d['cell_type'])}\n{locus(d['marker'])}\nI²={d['i2']:.0f}%"
                        for d in leads], fontsize=7, rotation=45, ha='right')
ax_mat.set_yticks([NC-1-i for i in range(NC)])
ax_mat.set_yticklabels([COHORT_LABEL[c] for c in COHORT_ORDER], fontsize=8)
ax_mat.tick_params(length=0)
for s in ['top','right','left','bottom']:
    ax_mat.spines[s].set_visible(False)

# ancestry side labels
ax_mat.text(-2.4, NC-1-(n_eur-1)/2, 'EUR', rotation=90, va='center', ha='center',
            fontsize=9, fontweight='bold', color=ANC_CLR['EUR'])
ax_mat.text(-2.4, NC-1-n_eur-(n_afr-1)/2, 'AFR', rotation=90, va='center', ha='center',
            fontsize=9, fontweight='bold', color=ANC_CLR['AFR'])
ax_mat.text(-2.4, 0, 'LAT', rotation=90, va='center', ha='center',
            fontsize=8, fontweight='bold', color=ANC_CLR['AMR'])

# ── legend ────────────────────────────────────────────────────────────────────
leg = [
    mpatches.Patch(color=tuple(C_SAME), label='Same direction as pooled'),
    mpatches.Patch(color=tuple(C_OPP),  label='Opposite direction'),
    mpatches.Patch(color=C_MISS,        label='Absent / allele mismatch'),
    Line2D([0],[0], marker='.', color='black', ls='none', markersize=8,
           label='Cohort P < 0.05'),
]
ax_mat.legend(handles=leg, fontsize=8, loc='upper left',
              bbox_to_anchor=(1.005, 1.0), frameon=True, framealpha=0.95,
              edgecolor='none', title='Effect direction', title_fontsize=8.5)

# saturation note
ax_mat.text(1.005, 0.30, "Colour saturation ∝ |β|\n(stronger effect = more saturated)",
            transform=ax_mat.transAxes, fontsize=7.5, va='top', color='#555')

# overall concordance annotation
tot_same = sum(d['same'] for d in leads)
tot_n    = sum(d['n'] for d in leads)
ax_mat.text(1.005, 0.15,
            f"Overall: {tot_same}/{tot_n}\n= {100*tot_same/tot_n:.0f}% concordant",
            transform=ax_mat.transAxes, fontsize=9, va='top', fontweight='bold',
            color='#2166AC')

for ext in ['png','svg','pdf']:
    out = OUTD / f"concordance_matrix.{ext}"
    fig.savefig(out, dpi=300 if ext=='png' else 150, bbox_inches='tight',
                facecolor='white')
    print(f"Saved {out}")
plt.close(fig)
