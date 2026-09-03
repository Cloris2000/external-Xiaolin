#!/usr/bin/env python3
"""
Three-panel ancestry heterogeneity figure v4.

Same as v3 except Panel B is a per-lead FOREST plot (β ± 95% CI for each
ancestry stratum) instead of the colour-coded dot plot — so the reader sees the
statistical uncertainty of each ancestry estimate.

  A (top)    : across-ancestry I² bar chart per lead SNP (shared x with B).
  B (bottom) : per-lead forest — Meta / EUR / AFR / LAT-AMR shown as points with
               vertical 95% CI whiskers (y = effect size β). Strata are colour
               coded and slightly x-offset within each lead's column.
  C          : showcase per-cohort forest plots for 6 representative loci.

Outputs: results/.../figures/ancestry_het_figure_v4.{png,svg,pdf}
(v3 is left untouched.)
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
COH  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_cohort_effects.tsv"
RSID = ROOT / "results/meta_sensitivity/ancestry_lead_effects/lead_rsids.tsv"
OUTD = ROOT / "results/meta_sensitivity/ancestry_lead_effects/figures"
OUTD.mkdir(parents=True, exist_ok=True)

# ── helpers ───────────────────────────────────────────────────────────────────
def f(v):
    try: return float(v)
    except: return None

def short_ct(ct):
    return ct.replace('Oligodendrocyte', 'Oligo')

def locus_label(marker):
    """chr7:12,285,140 → 7:12.3M"""
    p = marker.split(':')
    chrom = p[0].replace('chr', '')
    pos   = int(p[1])
    return f"{chrom}:{pos/1e6:.1f}M"

# ── load data ─────────────────────────────────────────────────────────────────
eff_rows = {(r['cell_type'], r['marker'], r['stratum']): r
            for r in csv.DictReader(open(EFF), delimiter='\t')}
het_all  = list(csv.DictReader(open(HET), delimiter='\t'))

# marker → rsID lookup (built by scripts/fetch_lead_rsids.py)
RSID_MAP = {}
if RSID.exists():
    RSID_MAP = {r['marker']: r['rsid']
                for r in csv.DictReader(open(RSID), delimiter='\t')}
def rsid_of(mk):
    rs = RSID_MAP.get(mk, '.')
    return rs if rs and rs != '.' else mk.split(':')[0].replace('chr', '') + ':' + mk.split(':')[1]

# Ordered lead list (Pooled rows only, preserving file order)
seen, leads = set(), []
for r in het_all:
    k = (r['cell_type'], r['marker'])
    if k not in seen and r['stratum'] == 'Pooled':
        seen.add(k)
        leads.append({'cell_type': r['cell_type'],
                      'marker':    r['marker'],
                      'i2':        f(r['i2']),
                      'het_p':     f(r['het_p']),
                      'beta':      f(r['beta']),
                      'p':         f(r['p'])})

N  = len(leads)   # 19 SNP columns
x  = np.arange(N)

# ── colours ───────────────────────────────────────────────────────────────────
HET_LOW  = "#BDBDBD"
HET_MED  = "#F4A460"
HET_HIGH = "#B22222"

def i2_colour(v):
    if v is None: return HET_LOW
    if v >= 75:   return HET_HIGH
    if v >= 30:   return HET_MED
    return HET_LOW

# ── per-cohort effects (harmonized) for the forest panel ──────────────────────
coh_rows = list(csv.DictReader(open(COH), delimiter='\t'))
coh_full = {}   # (ct, mk) -> list of dict(cohort, anc, beta, se)
for r in coh_rows:
    if r['beta'] and r['status'] in ('ok', 'flipped'):
        coh_full.setdefault((r['cell_type'], r['marker']), []).append(
            {'cohort': r['cohort'], 'anc': r['ancestry'],
             'beta': f(r['beta']), 'se': f(r['se'])})

# cohort display order + ancestry colours
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
META_CLR = {'EUR':'#2166AC','AFR':'#D95F02','Pooled':'#333333'}
# Panel C shows each cohort's full-sample contribution to the pooled meta; some
# cohorts (AMP-AD diverse) are ancestry-mixed, so cohorts use one neutral colour.
COH_CLR = '#5B7FA6'

# ── Panel B strata (forest) ───────────────────────────────────────────────────
STRATA      = ['Pooled', 'EUR', 'AFR', 'AMR']
STRATA_LBL  = {'Pooled': 'Meta', 'EUR': 'EUR', 'AFR': 'AFR', 'AMR': 'LAT/AMR'}
STRATA_CLR  = {'Pooled': '#333333', 'EUR': '#2166AC', 'AFR': '#D95F02', 'AMR': '#1B9E77'}
STRATA_OFF  = {'Pooled': -0.27, 'EUR': -0.09, 'AFR': 0.09, 'AMR': 0.27}

def is_snp(mk):
    p = mk.split(':'); return len(p) >= 4 and all(len(a) == 1 for a in p[2:])

# ── showcase loci: SNP leads chosen to represent under-sampled cohorts ─────────
TARGET_GTEX    = 'GTEx_v10'
TARGET_GVEX    = 'GVEX'
TARGET_DIVERSE = {'AMP_AD_Rush', 'AMP_AD_Mayo'}

def coh_set(l):
    return {r['cohort'] for r in coh_full.get((l['cell_type'], l['marker']), [])}

def show_score(l):
    cs = coh_set(l)
    cov = (2 * (TARGET_GTEX in cs)
           + 2 * (len(cs & TARGET_DIVERSE) > 0)
           + 1 * (TARGET_GVEX in cs))
    return (cov, len(cs))

best_by_mk = {}
for l in leads:
    if not is_snp(l['marker']):
        continue
    mk = l['marker']
    if mk not in best_by_mk or len(coh_set(l)) > len(coh_set(best_by_mk[mk])):
        best_by_mk[mk] = l
showcase = sorted(best_by_mk.values(), key=lambda l: (-show_score(l)[0], -show_score(l)[1]))[:6]
showcase.sort(key=lambda l: (int(l['marker'].split(':')[0].replace('chr','')),
                             int(l['marker'].split(':')[1])))

# ── across-ancestry heterogeneity (Cochran's Q / I² over EUR/AFR/LAT estimates) ─
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
        return None, k
    w = [1/se**2 for _, se in est]
    bbar = sum(wi*bi for (bi, _), wi in zip(est, w)) / sum(w)
    Q = sum(wi*(bi-bbar)**2 for (bi, _), wi in zip(est, w))
    df = k - 1
    i2 = max(0.0, (Q-df)/Q*100) if Q > 0 else 0.0
    return i2, k

anc_i2 = [across_ancestry_i2(l['cell_type'], l['marker']) for l in leads]

# ── column x-axis labels (shared A+B) ────────────────────────────────────────
col_labels = [f"{short_ct(r['cell_type'])}\n{locus_label(r['marker'])}"
              for r in leads]

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 16))
fig.patch.set_facecolor('white')

outer = gridspec.GridSpec(
    2, 1, figure=fig, height_ratios=[1.25, 2.5], hspace=0.28,
    left=0.075, right=0.90, top=0.96, bottom=0.055)

top_gs = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=outer[0], height_ratios=[1.6, 2.2], hspace=0.10)
ax_a = fig.add_subplot(top_gs[0])
ax_b = fig.add_subplot(top_gs[1], sharex=ax_a)

bot_gs = gridspec.GridSpecFromSubplotSpec(
    2, 3, subplot_spec=outer[1], hspace=0.30, wspace=0.42)
ax_c = [fig.add_subplot(bot_gs[i]) for i in range(6)]

# ── Panel A – across-ANCESTRY I² heterogeneity bars ───────────────────────────
for j, (i2v, k) in enumerate(anc_i2):
    if i2v is None:
        continue
    ax_a.bar(j, i2v, width=0.7, color=i2_colour(i2v), linewidth=0, zorder=2)

ax_a.set_xlim(-0.6, N - 0.4)
ax_a.set_ylim(0, 105)
ax_a.set_ylabel("Across-ancestry\nheterogeneity (I², %)", fontsize=9, labelpad=6)
ax_a.tick_params(axis='y', labelsize=8)
ax_a.tick_params(axis='x', bottom=False, labelbottom=False)  # shared with B
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(-0.065, 1.02, 'A', transform=ax_a.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')

leg_a_cols = [
    mpatches.Patch(color=HET_LOW,  label='I² < 30%'),
    mpatches.Patch(color=HET_MED,  label='30–75%'),
    mpatches.Patch(color=HET_HIGH, label='≥ 75%'),
]
ax_a.legend(handles=leg_a_cols, fontsize=8, loc='upper right',
            frameon=True, framealpha=0.9, edgecolor='none', handlelength=1.0)

# ── Panel B – per-lead ancestry forest (β ± 95% CI) ───────────────────────────
ax_b.axhline(0, color='grey', lw=0.8, ls='--', zorder=1)
all_lo, all_hi = [], []
for j, lead in enumerate(leads):
    for s in STRATA:
        e = eff_rows.get((lead['cell_type'], lead['marker'], s))
        if e is None:
            continue
        b = f(e['beta'])
        if b is None:
            continue
        lo = f(e.get('ci_lo')); hi = f(e.get('ci_hi'))
        if lo is None or hi is None:
            se = f(e['se'])
            if se is None:
                continue
            lo, hi = b - 1.96 * se, b + 1.96 * se
        xx = j + STRATA_OFF[s]
        ax_b.errorbar(xx, b, yerr=[[b - lo], [hi - b]],
                      fmt='D' if s == 'Pooled' else 'o',
                      color=STRATA_CLR[s],
                      markersize=4.8 if s == 'Pooled' else 4.0,
                      markeredgewidth=0, elinewidth=1.0, capsize=0, zorder=3)
        all_lo.append(lo); all_hi.append(hi)

lo_b = min(all_lo); hi_b = max(all_hi)
pad  = 0.06 * (hi_b - lo_b)
ax_b.set_xlim(-0.6, N - 0.4)
ax_b.set_ylim(lo_b - pad, hi_b + pad)
ax_b.set_ylabel("Effect size (β, 95% CI)", fontsize=9, labelpad=6)
ax_b.tick_params(axis='y', labelsize=8)
ax_b.set_xticks(x)
ax_b.set_xticklabels(col_labels, fontsize=7.5, rotation=45, ha='right')
ax_b.tick_params(axis='x', length=3)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.text(-0.065, 1.04, 'B', transform=ax_b.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')

# faint vertical guides between lead columns
for j in range(N - 1):
    ax_b.axvline(j + 0.5, color='#EEEEEE', lw=0.5, zorder=0)

leg_b = [Line2D([0], [0], marker='D' if s == 'Pooled' else 'o',
                color=STRATA_CLR[s], ls='none', markersize=6, label=STRATA_LBL[s])
         for s in STRATA]
ax_b.legend(handles=leg_b, fontsize=8, loc='upper left',
            bbox_to_anchor=(1.005, 1.0), frameon=True, framealpha=0.9,
            edgecolor='none', title='Ancestry', title_fontsize=8)

# ── Panel C – showcase per-cohort forest plots ────────────────────────────────
allb = []
for lead in showcase:
    for rec in coh_full.get((lead['cell_type'], lead['marker']), []):
        allb.append(abs(rec['beta']) + 1.96 * rec['se'])
XR = (np.percentile(allb, 98) if allb else 0.6) * 1.05

for pi, lead in enumerate(showcase):
    ax = ax_c[pi]
    ct, mk = lead['cell_type'], lead['marker']
    recs = sorted(coh_full.get((ct, mk), []), key=lambda r: COHORT_RANK.get(r['cohort'], 99))

    rows = []   # (label, beta, se, colour, is_meta)
    for s, lbl in [('Pooled', 'Meta')]:
        e = eff_rows.get((ct, mk, s))
        if e is None: continue
        b = f(e['beta']); se = f(e['se'])
        if b is None or se is None: continue
        rows.append((lbl, b, se, META_CLR[s], True))
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
                    markersize=6.5 if is_meta else 5.5, markeredgewidth=0,
                    linewidth=1.6 if is_meta else 1.1,
                    capsize=0, zorder=3)

    ax.set_xlim(-XR, XR)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.5)
    ax.set_xticks([t for t in ax.get_xticks() if -XR <= t <= XR])
    ax.tick_params(axis='x', labelsize=9.5)
    ax.set_xlabel("β (95% CI)", fontsize=10.5)
    ax.set_title(f"{short_ct(ct)}  ({rsid_of(mk)})",
                 fontsize=11.5, fontweight='bold', pad=5, loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Panel C label — aligned with the 'A'/'B' labels (same figure x-position)
LABEL_X = 0.075 + (-0.065) * (0.90 - 0.075)
_pos_c0 = ax_c[0].get_position()
fig.text(LABEL_X, _pos_c0.y1 + 0.008, 'C',
         fontsize=14, fontweight='bold', va='bottom', ha='left')

leg_c = [
    Line2D([0],[0], marker='o', color=COH_CLR, ls='none', markersize=7, label='cohort'),
    Line2D([0],[0], marker='D', color=META_CLR['Pooled'], ls='none', markersize=7, label='meta-analysis'),
]
fig.legend(handles=leg_c, fontsize=9.5, loc='lower center',
           frameon=True, framealpha=0.9, edgecolor='none',
           bbox_to_anchor=(0.5, 0.008), ncol=2)

# ── export ────────────────────────────────────────────────────────────────────
for ext in ['png', 'svg', 'pdf']:
    out = OUTD / f"ancestry_het_figure_v4.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")

plt.close(fig)
