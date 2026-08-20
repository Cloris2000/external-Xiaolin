#!/usr/bin/env python3
"""
Three-panel ancestry heterogeneity figure v3.

Left column (Panels A + B, merged / shared x-axis):
  A (top)    : I² bar chart per lead SNP.
               Y-axis = "Heterogeneity (I², %)".
               I²=30/75% reference lines annotated OUTSIDE the plot area.
               Significance stars for Q-test; star legend inside panel.
  B (bottom) : 4-row beta heatmap [Pooled | EUR | AFR | LAT/AMR].
               Diverging RdBu_r colour: red = negative, blue = positive, grey = missing.
               X-axis tick labels (SNP locus) shown only here (shared with A).

Right column (Panel C — forest plots):
  Row labels = "CellType  chr:pos" to avoid duplicate cell-type names.
  Shows EUR / AFR / Pooled beta ± 95% CI.

Panel B and C are NOT redundant:
  B = colour-coded magnitude overview (no uncertainty).
  C = beta ± 95% CI with error bars (statistical precision).
"""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable

# ── paths ─────────────────────────────────────────────────────────────────────
ROOT = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")
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
CLR_F    = {'EUR': '#2166AC', 'AFR': '#D95F02', 'Pooled': '#333333'}

def i2_colour(v):
    if v is None: return HET_LOW
    if v >= 75:   return HET_HIGH
    if v >= 30:   return HET_MED
    return HET_LOW

# Diverging colourmap for heatmap  (red = negative, blue = positive)
CMAP = plt.get_cmap('RdBu')   # blue=high, red=low  → perfect for beta
VMAX = 0.42                    # symmetric range
VNORM = mcolors.TwoSlopeNorm(vcenter=0, vmin=-VMAX, vmax=VMAX)
CLR_MISSING = "#D4D4D4"        # grey for absent data

# ── heatmap data matrix ───────────────────────────────────────────────────────
HMAP_STRATA  = ['Pooled', 'EUR', 'AFR', 'AMR']
HMAP_LABELS  = ['Meta', 'EUR', 'AFR', 'LAT/AMR']
NR           = len(HMAP_STRATA)  # 4 rows

beta_mat = np.full((NR, N), np.nan)
for j, lead in enumerate(leads):
    for i, s in enumerate(HMAP_STRATA):
        e = eff_rows.get((lead['cell_type'], lead['marker'], s))
        if e is None: continue
        bv = f(e['beta'])
        if bv is not None:
            beta_mat[i, j] = bv

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
ANC_CLR = {'EUR':'#2166AC','AFR':'#D95F02','AMR':'#1B9E77'}
META_CLR = {'EUR':'#2166AC','AFR':'#D95F02','Pooled':'#333333'}
# Panel C shows each cohort's full-sample contribution to the pooled meta; some
# cohorts (AMP-AD diverse) are ancestry-mixed, so cohorts use one neutral colour.
COH_CLR = '#5B7FA6'

def is_snp(mk):
    p = mk.split(':'); return len(p) >= 4 and all(len(a) == 1 for a in p[2:])

# ── showcase loci: SNP leads chosen to represent under-sampled cohorts ─────────
# Priority: guarantee GTEx_v10, GVEX and the AMP-AD diverse cohorts (Rush AFR,
# Mayo-LAT AMR) are represented, then fall back to best-powered loci.
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

# dedupe by marker (a variant can lead in several cell types) keeping most cohorts
best_by_mk = {}
for l in leads:
    if not is_snp(l['marker']):
        continue
    mk = l['marker']
    if mk not in best_by_mk or len(coh_set(l)) > len(coh_set(best_by_mk[mk])):
        best_by_mk[mk] = l
showcase = sorted(best_by_mk.values(), key=lambda l: (-show_score(l)[0], -show_score(l)[1]))[:6]
# display order: by chromosome/position for a stable arrangement
showcase.sort(key=lambda l: (int(l['marker'].split(':')[0].replace('chr','')),
                             int(l['marker'].split(':')[1])))

# ── across-ancestry heterogeneity (Cochran's Q / I² over EUR/AFR/LAT estimates) ─
def across_ancestry_i2(ct, mk):
    """Fixed-effect Q across the ancestry meta point estimates (EUR/AFR/LAT)."""
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
# A (I² bars) + B (ancestry β dots) stacked & sharing x on top;
# C = grid of showcase per-cohort forest plots at the bottom.
fig = plt.figure(figsize=(14, 16))
fig.patch.set_facecolor('white')

outer = gridspec.GridSpec(
    2, 1, figure=fig, height_ratios=[1.1, 2.5], hspace=0.26,
    left=0.075, right=0.93, top=0.96, bottom=0.055)

top_gs = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=outer[0], height_ratios=[2.4, 1.3], hspace=0.10)
ax_a = fig.add_subplot(top_gs[0])
ax_b = fig.add_subplot(top_gs[1], sharex=ax_a)

bot_gs = gridspec.GridSpecFromSubplotSpec(
    2, 3, subplot_spec=outer[1], hspace=0.30, wspace=0.42)
ax_c = [fig.add_subplot(bot_gs[i]) for i in range(6)]

# ── Panel A – across-ANCESTRY I² heterogeneity bars ───────────────────────────
# I² computed across the EUR/AFR/LAT meta estimates (aligns with Panel B).
for j, (i2v, k) in enumerate(anc_i2):
    if i2v is None:
        # only one ancestry available → cannot compute; leave the bar empty
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

# ── Panel B – ancestry β dots ─────────────────────────────────────────────────
# Rows top→bottom: Meta, EUR, AFR, LAT/AMR (row i drawn at y = NR-1-i).
DOT_SIZE = 240
for i in range(NR):
    yy = NR - 1 - i
    for j in range(N):
        bv = beta_mat[i, j]
        if np.isnan(bv):
            ax_b.scatter(j, yy, s=70, facecolors='none',
                         edgecolors=CLR_MISSING, linewidths=1.0, marker='o', zorder=2)
        else:
            ax_b.scatter(j, yy, s=DOT_SIZE, c=[CMAP(VNORM(bv))],
                         marker='o', linewidths=0.4, edgecolors='#666', zorder=3)

ax_b.set_xlim(-0.6, N - 0.4)
ax_b.set_ylim(-0.6, NR - 0.4)
ax_b.axhline(NR - 1.5, color='black', lw=1.1, zorder=4)   # line under Meta

ax_b.set_xticks(x)
ax_b.set_xticklabels(col_labels, fontsize=7.5, rotation=45, ha='right')
ax_b.tick_params(axis='x', length=3)
ax_b.set_yticks([NR - 1 - i for i in range(NR)])
ax_b.set_yticklabels(HMAP_LABELS, fontsize=9)
ax_b.tick_params(axis='y', length=0)
for s in ['top', 'right', 'left', 'bottom']:
    ax_b.spines[s].set_visible(False)
ax_b.text(-0.065, 1.06, 'B', transform=ax_b.transAxes,
          fontsize=14, fontweight='bold', va='bottom', ha='left')

# Colourbar for beta — vertical, on the right-hand side of Panel B
sm = ScalarMappable(cmap=CMAP, norm=VNORM); sm.set_array([])
cax = ax_b.inset_axes([1.015, 0.05, 0.013, 0.9])
cb  = fig.colorbar(sm, cax=cax)
cb.set_label('Effect size (β)', fontsize=9, labelpad=3)
cb.ax.tick_params(labelsize=8)

# ── Panel C – showcase per-cohort forest tables ───────────────────────────────
# Each mini-panel: left = forest (β ± 95% CI, box area ∝ inverse-variance weight),
# right = numeric columns  β (95% CI)  and  Weight (% of pooled meta, SCHEME STDERR).
allb = []
for lead in showcase:
    for rec in coh_full.get((lead['cell_type'], lead['marker']), []):
        allb.append(abs(rec['beta']) + 1.96 * rec['se'])
XR = (np.percentile(allb, 98) if allb else 0.6) * 1.05

def ci_txt(b, se):
    return f"{b:+.2f} ({b-1.96*se:+.2f}, {b+1.96*se:+.2f})"

for pi, lead in enumerate(showcase):
    ax = ax_c[pi]
    ct, mk = lead['cell_type'], lead['marker']
    recs = sorted(coh_full.get((ct, mk), []), key=lambda r: COHORT_RANK.get(r['cohort'], 99))

    # inverse-variance weights across the contributing cohorts (= pooled meta weights)
    inv = {id(rec): (1.0 / rec['se']**2 if rec['se'] and rec['se'] > 0 else 0.0)
           for rec in recs}
    wsum = sum(inv.values()) or 1.0

    # rows: cohorts (top) then the pooled meta summary (bottom)
    rows = []   # (label, beta, se, colour, is_meta, weight_pct)
    for s, lbl in [('Pooled', 'Meta')]:
        e = eff_rows.get((ct, mk, s))
        if e is None: continue
        b = f(e['beta']); se = f(e['se'])
        if b is None or se is None: continue
        rows.append((lbl, b, se, META_CLR[s], True, None))
    for rec in reversed(recs):
        rows.append((COHORT_LABEL.get(rec['cohort'], rec['cohort']),
                     rec['beta'], rec['se'], COH_CLR,
                     False, 100.0 * inv[id(rec)] / wsum))

    yy = np.arange(len(rows))
    ax.axvline(0, color='grey', lw=0.6, ls='--', zorder=1)
    n_meta = sum(1 for r in rows if r[4])
    if n_meta:
        ax.axhline(n_meta - 0.5, color='#CCCCCC', lw=0.8, zorder=1)

    # uniform circles for cohorts; small diamond for meta
    for yi, (lbl, b, se, col, is_meta, wpct) in zip(yy, rows):
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
    # left-aligned title over the forest region so it never overlaps the CI column
    ax.set_title(f"{short_ct(ct)}  ({rsid_of(mk)})",
                 fontsize=11.5, fontweight='bold', pad=5, loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Panel C label — aligned with the 'A'/'B' labels (same figure x-position)
LABEL_X = 0.075 + (-0.065) * (0.93 - 0.075)   # matches ax_a/ax_b transAxes x=-0.065
_pos_c0 = ax_c[0].get_position()
fig.text(LABEL_X, _pos_c0.y1 + 0.008, 'C',
         fontsize=14, fontweight='bold', va='bottom', ha='left')

# Shared legend at the bottom of Panel C
leg_c = [
    Line2D([0],[0], marker='o', color=COH_CLR, ls='none', markersize=7, label='cohort'),
    Line2D([0],[0], marker='D', color=META_CLR['Pooled'], ls='none', markersize=7, label='meta-analysis'),
]
fig.legend(handles=leg_c, fontsize=9.5, loc='lower center',
           frameon=True, framealpha=0.9, edgecolor='none',
           bbox_to_anchor=(0.5, 0.008), ncol=2)

# ── export ────────────────────────────────────────────────────────────────────
for ext in ['png', 'svg', 'pdf']:
    out = OUTD / f"ancestry_het_figure_v3.{ext}"
    fig.savefig(out, dpi=300 if ext == 'png' else 150,
                bbox_inches='tight', facecolor='white')
    print(f"Saved {out}")

plt.close(fig)
