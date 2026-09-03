#!/usr/bin/env python3
"""
Meta-analysis sensitivity checks for the 19 genome-wide-significant lead SNPs.

For every lead (cell type × marker) we recompute the pooled effect from the
individual-cohort β/SE (status ok/flipped) using two models and a leave-one-out
scan:

  * Fixed-effect inverse-variance weighting (IVW) — matches METAL SCHEME STDERR.
  * Random-effects DerSimonian–Laird (DL) — accounts for between-cohort variance.
  * Leave-one-cohort-out (LOO) fixed-effect scan — largest deviation & worst P.

Outputs:
  results/.../ancestry_lead_effects/lead_meta_sensitivity.tsv
"""
import csv
import math
from pathlib import Path
from collections import defaultdict

ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
COH  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_cohort_effects.tsv"
HET  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_het_stats.tsv"
OUT  = ROOT / "results/meta_sensitivity/ancestry_lead_effects/lead_meta_sensitivity.tsv"

SQRT2 = math.sqrt(2.0)

def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

def norm_sf(z):
    """Two-sided p from |z| via erfc (no scipy dependency)."""
    return math.erfc(abs(z) / SQRT2)

def fixed_effect(bs, ses):
    """Inverse-variance fixed-effect meta. Returns (beta, se, p, k)."""
    w = [1.0 / (s * s) for s in ses]
    sw = sum(w)
    beta = sum(wi * bi for wi, bi in zip(w, bs)) / sw
    se = math.sqrt(1.0 / sw)
    p = norm_sf(beta / se)
    return beta, se, p, len(bs)

def dl_random_effects(bs, ses):
    """DerSimonian–Laird random-effects meta. Returns (beta, se, p, tau2, i2)."""
    k = len(bs)
    w = [1.0 / (s * s) for s in ses]
    sw = sum(w)
    b_fe = sum(wi * bi for wi, bi in zip(w, bs)) / sw
    Q = sum(wi * (bi - b_fe) ** 2 for wi, bi in zip(w, bs))
    if k > 1:
        c = sw - sum(wi * wi for wi in w) / sw
        tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
        i2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    else:
        tau2, i2 = 0.0, 0.0
    wr = [1.0 / (s * s + tau2) for s in ses]
    swr = sum(wr)
    beta = sum(wi * bi for wi, bi in zip(wr, bs)) / swr
    se = math.sqrt(1.0 / swr)
    p = norm_sf(beta / se)
    return beta, se, p, tau2, i2

# ── load per-cohort effects ────────────────────────────────────────────────────
coh = defaultdict(list)  # (ct, marker) -> [(cohort, beta, se), ...]
for r in csv.DictReader(open(COH), delimiter='\t'):
    if r['status'] in ('ok', 'flipped'):
        b, s = f(r['beta']), f(r['se'])
        if b is not None and s is not None and s > 0:
            coh[(r['cell_type'], r['marker'])].append((r['cohort'], b, s))

# lead order + METAL pooled stats for reference
leads = []
for r in csv.DictReader(open(HET), delimiter='\t'):
    if r['stratum'] == 'Pooled':
        leads.append({'cell_type': r['cell_type'], 'marker': r['marker'],
                      'metal_beta': f(r['beta']), 'metal_p': f(r['p']),
                      'metal_i2': f(r['i2'])})

# ── compute ────────────────────────────────────────────────────────────────────
cols = ['cell_type', 'marker', 'k_cohorts',
        'fe_beta', 'fe_se', 'fe_p',
        're_beta', 're_se', 're_p', 're_tau2', 'i2',
        'loo_beta_min', 'loo_beta_max', 'loo_p_max', 'loo_worst_cohort',
        'sign_stable', 'gw_stable_fe', 'gw_stable_re', 'gw_stable_loo']
GW = 5e-8
rows = []
for L in leads:
    key = (L['cell_type'], L['marker'])
    entries = coh.get(key, [])
    if len(entries) < 2:
        continue
    cohorts = [e[0] for e in entries]
    bs = [e[1] for e in entries]
    ses = [e[2] for e in entries]

    fe_b, fe_se, fe_p, k = fixed_effect(bs, ses)
    re_b, re_se, re_p, tau2, i2 = dl_random_effects(bs, ses)

    # leave-one-cohort-out fixed-effect scan
    loo = []
    for i in range(k):
        bb = bs[:i] + bs[i+1:]
        ss = ses[:i] + ses[i+1:]
        lb, lse, lp, _ = fixed_effect(bb, ss)
        loo.append((cohorts[i], lb, lp))
    loo_bmin = min(x[1] for x in loo)
    loo_bmax = max(x[1] for x in loo)
    worst = max(loo, key=lambda x: x[2])  # largest (least significant) p
    loo_pmax = worst[2]

    ref_sign = math.copysign(1, fe_b)
    sign_stable = all(math.copysign(1, x[1]) == ref_sign for x in loo)

    rows.append({
        'cell_type': L['cell_type'], 'marker': L['marker'], 'k_cohorts': k,
        'fe_beta': f"{fe_b:.4f}", 'fe_se': f"{fe_se:.4f}", 'fe_p': f"{fe_p:.3e}",
        're_beta': f"{re_b:.4f}", 're_se': f"{re_se:.4f}", 're_p': f"{re_p:.3e}",
        're_tau2': f"{tau2:.4f}", 'i2': f"{i2:.1f}",
        'loo_beta_min': f"{loo_bmin:.4f}", 'loo_beta_max': f"{loo_bmax:.4f}",
        'loo_p_max': f"{loo_pmax:.3e}", 'loo_worst_cohort': worst[0],
        'sign_stable': 'yes' if sign_stable else 'no',
        'gw_stable_fe': 'yes' if fe_p < GW else 'no',
        'gw_stable_re': 'yes' if re_p < GW else 'no',
        'gw_stable_loo': 'yes' if loo_pmax < GW else 'no',
    })

with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
    w.writeheader()
    w.writerows(rows)

# ── console summary ─────────────────────────────────────────────────────────────
n = len(rows)
re_gw  = sum(r['gw_stable_re'] == 'yes' for r in rows)
loo_gw = sum(r['gw_stable_loo'] == 'yes' for r in rows)
sgn    = sum(r['sign_stable'] == 'yes' for r in rows)
print(f"Leads analysed:              {n}")
print(f"Sign-stable under LOO:       {sgn}/{n}")
print(f"GW-significant (random-eff): {re_gw}/{n}")
print(f"GW-significant (worst LOO):  {loo_gw}/{n}")
print(f"Wrote {OUT}")
