#!/usr/bin/env python3
"""
Extract pooled 15-cohort lead SNPs and look up effects in ancestry strata.

For each cell type:
  1. Pull variants with P <= --p-thresh from the pooled 15-cohort METAL .tbl
  2. Greedy distance-clump to independent leads (--clump-kb)
  3. Look up each lead in:
       - ancestry EUR meta (METAL)
       - ancestry AFR meta (METAL)
       - Mayo AMR single-cohort GWAS (REGENIE; Latino / Hispanic-first AMR)
       - pooled 15-cohort meta (already known; re-emitted for plotting)
  4. Harmonize betas to the pooled Allele1 effect allele

Outputs (default: results/meta_sensitivity/ancestry_lead_effects/):
  ancestry_lead_effects.tsv   — long table for forest plots
  ancestry_leads.tsv          — one row per lead (clumped)
  ancestry_lead_summary.tsv   — per cell type: n_gw, n_leads, n_with_AFR/AMR

Run after ancestry EUR/AFR metas finish. Missing strata are allowed (NA rows).
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


CELL_TYPES = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
    "L5.ET", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", "Oligodendrocyte",
    "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC",
]

# Approximate analysis N from keep-lists (docs/ancestry_specific/ancestry_pairs_n50.tsv).
# EUR excludes Omni5M (dropped from EUR meta). Labels are for annotation only.
DEFAULT_N = {
    "EUR": 2818,   # 11 EUR cohorts in ancestry EUR meta
    "AFR": 270,    # 4 AFR cohorts
    "AMR": 178,    # Mayo AMR only
    "Pooled": None,
}

STRATUM_ORDER = ["EUR", "AFR", "AMR", "Pooled"]
STRATUM_LABEL = {
    "EUR": "EUR meta (White)",
    "AFR": "AFR meta (African American)",
    "AMR": "Mayo AMR (Latino; single cohort)",
    "Pooled": "Pooled 15-cohort meta",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", default=".",
                   help="Nextflow project root (default: .)")
    p.add_argument("--pooled-dir", default=None,
                   help="Pooled meta dir (default: results/meta_analysis_15cohorts)")
    p.add_argument("--eur-dir", default=None,
                   help="EUR ancestry meta dir (default: results/meta_sensitivity/ancestry_specific_EUR)")
    p.add_argument("--afr-dir", default=None,
                   help="AFR ancestry meta dir (default: results/meta_sensitivity/ancestry_specific_AFR)")
    p.add_argument("--amr-pattern", default=None,
                   help="REGENIE path with {cell_type}; default: "
                        "results/AMP_AD_Mayo_AMR/regenie_step2/"
                        "AMP_AD_Mayo_AMR_{cell_type}_step2_{cell_type}.regenie")
    p.add_argument("--output-dir", default=None,
                   help="Output directory")
    p.add_argument("--cell-types", nargs="+", default=None,
                   help="Subset of cell types (default: all 19)")
    p.add_argument("--p-thresh", type=float, default=5e-8,
                   help="Pooled lead discovery threshold (default: 5e-8)")
    p.add_argument("--clump-kb", type=int, default=500,
                   help="Greedy clump window in kb (default: 500)")
    p.add_argument("--max-leads-per-ct", type=int, default=0,
                   help="Keep only top N leads per cell type after clumping (0 = all)")
    p.add_argument("--markers-tsv", default=None,
                   help="Optional TSV with columns cell_type,marker to force lead set")
    p.add_argument("--prefer-snps", action="store_true", default=True,
                   help="If lead is an indel, prefer best SNP in the clump window (default)")
    p.add_argument("--no-prefer-snps", action="store_false", dest="prefer_snps",
                   help="Keep most significant variant even if indel")
    return p.parse_args()


def to_float(val):
    if val is None:
        return None
    s = str(val).strip()
    if s == "" or s.upper() in ("NA", "NAN", "."):
        return None
    try:
        v = float(s)
        return None if (math.isnan(v) or math.isinf(v)) else v
    except ValueError:
        return None


def norm_marker(marker: str) -> str:
    """Uppercase alleles in chr:pos:a1:a2 IDs for cross-file matching."""
    parts = marker.strip().split(":")
    if len(parts) < 4:
        return marker.strip()
    chrom, pos, *alleles = parts
    return ":".join([chrom, pos] + [a.upper() for a in alleles])


def parse_chr_pos(marker: str):
    parts = marker.split(":")
    if len(parts) < 2:
        return None, None
    chrom = parts[0].replace("chr", "")
    try:
        return chrom, int(parts[1])
    except ValueError:
        return chrom, None


def find_tbl(meta_dir: Path, cell_type: str):
    """Prefer non-empty *1.tbl (METAL ANALYZE HETEROGENEITY), else largest .tbl."""
    if not meta_dir.exists():
        return None
    cands = [p for p in meta_dir.glob(f"{cell_type}_meta_analysis_*.tbl") if p.stat().st_size > 0]
    if not cands:
        return None
    cands.sort(key=lambda p: (0 if p.name.endswith("1.tbl") else 1, -p.stat().st_size))
    return cands[0]


def is_snp(marker: str) -> bool:
    parts = marker.split(":")
    if len(parts) < 4:
        return False
    return all(len(a) == 1 for a in parts[2:])


def clump_leads(hits, window_bp: int, prefer_snps: bool = True):
    """Greedy clump by P; optionally replace indel index with best SNP in-window."""
    ordered = sorted(hits, key=lambda h: (h["p"] if h["p"] is not None else 1.0, h["marker"]))
    kept = []
    by_chrom = defaultdict(list)  # chrom -> list of positions kept
    for h in ordered:
        chrom, pos = parse_chr_pos(h["marker"])
        if pos is None:
            kept.append(h)
            continue
        if any(abs(pos - p0) <= window_bp for p0 in by_chrom[chrom]):
            continue
        chosen = h
        if prefer_snps and not is_snp(h["marker"]):
            # Best SNP among hits within the same window (same chrom)
            snps = []
            for h2 in ordered:
                c2, p2 = parse_chr_pos(h2["marker"])
                if c2 == chrom and p2 is not None and abs(p2 - pos) <= window_bp and is_snp(h2["marker"]):
                    snps.append(h2)
            if snps:
                chosen = min(snps, key=lambda x: x["p"] if x["p"] is not None else 1.0)
                _, pos = parse_chr_pos(chosen["marker"])
        by_chrom[chrom].append(pos)
        kept.append(chosen)
    return kept


def load_forced_markers(path: Path):
    """Return dict cell_type -> list of marker strings."""
    out = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ct = row["cell_type"].strip()
            m = row["marker"].strip()
            if ct and m:
                out[ct].append(m)
    return out


def extract_pooled_hits(tbl: Path, p_thresh: float):
    hits = []
    with open(tbl, "rb") as fh:
        header = fh.readline().decode("utf-8", errors="replace").rstrip("\n").split("\t")
        idx = {c: header.index(c) for c in
               ["MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value"]}
        for line in fh:
            parts = line.rstrip(b"\n").split(b"\t")
            if len(parts) <= max(idx.values()):
                continue
            p = to_float(parts[idx["P-value"]].decode())
            if p is None or p > p_thresh:
                continue
            marker = parts[idx["MarkerName"]].decode().strip()
            hits.append({
                "marker": marker,
                "marker_key": norm_marker(marker),
                "allele1": parts[idx["Allele1"]].decode().strip().upper(),
                "allele2": parts[idx["Allele2"]].decode().strip().upper(),
                "freq1": to_float(parts[idx["Freq1"]].decode()),
                "beta": to_float(parts[idx["Effect"]].decode()),
                "se": to_float(parts[idx["StdErr"]].decode()),
                "p": p,
            })
    return hits


def lookup_metal_markers(tbl: Path, want_keys: set):
    """Return marker_key -> dict with alleles/beta/se/p/freq (Allele1 = effect)."""
    found = {}
    if tbl is None or not want_keys:
        return found
    with open(tbl, "rb") as fh:
        header = fh.readline().decode("utf-8", errors="replace").rstrip("\n").split("\t")
        try:
            idx = {c: header.index(c) for c in
                   ["MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value"]}
        except ValueError as exc:
            raise SystemExit(f"Bad METAL header in {tbl}: {exc}") from exc
        for line in fh:
            parts = line.rstrip(b"\n").split(b"\t")
            if len(parts) <= max(idx.values()):
                continue
            key = norm_marker(parts[idx["MarkerName"]].decode())
            if key not in want_keys:
                continue
            found[key] = {
                "allele1": parts[idx["Allele1"]].decode().strip().upper(),
                "allele2": parts[idx["Allele2"]].decode().strip().upper(),
                "freq1": to_float(parts[idx["Freq1"]].decode()),
                "beta": to_float(parts[idx["Effect"]].decode()),
                "se": to_float(parts[idx["StdErr"]].decode()),
                "p": to_float(parts[idx["P-value"]].decode()),
            }
            if len(found) == len(want_keys):
                break
    return found


def lookup_regenie_markers(path: Path, want_keys: set):
    """REGENIE: BETA is for ALLELE1; ID should match MarkerName style."""
    found = {}
    if path is None or not path.exists() or not want_keys:
        return found
    with open(path, "rb") as fh:
        header = fh.readline().decode("utf-8", errors="replace").rstrip("\n").split()
        try:
            idx = {c: header.index(c) for c in
                   ["ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"]}
        except ValueError as exc:
            raise SystemExit(f"Bad REGENIE header in {path}: {exc}") from exc
        for line in fh:
            parts = line.rstrip(b"\n").split()
            if len(parts) <= max(idx.values()):
                continue
            key = norm_marker(parts[idx["ID"]].decode())
            if key not in want_keys:
                continue
            log10p = to_float(parts[idx["LOG10P"]].decode())
            p = 10 ** (-log10p) if log10p is not None else None
            found[key] = {
                "allele1": parts[idx["ALLELE1"]].decode().strip().upper(),
                "allele2": parts[idx["ALLELE0"]].decode().strip().upper(),
                "freq1": to_float(parts[idx["A1FREQ"]].decode()),
                "beta": to_float(parts[idx["BETA"]].decode()),
                "se": to_float(parts[idx["SE"]].decode()),
                "p": p,
                "n": int(float(parts[idx["N"]].decode()))
                if to_float(parts[idx["N"]].decode()) is not None else None,
            }
            if len(found) == len(want_keys):
                break
    return found


def harmonize(rec, target_a1: str, target_a2: str):
    """Align beta/freq to target_a1 as effect allele. Returns (beta, se, p, freq, status)."""
    if rec is None:
        return None, None, None, None, "missing"
    a1, a2 = rec["allele1"], rec["allele2"]
    beta, se, p, freq = rec["beta"], rec["se"], rec["p"], rec["freq1"]
    if a1 == target_a1 and a2 == target_a2:
        return beta, se, p, freq, "ok"
    if a1 == target_a2 and a2 == target_a1:
        flip_freq = (1.0 - freq) if freq is not None else None
        flip_beta = (-beta) if beta is not None else None
        return flip_beta, se, p, flip_freq, "flipped"
    return None, None, p, None, "allele_mismatch"


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    args = parse_args()
    project = Path(args.project_dir).resolve()
    pooled_dir = Path(args.pooled_dir) if args.pooled_dir else project / "results/meta_analysis_15cohorts"
    eur_dir = Path(args.eur_dir) if args.eur_dir else project / "results/meta_sensitivity/ancestry_specific_EUR"
    afr_dir = Path(args.afr_dir) if args.afr_dir else project / "results/meta_sensitivity/ancestry_specific_AFR"
    amr_pat = args.amr_pattern or str(
        project / "results/AMP_AD_Mayo_AMR/regenie_step2/"
        "AMP_AD_Mayo_AMR_{cell_type}_step2_{cell_type}.regenie"
    )
    out_dir = Path(args.output_dir) if args.output_dir else project / "results/meta_sensitivity/ancestry_lead_effects"
    cell_types = args.cell_types or CELL_TYPES
    window_bp = args.clump_kb * 1000
    forced = load_forced_markers(Path(args.markers_tsv)) if args.markers_tsv else None

    print(f"Pooled: {pooled_dir}")
    print(f"EUR:    {eur_dir}  [{'OK' if eur_dir.exists() else 'MISSING'}]")
    print(f"AFR:    {afr_dir}  [{'OK' if afr_dir.exists() else 'MISSING'}]")
    print(f"Output: {out_dir}")

    lead_rows = []
    effect_rows = []
    summary_rows = []

    for ct in cell_types:
        pooled_tbl = find_tbl(pooled_dir, ct)
        if pooled_tbl is None:
            print(f"[{ct}] SKIP — no pooled .tbl")
            continue

        raw = []
        if forced is not None:
            want = {norm_marker(m) for m in forced.get(ct, [])}
            if not want:
                print(f"[{ct}] no forced markers")
                continue
            pooled_map = lookup_metal_markers(pooled_tbl, want)
            leads = []
            for m in forced.get(ct, []):
                key = norm_marker(m)
                rec = pooled_map.get(key)
                if rec is None:
                    print(f"  [{ct}] forced marker missing in pooled: {m}")
                    continue
                leads.append({
                    "marker": m,
                    "marker_key": key,
                    "allele1": rec["allele1"],
                    "allele2": rec["allele2"],
                    "freq1": rec["freq1"],
                    "beta": rec["beta"],
                    "se": rec["se"],
                    "p": rec["p"],
                })
        else:
            raw = extract_pooled_hits(pooled_tbl, args.p_thresh)
            leads = clump_leads(raw, window_bp, prefer_snps=args.prefer_snps)
            # Deduplicate if indel→SNP remapping collapsed loci
            seen = set()
            uniq = []
            for h in leads:
                if h["marker_key"] in seen:
                    continue
                seen.add(h["marker_key"])
                uniq.append(h)
            leads = uniq
            if args.max_leads_per_ct and len(leads) > args.max_leads_per_ct:
                leads = leads[: args.max_leads_per_ct]

        n_raw = len(raw) if forced is None else len(leads)
        print(f"[{ct}] pooled P<={args.p_thresh:g}: {n_raw} -> {len(leads)} leads  ({pooled_tbl.name})")
        if not leads:
            summary_rows.append({
                "cell_type": ct, "n_pooled_below_thresh": 0, "n_leads": 0,
                "n_with_EUR": 0, "n_with_AFR": 0, "n_with_AMR": 0,
            })
            continue

        want_keys = {h["marker_key"] for h in leads}
        eur_map = lookup_metal_markers(find_tbl(eur_dir, ct), want_keys)
        afr_map = lookup_metal_markers(find_tbl(afr_dir, ct), want_keys)
        amr_path = Path(amr_pat.format(cell_type=ct))
        amr_map = lookup_regenie_markers(amr_path, want_keys)

        n_eur = n_afr = n_amr = 0
        for rank, h in enumerate(leads, start=1):
            key = h["marker_key"]
            chrom, pos = parse_chr_pos(h["marker"])
            lead_rows.append({
                "cell_type": ct,
                "lead_rank": rank,
                "marker": h["marker"],
                "chrom": chrom,
                "pos": pos,
                "effect_allele": h["allele1"],
                "other_allele": h["allele2"],
                "pooled_freq1": h["freq1"],
                "pooled_beta": h["beta"],
                "pooled_se": h["se"],
                "pooled_p": h["p"],
            })

            stratum_recs = {
                "Pooled": {
                    "allele1": h["allele1"], "allele2": h["allele2"],
                    "freq1": h["freq1"], "beta": h["beta"], "se": h["se"], "p": h["p"],
                },
                "EUR": eur_map.get(key),
                "AFR": afr_map.get(key),
                "AMR": amr_map.get(key),
            }

            for stratum in STRATUM_ORDER:
                rec = stratum_recs[stratum]
                beta, se, p, freq, status = harmonize(rec, h["allele1"], h["allele2"])
                n_val = None
                if status in ("ok", "flipped"):
                    if stratum == "AMR" and rec is not None:
                        n_val = rec.get("n", DEFAULT_N["AMR"])
                    elif stratum in DEFAULT_N:
                        n_val = DEFAULT_N[stratum]
                    if stratum == "EUR":
                        n_eur += 1
                    elif stratum == "AFR":
                        n_afr += 1
                    elif stratum == "AMR":
                        n_amr += 1
                ci_lo = (beta - 1.96 * se) if (beta is not None and se is not None) else None
                ci_hi = (beta + 1.96 * se) if (beta is not None and se is not None) else None
                effect_rows.append({
                    "cell_type": ct,
                    "lead_rank": rank,
                    "marker": h["marker"],
                    "effect_allele": h["allele1"],
                    "other_allele": h["allele2"],
                    "stratum": stratum,
                    "stratum_label": STRATUM_LABEL[stratum],
                    "n": n_val,
                    "freq_effect": freq,
                    "beta": beta,
                    "se": se,
                    "ci_lo": ci_lo,
                    "ci_hi": ci_hi,
                    "p": p,
                    "harmonize_status": status,
                })

        summary_rows.append({
            "cell_type": ct,
            "n_pooled_below_thresh": n_raw,
            "n_leads": len(leads),
            "n_with_EUR": n_eur,
            "n_with_AFR": n_afr,
            "n_with_AMR": n_amr,
        })

    write_tsv(
        out_dir / "ancestry_leads.tsv",
        lead_rows,
        ["cell_type", "lead_rank", "marker", "chrom", "pos",
         "effect_allele", "other_allele", "pooled_freq1", "pooled_beta", "pooled_se", "pooled_p"],
    )
    write_tsv(
        out_dir / "ancestry_lead_effects.tsv",
        effect_rows,
        ["cell_type", "lead_rank", "marker", "effect_allele", "other_allele",
         "stratum", "stratum_label", "n", "freq_effect", "beta", "se", "ci_lo", "ci_hi",
         "p", "harmonize_status"],
    )
    write_tsv(
        out_dir / "ancestry_lead_summary.tsv",
        summary_rows,
        ["cell_type", "n_pooled_below_thresh", "n_leads",
         "n_with_EUR", "n_with_AFR", "n_with_AMR"],
    )
    print(f"\nWrote {len(lead_rows)} leads / {len(effect_rows)} effect rows -> {out_dir}")


if __name__ == "__main__":
    main()
