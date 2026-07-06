#!/usr/bin/env python3
"""
Compare pooled 15-cohort meta-analysis hits with stratified meta-analysis runs.

For each cell type, identifies genome-wide and suggestive hits in the pooled meta
and each stratum, then reports overlap, effect-direction concordance, and per-stratum
hit counts.

Outputs (default: results/meta_sensitivity/stratified_comparison/):
  stratified_meta_summary.tsv           — per stratum × cell type: n_gw, n_suggestive
  stratified_hit_overlap.tsv            — pooled top hits × stratum presence/significance
  stratified_effect_concordance.tsv     — overlapping hits: beta sign agreement
  stratified_unique_hits.tsv            — GW/suggestive hits unique to each stratum
"""

import argparse
import csv
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


CELL_TYPES = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
    "L5.ET", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", "Oligodendrocyte",
    "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC",
]

DEFAULT_STRATA = {
    "pooled_15cohorts": "results/meta_analysis_15cohorts",
    "AD_neurological": "results/meta_sensitivity/stratified_ad_neurological_15cohorts",
    "neurologically_normal": "results/meta_sensitivity/stratified_neurologically_normal_15cohorts",
    "DLPFC_only": "results/meta_sensitivity/stratified_dlpfc_only_15cohorts",
    "EUR_homogeneous": "results/meta_sensitivity/stratified_eur_homogeneous_15cohorts",
    "psychiatric_mixed": "results/meta_sensitivity/stratified_psychiatric_mixed_15cohorts",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", default=".",
                   help="Nextflow project root (default: current directory)")
    p.add_argument("--pooled-dir", default=None,
                   help="Pooled meta results directory (default: results/meta_analysis_15cohorts)")
    p.add_argument("--strata-tsv", default=None,
                   help="Optional TSV with columns: stratum_label, results_dir")
    p.add_argument("--output-dir", default=None,
                   help="Output directory (default: results/meta_sensitivity/stratified_comparison)")
    p.add_argument("--gw-p-thresh", type=float, default=5e-8)
    p.add_argument("--suggestive-p-thresh", type=float, default=1e-5)
    p.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1),
                   help="Parallel cell-type workers (default: min(4, n_cpus))")
    return p.parse_args()


def to_float(val):
    if val is None:
        return None
    s = str(val).strip()
    if s == "" or s.upper() in ("NA", "NAN"):
        return None
    try:
        v = float(s)
        return None if (math.isnan(v) or math.isinf(v)) else v
    except ValueError:
        return None


def find_tbl(meta_dir, cell_type):
    meta_dir = Path(meta_dir)
    if not meta_dir.exists():
        return None
    matches = sorted(
        meta_dir.glob(f"{cell_type}_meta_analysis_*.tbl"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return matches[0] if matches else None


def load_hits(meta_dir, cell_type, gw_thresh, sug_thresh):
    """Return dict marker -> {p, beta, het_isq, is_gw, is_suggestive}"""
    tbl = find_tbl(meta_dir, cell_type)
    hits = {}
    if tbl is None:
        return hits

    with open(tbl, "rb") as fh:
        header = fh.readline().decode("utf-8", errors="replace").rstrip("\n").split("\t")
        try:
            mi = header.index("MarkerName")
            pi = header.index("P-value")
            bi = header.index("Effect")
            hi = header.index("HetISq")
        except ValueError as exc:
            raise SystemExit(f"Unexpected METAL header in {tbl}: {exc}") from exc

        for line in fh:
            parts = line.rstrip(b"\n").split(b"\t")
            if len(parts) <= max(mi, pi, bi, hi):
                continue
            p = to_float(parts[pi].decode("utf-8", errors="replace"))
            if p is None or p > sug_thresh:
                continue
            marker = parts[mi].decode("utf-8", errors="replace").strip()
            beta = to_float(parts[bi].decode("utf-8", errors="replace"))
            hits[marker] = {
                "p": p,
                "beta": beta,
                "het_isq": to_float(parts[hi].decode("utf-8", errors="replace")),
                "is_gw": p <= gw_thresh,
                "is_suggestive": True,
            }
    return hits


def load_cell_type_hits(ct, strata, gw_thresh, sug_thresh):
    return ct, {
        label: load_hits(meta_dir, ct, gw_thresh, sug_thresh) if meta_dir.exists() else {}
        for label, meta_dir in strata.items()
    }


def beta_sign(beta):
    if beta is None or beta == 0:
        return 0
    return 1 if beta > 0 else -1


def load_strata(project_dir, pooled_dir, strata_tsv):
    project_dir = Path(project_dir)
    strata = {}

    if strata_tsv:
        with open(strata_tsv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                label = row["stratum_label"].strip()
                rel = row["results_dir"].strip()
                strata[label] = project_dir / rel
    else:
        for label, rel in DEFAULT_STRATA.items():
            strata[label] = project_dir / rel

    if pooled_dir:
        strata["pooled_15cohorts"] = Path(pooled_dir)
    elif "pooled_15cohorts" in strata:
        strata["pooled_15cohorts"] = project_dir / strata["pooled_15cohorts"] if not Path(strata["pooled_15cohorts"]).is_absolute() else Path(strata["pooled_15cohorts"])

    return strata


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    project_dir = Path(args.project_dir).resolve()
    output_dir = Path(args.output_dir) if args.output_dir else project_dir / "results/meta_sensitivity/stratified_comparison"
    strata = load_strata(project_dir, args.pooled_dir, args.strata_tsv)

    pooled_label = "pooled_15cohorts"
    if pooled_label not in strata:
        raise SystemExit(f"Missing pooled meta directory (expected key '{pooled_label}')")

    stratified_labels = [k for k in sorted(strata.keys()) if k != pooled_label]
    pooled_dir = strata[pooled_label]

    print("Strata to compare:")
    for label, path in strata.items():
        status = "OK" if path.exists() else "MISSING"
        print(f"  [{status}] {label}: {path}")

    # Load suggestive hits per stratum × cell type (parallelized by cell type).
    all_hits = {label: {} for label in strata}
    workers = max(1, args.workers)
    print(f"Loading suggestive hits across {len(CELL_TYPES)} cell types ({workers} workers)...")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(load_cell_type_hits, ct, strata, args.gw_p_thresh, args.suggestive_p_thresh)
            for ct in CELL_TYPES
        ]
        for fut in as_completed(futures):
            ct, per_label = fut.result()
            for label, hits in per_label.items():
                all_hits[label][ct] = hits
            print(f"  done: {ct}", flush=True)

    # Summary per stratum × cell type
    summary_rows = []
    for label in strata:
        for ct in CELL_TYPES:
            hits = all_hits[label].get(ct, {})
            n_gw = sum(1 for h in hits.values() if h["is_gw"])
            n_sug = sum(1 for h in hits.values() if h["is_suggestive"])
            summary_rows.append({
                "stratum": label,
                "cell_type": ct,
                "n_gw_hits": n_gw,
                "n_suggestive_hits": n_sug,
                "results_dir": str(strata[label]),
                "dir_exists": int(strata[label].exists()),
            })

    # Pooled suggestive hits vs strata
    overlap_rows = []
    concordance_rows = []
    unique_rows = []

    for ct in CELL_TYPES:
        pooled_hits = all_hits[pooled_label].get(ct, {})
        for marker, pooled in pooled_hits.items():
            row = {
                "cell_type": ct,
                "marker": marker,
                "pooled_p": pooled["p"],
                "pooled_beta": pooled["beta"] if pooled["beta"] is not None else "",
                "pooled_het_isq": pooled["het_isq"] if pooled["het_isq"] is not None else "",
                "pooled_is_gw": int(pooled["is_gw"]),
            }
            n_strata_gw = 0
            n_strata_sug = 0
            n_same_sign = 0
            n_with_beta = 0

            for slabel in stratified_labels:
                sh = all_hits[slabel].get(ct, {}).get(marker)
                prefix = slabel
                if sh is None:
                    row[f"{prefix}_present"] = 0
                    row[f"{prefix}_p"] = ""
                    row[f"{prefix}_beta"] = ""
                    row[f"{prefix}_is_gw"] = 0
                    row[f"{prefix}_is_suggestive"] = 0
                    continue

                row[f"{prefix}_present"] = 1
                row[f"{prefix}_p"] = sh["p"]
                row[f"{prefix}_beta"] = sh["beta"] if sh["beta"] is not None else ""
                row[f"{prefix}_is_gw"] = int(sh["is_gw"])
                row[f"{prefix}_is_suggestive"] = int(sh["is_suggestive"])

                if sh["is_gw"]:
                    n_strata_gw += 1
                if sh["is_suggestive"]:
                    n_strata_sug += 1
                if pooled["beta"] is not None and sh["beta"] is not None:
                    n_with_beta += 1
                    if beta_sign(pooled["beta"]) == beta_sign(sh["beta"]):
                        n_same_sign += 1

            row["n_strata_gw"] = n_strata_gw
            row["n_strata_suggestive"] = n_strata_sug
            row["replicated_gw_in_any_stratum"] = int(n_strata_gw > 0)
            overlap_rows.append(row)

            if n_with_beta > 0:
                concordance_rows.append({
                    "cell_type": ct,
                    "marker": marker,
                    "pooled_p": pooled["p"],
                    "pooled_is_gw": int(pooled["is_gw"]),
                    "n_strata_with_hit": sum(1 for sl in stratified_labels if all_hits[sl].get(ct, {}).get(marker)),
                    "n_strata_same_beta_sign": n_same_sign,
                    "pct_same_sign": round(100.0 * n_same_sign / n_with_beta, 2),
                    "n_strata_gw": n_strata_gw,
                })

        # Hits unique to each stratum (GW only)
        for slabel in stratified_labels:
            for marker, sh in all_hits[slabel].get(ct, {}).items():
                if not sh["is_gw"]:
                    continue
                if marker in pooled_hits and pooled_hits[marker]["is_gw"]:
                    continue
                unique_rows.append({
                    "stratum": slabel,
                    "cell_type": ct,
                    "marker": marker,
                    "p": sh["p"],
                    "beta": sh["beta"] if sh["beta"] is not None else "",
                    "het_isq": sh["het_isq"] if sh["het_isq"] is not None else "",
                    "unique_to_stratum": 1,
                    "also_gw_in_pooled": 0,
                })

    overlap_fields = ["cell_type", "marker", "pooled_p", "pooled_beta", "pooled_het_isq", "pooled_is_gw",
                      "n_strata_gw", "n_strata_suggestive", "replicated_gw_in_any_stratum"]
    for sl in stratified_labels:
        overlap_fields.extend([f"{sl}_present", f"{sl}_p", f"{sl}_beta", f"{sl}_is_gw", f"{sl}_is_suggestive"])

    write_tsv(output_dir / "stratified_meta_summary.tsv", summary_rows,
              ["stratum", "cell_type", "n_gw_hits", "n_suggestive_hits", "results_dir", "dir_exists"])
    write_tsv(output_dir / "stratified_hit_overlap.tsv", overlap_rows, overlap_fields)
    write_tsv(output_dir / "stratified_effect_concordance.tsv", concordance_rows,
              ["cell_type", "marker", "pooled_p", "pooled_is_gw", "n_strata_with_hit",
               "n_strata_same_beta_sign", "pct_same_sign", "n_strata_gw"])
    write_tsv(output_dir / "stratified_unique_hits.tsv", unique_rows,
              ["stratum", "cell_type", "marker", "p", "beta", "het_isq", "unique_to_stratum", "also_gw_in_pooled"])

    print(f"\nWrote comparison tables to: {output_dir}")


if __name__ == "__main__":
    main()
