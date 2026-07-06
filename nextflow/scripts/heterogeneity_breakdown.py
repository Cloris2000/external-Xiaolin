#!/usr/bin/env python3
"""
Heterogeneity breakdown across cohort, diagnosis, ancestry, and analysis-type groups.

For each top hit (GW p<5e-8 or suggestive p<1e-5) in the meta-analysis METAL .tbl files,
this script:
  1. Extracts per-cohort BETA/SE from individual REGENIE sumstat files.
  2. Computes within-group Cochran's Q and I² for each grouping dimension:
       - diagnosis_group  (AD_neurological / psychiatric_mixed / neurologically_normal)
       - ancestry_class   (EUR_homogeneous / AFR_enriched / mixed_ancestry)
       - analysis_type    (derived from tissue_filter + build columns)
  3. Outputs three TSVs:
       per_cohort_effects.tsv       -- raw per-cohort betas for all top hits (forest plots)
       group_i2_summary.tsv         -- group-level I²/Q per cell type × variant × dimension
       overall_group_i2.tsv         -- genome-wide I² distribution per group (% I²>50, median)

Parallel modes (recommended):
  --cell-type Astrocyte   Process one cell type (SLURM array worker)
  --merge-only            Merge per-cell-type outputs into final TSVs
"""

import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


CELL_TYPES = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
    "L5.ET", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", "Oligodendrocyte",
    "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC",
]

EFFECT_FIELDS = [
    "cell_type", "marker", "meta_P_value", "meta_HetISq", "is_gw",
    "cohort", "beta", "se", "n", "effect_allele",
    "diagnosis_group", "ancestry_class", "analysis_type",
]

GROUP_I2_FIELDS = [
    "cell_type", "marker", "P_value", "is_gw", "meta_HetISq",
    "dimension", "group", "n_cohorts_in_group",
    "within_group_I2", "within_group_Q", "within_group_HetPVal",
]

OVERALL_FIELDS = [
    "cell_type", "dimension", "group",
    "n_suggestive_variants", "median_HetISq", "pct75_HetISq", "pct_variants_I2_ge50",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cell-type", default=None, help="Process one cell type (worker mode)")
    p.add_argument("--merge-only", action="store_true", help="Merge per-cell-type outputs")
    p.add_argument("--meta-tbls", nargs="+", default=None, help="METAL .tbl files (one per cell type)")
    p.add_argument("--meta-tbl-dir", default=None, help="Directory with METAL .tbl files (worker mode)")
    p.add_argument("--cohort-sumstat-pattern", default=None,
                   help="Path pattern with {cohort} and {cell_type} placeholders")
    p.add_argument("--cohort-metadata", default=None,
                   help="docs/meta_cohort_metadata.tsv")
    p.add_argument("--output-dir", required=True, help="Directory for output TSVs")
    p.add_argument("--cohorts", nargs="+", default=None,
                   help="Ordered cohort list matching METAL Direction string")
    p.add_argument("--gw-p-thresh", type=float, default=5e-8)
    p.add_argument("--suggestive-p-thresh", type=float, default=1e-5)
    p.add_argument("--sumstat-results-dir", default="results",
                   help="Root results directory for resolving {cohort} pattern")
    return p.parse_args()


def to_float(val):
    if val is None:
        return None
    s = str(val).strip()
    if s == "" or s.upper() in ("NA", "NAN", "INF", "-INF"):
        return None
    try:
        v = float(s)
        return None if (math.isnan(v) or math.isinf(v)) else v
    except ValueError:
        return None


def infer_cell_type(path_str):
    name = Path(path_str).name
    suffix = "_meta_analysis_"
    if suffix in name:
        return name.split(suffix, 1)[0]
    return name.replace(".tbl", "")


def find_tbl_for_cell_type(cell_type, meta_tbl_dir, meta_tbls=None):
    if meta_tbls:
        for path in meta_tbls:
            if infer_cell_type(path) == cell_type and not str(path).endswith("1.tbl"):
                return str(path)
    if meta_tbl_dir:
        for path in sorted(Path(meta_tbl_dir).glob(f"{cell_type}_meta_analysis_*.tbl")):
            if not path.name.endswith("1.tbl"):
                return str(path)
    return None


def cochran_q_i2(betas, ses):
    weights = []
    for se in ses:
        if se is None or se <= 0:
            return None, None, None, None
        weights.append(1.0 / (se * se))

    n = len(betas)
    if n < 2:
        return None, None, None, None

    w_sum = sum(weights)
    beta_bar = sum(w * b for w, b in zip(weights, betas)) / w_sum
    Q = sum(w * (b - beta_bar) ** 2 for w, b in zip(weights, betas))
    df = n - 1
    I2 = max(0.0, 100.0 * (Q - df) / Q) if Q > 0 else 0.0

    try:
        import scipy.stats as st
        pval = float(1.0 - st.chi2.cdf(Q, df))
    except ImportError:
        pval = _chi2_sf(Q, df)

    return Q, df, pval, I2


def _chi2_sf(x, k):
    if x <= 0:
        return 1.0
    a = k / 2.0
    try:
        return 1.0 - _regularised_lower_gamma(a, x / 2.0)
    except Exception:
        return float("nan")


def _regularised_lower_gamma(a, x):
    if x == 0:
        return 0.0
    lga = math.lgamma(a)
    term = math.exp(a * math.log(x) - x - lga) / a
    total = term
    for i in range(1, 300):
        term *= x / (a + i)
        total += term
        if abs(term) < 1e-12 * abs(total):
            break
    return min(total, 1.0)


def load_cohort_metadata(path):
    meta = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cohort = row["cohort"].strip()
            meta[cohort] = {k: v.strip() for k, v in row.items()}
    return meta


def derive_analysis_type(meta_row):
    tissue = meta_row.get("tissue_filter", "").lower()
    build = meta_row.get("build", "")
    platform = meta_row.get("platform_group", "").lower()

    if "frontal cortex" in tissue or "ba9" in tissue:
        if build == "hg38" and "gtex" in platform:
            return "control_brain_WGS"
        return "control_brain"
    if "temporal" in tissue or "superior temporal" in tissue:
        return "non_DLPFC_brain"
    if "dlpfc" in tissue or "prefrontal" in tissue or "dorsolateral" in tissue:
        if build == "hg38" and meta_row.get("imputation_metric", "") == "not_applicable":
            return "DLPFC_WGS"
        return "DLPFC_array_imputed"
    return "other"


def build_group_maps(cohort_metadata, cohort_list):
    dims = {"diagnosis_group": {}, "ancestry_class": {}, "analysis_type": {}}
    for cohort in cohort_list:
        m = cohort_metadata.get(cohort, {})
        dims["diagnosis_group"][cohort] = m.get("diagnosis_group", "unknown")
        dims["ancestry_class"][cohort] = m.get("ancestry_class", "unknown")
        dims["analysis_type"][cohort] = derive_analysis_type(m)
    return dims


def parse_metal_top_hits(tbl_paths, gw_thresh, sug_thresh, cell_type=None):
    top_hits = {}
    for path in tbl_paths:
        path = str(path)
        if Path(path).name.endswith("1.tbl"):
            continue
        ct = infer_cell_type(path)
        if cell_type and ct != cell_type:
            continue
        hits = []
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                p = to_float(row.get("P-value"))
                if p is None or p > sug_thresh:
                    continue
                hits.append({
                    "MarkerName": row.get("MarkerName", "").strip(),
                    "Allele1": row.get("Allele1", "").strip().upper(),
                    "Allele2": row.get("Allele2", "").strip().upper(),
                    "meta_effect": to_float(row.get("Effect")),
                    "meta_se": to_float(row.get("StdErr")),
                    "P_value": p,
                    "HetISq": to_float(row.get("HetISq")),
                    "HetPVal": to_float(row.get("HetPVal")),
                    "Direction": row.get("Direction", "").strip(),
                    "is_gw": p <= gw_thresh,
                })
        top_hits[ct] = hits
    return top_hits


def resolve_sumstat_path(cohort, cell_type, sumstat_pattern, results_dir):
    path = sumstat_pattern.format(cohort=cohort, cell_type=cell_type)
    if not os.path.isabs(path) and not path.startswith("results/"):
        path = os.path.join(results_dir, cohort, "regenie_step2",
                            f"{cohort}_{cell_type}_step2.regenie.raw_p")
    return path


def batch_grep_markers(sumstat_path, markers):
    path = Path(sumstat_path)
    results = {}
    if not path.exists() or not markers:
        return results

    with tempfile.NamedTemporaryFile("w", suffix=".markers", delete=False) as tmp:
        for m in sorted(markers):
            tmp.write(m + "\n")
        tmp_path = tmp.name

    try:
        proc = subprocess.run(
            ["grep", "-F", "-f", tmp_path, str(path)],
            capture_output=True, text=True, timeout=900,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return results
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    if proc.returncode not in (0, 1):
        return results

    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) < 10 or parts[0] == "CHROM":
            continue
        vid = parts[2]
        if vid in markers and vid not in results:
            results[vid] = {
                "beta": to_float(parts[8]),
                "se": to_float(parts[9]),
                "n": int(to_float(parts[10])) if len(parts) > 10 and to_float(parts[10]) else None,
                "allele1": parts[4].upper(),
            }
    return results


def collect_per_cohort_effects(cell_type, hits, cohort_list, sumstat_pattern, results_dir):
    effects = defaultdict(dict)
    markers = {h["MarkerName"] for h in hits}
    hit_by_marker = {h["MarkerName"]: h for h in hits}

    for cohort in cohort_list:
        sp = resolve_sumstat_path(cohort, cell_type, sumstat_pattern, results_dir)
        found = batch_grep_markers(sp, markers)
        if not found:
            print(f"  [WARN] No data for {cohort} / {cell_type}: {sp}", file=sys.stderr)
            continue
        for marker, entry in found.items():
            if entry["beta"] is None or entry["se"] is None:
                continue
            hit = hit_by_marker[marker]
            beta = entry["beta"]
            if entry["allele1"] and hit["Allele1"]:
                if entry["allele1"].upper() != hit["Allele1"].upper():
                    beta = -beta
            effects[marker][cohort] = {
                "beta": beta,
                "se": entry["se"],
                "n": entry["n"],
                "allele1": entry["allele1"],
            }
    return effects


def compute_group_i2(effects, cell_type, hits, cohort_list, group_maps):
    rows = []
    for hit in hits:
        marker = hit["MarkerName"]
        cohort_data = effects.get(marker, {})

        for dim_name, cohort_to_group in group_maps.items():
            groups = defaultdict(list)
            for cohort in cohort_list:
                group_label = cohort_to_group.get(cohort, "unknown")
                entry = cohort_data.get(cohort)
                if entry is None:
                    continue
                groups[group_label].append((entry["beta"], entry["se"]))

            for group_label, pairs in groups.items():
                betas = [b for b, s in pairs]
                ses = [s for b, s in pairs]
                Q, df, pval, I2 = cochran_q_i2(betas, ses)
                rows.append({
                    "cell_type": cell_type,
                    "marker": marker,
                    "P_value": hit["P_value"],
                    "is_gw": int(hit["is_gw"]),
                    "meta_HetISq": "" if hit["HetISq"] is None else hit["HetISq"],
                    "dimension": dim_name,
                    "group": group_label,
                    "n_cohorts_in_group": len(betas),
                    "within_group_I2": "" if I2 is None else round(I2, 2),
                    "within_group_Q": "" if Q is None else round(Q, 4),
                    "within_group_HetPVal": "" if pval is None else round(pval, 6),
                })
    return rows


def compute_overall_group_i2(tbl_path, cell_type, cohort_list, group_maps, gw_thresh, sug_thresh):
    group_i2_values = defaultdict(list)

    with open(tbl_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            p = to_float(row.get("P-value"))
            i2 = to_float(row.get("HetISq"))
            direction = row.get("Direction", "").strip()

            if i2 is None or p is None or p > sug_thresh:
                continue

            present_cohorts = []
            for idx, char in enumerate(direction):
                if char in ("+", "-") and idx < len(cohort_list):
                    present_cohorts.append(cohort_list[idx])
            if len(present_cohorts) < 2:
                continue

            for dim_name, cohort_to_group in group_maps.items():
                groups_present = defaultdict(int)
                for cohort in present_cohorts:
                    g = cohort_to_group.get(cohort, "unknown")
                    groups_present[g] += 1
                for g in groups_present:
                    group_i2_values[(cell_type, dim_name, g)].append(i2)

    rows = []
    for (ct, dim_name, group), i2_list in sorted(group_i2_values.items()):
        if not i2_list:
            continue
        n = len(i2_list)
        sorted_i2 = sorted(i2_list)
        rows.append({
            "cell_type": ct,
            "dimension": dim_name,
            "group": group,
            "n_suggestive_variants": n,
            "median_HetISq": round(sorted_i2[n // 2], 2),
            "pct75_HetISq": round(sorted_i2[int(0.75 * n)], 2),
            "pct_variants_I2_ge50": round(100.0 * sum(1 for v in i2_list if v >= 50) / n, 2),
        })
    return rows


def build_effect_rows(cell_type, hits, effects, cohort_list, cohort_metadata):
    rows = []
    for hit in hits:
        marker = hit["MarkerName"]
        for cohort in cohort_list:
            entry = effects.get(marker, {}).get(cohort)
            if entry is None:
                continue
            m = cohort_metadata.get(cohort, {})
            rows.append({
                "cell_type": cell_type,
                "marker": marker,
                "meta_P_value": hit["P_value"],
                "meta_HetISq": "" if hit["HetISq"] is None else hit["HetISq"],
                "is_gw": int(hit["is_gw"]),
                "cohort": cohort,
                "beta": entry["beta"],
                "se": entry["se"],
                "n": entry["n"] if entry["n"] is not None else "",
                "effect_allele": entry["allele1"],
                "diagnosis_group": m.get("diagnosis_group", ""),
                "ancestry_class": m.get("ancestry_class", ""),
                "analysis_type": derive_analysis_type(m),
            })
    return rows


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def process_cell_type(cell_type, args):
    tbl_path = find_tbl_for_cell_type(cell_type, args.meta_tbl_dir, args.meta_tbls)
    if not tbl_path:
        sys.exit(f"No METAL .tbl found for cell type: {cell_type}")

    cohort_metadata = load_cohort_metadata(args.cohort_metadata)
    cohort_list = args.cohorts
    group_maps = build_group_maps(cohort_metadata, cohort_list)

    top_hits = parse_metal_top_hits([tbl_path], args.gw_p_thresh, args.suggestive_p_thresh, cell_type)
    hits = top_hits.get(cell_type, [])
    print(f"[{cell_type}] {len(hits)} suggestive hits from {tbl_path}", file=sys.stderr)

    effects = collect_per_cohort_effects(
        cell_type, hits, cohort_list, args.cohort_sumstat_pattern, args.sumstat_results_dir,
    )

    effect_rows = build_effect_rows(cell_type, hits, effects, cohort_list, cohort_metadata)
    group_i2_rows = compute_group_i2(effects, cell_type, hits, cohort_list, group_maps)
    overall_rows = compute_overall_group_i2(
        tbl_path, cell_type, cohort_list, group_maps, args.gw_p_thresh, args.suggestive_p_thresh,
    )

    ct_dir = Path(args.output_dir) / "per_celltype"
    ct_dir.mkdir(parents=True, exist_ok=True)

    write_tsv(str(ct_dir / f"{cell_type}_per_cohort_effects.tsv"), effect_rows, EFFECT_FIELDS)
    write_tsv(str(ct_dir / f"{cell_type}_group_i2_summary.tsv"), group_i2_rows, GROUP_I2_FIELDS)
    write_tsv(str(ct_dir / f"{cell_type}_overall_group_i2.tsv"), overall_rows, OVERALL_FIELDS)

    print(f"[{cell_type}] wrote {len(effect_rows)} effect rows", file=sys.stderr)


def merge_outputs(output_dir):
    ct_dir = Path(output_dir) / "per_celltype"
    if not ct_dir.exists():
        sys.exit(f"Missing {ct_dir} — run array workers first")

    all_effects, all_group, all_overall = [], [], []

    for ct in CELL_TYPES:
        for suffix, collector, fields in [
            ("_per_cohort_effects.tsv", all_effects, EFFECT_FIELDS),
            ("_group_i2_summary.tsv", all_group, GROUP_I2_FIELDS),
            ("_overall_group_i2.tsv", all_overall, OVERALL_FIELDS),
        ]:
            path = ct_dir / f"{ct}{suffix}"
            if not path.exists():
                print(f"WARN: missing {path}", file=sys.stderr)
                continue
            with open(path, newline="") as fh:
                collector.extend(list(csv.DictReader(fh, delimiter="\t")))

    write_tsv(str(Path(output_dir) / "per_cohort_effects.tsv"), all_effects, EFFECT_FIELDS)
    write_tsv(str(Path(output_dir) / "group_i2_summary.tsv"), all_group, GROUP_I2_FIELDS)
    write_tsv(str(Path(output_dir) / "overall_group_i2.tsv"), all_overall, OVERALL_FIELDS)

    print(f"Merged {len(all_effects)} effect rows, {len(all_group)} group-I2 rows", file=sys.stderr)
    print(f"Wrote {output_dir}", file=sys.stderr)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if args.merge_only:
        merge_outputs(args.output_dir)
        return

    if not args.cell_type:
        sys.exit("Specify --cell-type for worker mode or --merge-only")

    if args.cell_type not in CELL_TYPES:
        sys.exit(f"Unknown cell type: {args.cell_type}")

    required = ["cohort_metadata", "cohort_sumstat_pattern", "cohorts"]
    for attr in required:
        if getattr(args, attr) is None:
            sys.exit(f"Worker mode requires --{attr.replace('_', '-')}")

    if args.meta_tbl_dir is None and not args.meta_tbls:
        sys.exit("Worker mode requires --meta-tbl-dir or --meta-tbls")

    process_cell_type(args.cell_type, args)


if __name__ == "__main__":
    main()
