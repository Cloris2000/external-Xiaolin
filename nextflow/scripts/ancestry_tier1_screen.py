#!/usr/bin/env python3
"""
Tier 1 ancestry masking screen (heterogeneity-focused).

Identifies pooled meta hits that may be ancestry-masked when pooling 15 cohorts.
Candidate selection (restored from original design):
  - Genome-wide hits (p < 5e-8)
  - Suggestive hits (p < 1e-5) with high heterogeneity (I² >= 50%)
  - Coloc priority loci (chr7 TMEM106B region, etc.)

Per-cohort BETA/SE are looked up from REGENIE sumstats. Run one cell type
per SLURM array task to avoid timeouts, then merge.

Modes:
  --cell-type Astrocyte   Process a single cell type (SLURM array worker)
  --merge-only            Merge per-cell-type outputs into final TSVs

Worker outputs (per cell type):
  {output_dir}/per_celltype/{cell_type}_candidates.tsv
  {output_dir}/per_celltype/{cell_type}_detail.tsv

Merged outputs:
  ancestry_masking_candidates.tsv
  ancestry_comparison_screened_hits.tsv
  ancestry_locus_detail.tsv
"""

import argparse
import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path


CELL_TYPES = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
    "L5.ET", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", "Oligodendrocyte",
    "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC",
]

EUR_COHORTS = [
    "ROSMAP", "ROSMAP_array", "Mayo", "MSBB",
    "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GTEx_v10", "NABEC",
]
HBCC_COHORTS = ["NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650"]
MIXED_COHORTS = ["GVEX", "AMP_AD_Rush", "AMP_AD_Mayo"]
ALL_COHORTS = EUR_COHORTS + HBCC_COHORTS + MIXED_COHORTS

CANDIDATE_FIELDS = [
    "cell_type", "marker", "masking_flag", "masking_score", "is_priority",
    "pooled_p", "pooled_beta", "pooled_i2", "pooled_het_p", "pooled_is_gw",
    "eur_meta_p", "eur_meta_beta", "eur_is_gw", "eur_is_suggestive",
    "hbcc_meta_p", "hbcc_meta_beta", "hbcc_is_gw", "hbcc_is_suggestive",
    "eur_mean_beta", "hbcc_mean_beta", "eur_n_cohorts", "hbcc_n_cohorts",
    "eur_hbcc_sign_discordant", "eur_only_significant", "hbcc_only_significant",
    "high_heterogeneity", "unique_to_stratum",
]

DETAIL_FIELDS = [
    "cell_type", "marker", "masking_flag", "masking_score",
    "cohort", "ancestry_group", "beta", "se", "pooled_p", "pooled_i2",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", default=".")
    p.add_argument("--cell-type", default=None, help="Process one cell type (worker mode)")
    p.add_argument("--merge-only", action="store_true", help="Merge per-cell-type outputs")
    p.add_argument("--hit-overlap-tsv",
                   default="results/meta_sensitivity/stratified_comparison_tier1/stratified_hit_overlap.tsv")
    p.add_argument("--unique-hits-tsv",
                   default="results/meta_sensitivity/stratified_comparison_tier1/stratified_unique_hits.tsv")
    p.add_argument("--coloc-tsv", default="results/coloc/coloc_all_results.tsv")
    p.add_argument("--output-dir", default="results/meta_sensitivity/ancestry_tier1")
    p.add_argument("--gw-p-thresh", type=float, default=5e-8)
    p.add_argument("--suggestive-p-thresh", type=float, default=1e-5)
    p.add_argument("--i2-thresh", type=float, default=50.0)
    return p.parse_args()


def to_float(val):
    if val is None:
        return None
    s = str(val).strip()
    if s in ("", "NA", "nan"):
        return None
    try:
        v = float(s)
        return None if math.isnan(v) or math.isinf(v) else v
    except ValueError:
        return None


def to_int(val):
    v = to_float(val)
    return int(v) if v is not None else 0


def load_priority_prefixes(coloc_tsv):
    priority = {ct: set() for ct in CELL_TYPES}
    path = Path(coloc_tsv)
    if not path.exists():
        return priority
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lid = row.get("locus_id", "").strip()
            parts = lid.split("_", 1)
            if len(parts) == 2 and parts[0] in priority:
                priority[parts[0]].add(parts[1])
    return priority


def marker_is_priority(cell_type, marker, priority_prefixes):
    prefixes = priority_prefixes.get(cell_type, set())
    parts = marker.split(":")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}" in prefixes
    return False


def is_candidate(row, cell_type, priority_prefixes, gw_thresh, sug_thresh, i2_thresh):
    """GW, or suggestive+high I², or coloc priority."""
    marker = row["marker"].strip()
    p = to_float(row.get("pooled_p"))
    i2 = to_float(row.get("pooled_het_isq"))
    is_gw = to_int(row.get("pooled_is_gw")) or (p is not None and p <= gw_thresh)
    is_priority = marker_is_priority(cell_type, marker, priority_prefixes)
    is_sug_high_het = (
        p is not None and p <= sug_thresh
        and i2 is not None and i2 >= i2_thresh
    )
    return is_gw or is_sug_high_het or is_priority


def load_candidates_for_cell_type(cell_type, overlap_tsv, unique_tsv, priority_prefixes, args):
    candidates = {}
    with open(overlap_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["cell_type"].strip() != cell_type:
                continue
            if not is_candidate(row, cell_type, priority_prefixes, args.gw_p_thresh,
                                  args.suggestive_p_thresh, args.i2_thresh):
                continue
            marker = row["marker"].strip()
            p = to_float(row.get("pooled_p"))
            i2 = to_float(row.get("pooled_het_isq"))
            candidates[marker] = {
                "cell_type": cell_type,
                "marker": marker,
                "pooled_p": p,
                "pooled_beta": to_float(row.get("pooled_beta")),
                "pooled_i2": i2,
                "pooled_het_p": None,
                "pooled_is_gw": int(p is not None and p <= args.gw_p_thresh),
                "eur_meta_p": to_float(row.get("EUR_homogeneous_p")),
                "eur_meta_beta": to_float(row.get("EUR_homogeneous_beta")),
                "eur_is_gw": to_int(row.get("EUR_homogeneous_is_gw")),
                "eur_is_suggestive": to_int(row.get("EUR_homogeneous_is_suggestive")),
                "hbcc_meta_p": to_float(row.get("HBCC_AFR_enriched_p")),
                "hbcc_meta_beta": to_float(row.get("HBCC_AFR_enriched_beta")),
                "hbcc_is_gw": to_int(row.get("HBCC_AFR_enriched_is_gw")),
                "hbcc_is_suggestive": to_int(row.get("HBCC_AFR_enriched_is_suggestive")),
                "is_priority": int(marker_is_priority(cell_type, marker, priority_prefixes)),
                "unique_to_stratum": "",
            }

    if Path(unique_tsv).exists():
        with open(unique_tsv, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["cell_type"].strip() != cell_type:
                    continue
                marker = row["marker"].strip()
                if marker not in candidates:
                    candidates[marker] = {
                        "cell_type": cell_type, "marker": marker,
                        "pooled_p": None, "pooled_beta": None, "pooled_i2": None,
                        "pooled_het_p": None, "pooled_is_gw": 1,
                        "eur_meta_p": to_float(row.get("p")),
                        "eur_meta_beta": to_float(row.get("beta")),
                        "eur_is_gw": 1, "eur_is_suggestive": 1,
                        "hbcc_meta_p": None, "hbcc_meta_beta": None,
                        "hbcc_is_gw": 0, "hbcc_is_suggestive": 0,
                        "is_priority": int(marker_is_priority(cell_type, marker, priority_prefixes)),
                        "unique_to_stratum": row.get("stratum", ""),
                    }
                else:
                    candidates[marker]["unique_to_stratum"] = row.get("stratum", "")

    return candidates


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
        vid, beta, se, a1 = parts[2], to_float(parts[8]), to_float(parts[9]), parts[4].upper()
        if vid in markers and vid not in results and beta is not None:
            results[vid] = {"beta": beta, "se": se, "allele1": a1}
    return results


def sign(v):
    if v is None or v == 0:
        return 0
    return 1 if v > 0 else -1


def mean_beta(vals):
    v = [x for x in vals if x is not None]
    return sum(v) / len(v) if v else None


def score_candidate(c, i2_thresh):
    eur_sig = bool(c["eur_is_gw"] or c["eur_is_suggestive"])
    hbcc_sig = bool(c["hbcc_is_gw"] or c["hbcc_is_suggestive"])
    high_het = c.get("pooled_i2") is not None and c["pooled_i2"] >= i2_thresh

    eur_sign = sign(c.get("eur_mean_beta"))
    hbcc_sign = sign(c.get("hbcc_mean_beta"))
    sign_discordant = int(eur_sign != 0 and hbcc_sign != 0 and eur_sign != hbcc_sign)

    eur_only = int(eur_sig and not hbcc_sig)
    hbcc_only = int(hbcc_sig and not eur_sig)

    score = 0
    if c["is_priority"]:
        score += 5
    if high_het:
        score += 2
    if sign_discordant:
        score += 3
    if high_het and sign_discordant:
        score += 2
    if eur_only or hbcc_only:
        score += 2
    if c.get("unique_to_stratum"):
        score += 2
    if c["pooled_is_gw"]:
        score += 1

    c.update({
        "high_heterogeneity": int(high_het),
        "eur_hbcc_sign_discordant": sign_discordant,
        "eur_only_significant": eur_only,
        "hbcc_only_significant": hbcc_only,
        "masking_score": score,
        "masking_flag": "yes" if score >= 4 else "maybe" if score >= 2 else "no",
    })
    return c


def process_cell_type(cell_type, args, project_dir, out_dir):
    overlap_path = project_dir / args.hit_overlap_tsv
    if not overlap_path.exists():
        sys.exit(f"Missing {overlap_path}")

    priority_prefixes = load_priority_prefixes(project_dir / args.coloc_tsv)
    candidates = load_candidates_for_cell_type(
        cell_type, overlap_path, project_dir / args.unique_hits_tsv,
        priority_prefixes, args,
    )
    print(f"[{cell_type}] {len(candidates)} candidates (GW / suggestive+I²>={args.i2_thresh} / coloc)", file=sys.stderr)

    markers = set(candidates.keys())
    cohort_data = {m: {} for m in markers}

    for cohort in ALL_COHORTS:
        sp = project_dir / "results" / cohort / "regenie_step2" / f"{cohort}_{cell_type}_step2.regenie.raw_p"
        found = batch_grep_markers(sp, markers)
        for marker, entry in found.items():
            cohort_data[marker][cohort] = entry

    candidate_rows = []
    detail_rows = []

    for marker, c in candidates.items():
        eur_betas, hbcc_betas = [], []
        for cohort, entry in cohort_data[marker].items():
            beta = entry["beta"]
            if cohort in EUR_COHORTS:
                eur_betas.append(beta)
            elif cohort in HBCC_COHORTS:
                hbcc_betas.append(beta)

        c["eur_mean_beta"] = mean_beta(eur_betas)
        c["hbcc_mean_beta"] = mean_beta(hbcc_betas)
        c["eur_n_cohorts"] = len(eur_betas)
        c["hbcc_n_cohorts"] = len(hbcc_betas)
        c = score_candidate(c, args.i2_thresh)
        candidate_rows.append(c)

        if c["masking_flag"] in ("yes", "maybe"):
            for cohort, entry in cohort_data[marker].items():
                grp = "EUR" if cohort in EUR_COHORTS else (
                    "AFR_enriched" if cohort in HBCC_COHORTS else "mixed")
                detail_rows.append({
                    "cell_type": cell_type, "marker": marker,
                    "masking_flag": c["masking_flag"], "masking_score": c["masking_score"],
                    "cohort": cohort, "ancestry_group": grp,
                    "beta": entry["beta"], "se": entry["se"],
                    "pooled_p": c["pooled_p"], "pooled_i2": c["pooled_i2"],
                })

    candidate_rows.sort(key=lambda r: (-r["masking_score"], r["pooled_p"] or 1))

    ct_dir = out_dir / "per_celltype"
    ct_dir.mkdir(parents=True, exist_ok=True)

    def write_tsv(path, rows, fields):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    write_tsv(ct_dir / f"{cell_type}_candidates.tsv", candidate_rows, CANDIDATE_FIELDS)
    write_tsv(ct_dir / f"{cell_type}_detail.tsv", detail_rows, DETAIL_FIELDS)

    n_flagged = sum(1 for r in candidate_rows if r["masking_flag"] in ("yes", "maybe"))
    print(f"[{cell_type}] flagged {n_flagged}/{len(candidate_rows)}", file=sys.stderr)


def merge_outputs(out_dir):
    ct_dir = out_dir / "per_celltype"
    if not ct_dir.exists():
        sys.exit(f"Missing {ct_dir} — run array workers first")

    all_candidates = []
    all_detail = []

    for ct in CELL_TYPES:
        cand_path = ct_dir / f"{ct}_candidates.tsv"
        det_path = ct_dir / f"{ct}_detail.tsv"
        if not cand_path.exists():
            print(f"WARN: missing {cand_path}", file=sys.stderr)
            continue
        with open(cand_path, newline="") as fh:
            all_candidates.extend(list(csv.DictReader(fh, delimiter="\t")))
        if det_path.exists():
            with open(det_path, newline="") as fh:
                all_detail.extend(list(csv.DictReader(fh, delimiter="\t")))

    all_candidates.sort(key=lambda r: (-int(r.get("masking_score", 0)),
                                       to_float(r.get("pooled_p")) or 1))
    flagged = [r for r in all_candidates if r.get("masking_flag") in ("yes", "maybe")]

    def write_tsv(path, rows, fields):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    write_tsv(out_dir / "ancestry_comparison_screened_hits.tsv", all_candidates, CANDIDATE_FIELDS)
    write_tsv(out_dir / "ancestry_masking_candidates.tsv", flagged, CANDIDATE_FIELDS)
    write_tsv(out_dir / "ancestry_locus_detail.tsv", all_detail, DETAIL_FIELDS)

    print(f"Merged {len(all_candidates)} screened hits, {len(flagged)} flagged", file=sys.stderr)
    print(f"Wrote {out_dir}", file=sys.stderr)


def main():
    args = parse_args()
    project_dir = Path(args.project_dir).resolve()
    out_dir = project_dir / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.merge_only:
        merge_outputs(out_dir)
        return

    if not args.cell_type:
        sys.exit("Specify --cell-type for worker mode or --merge-only")

    if args.cell_type not in CELL_TYPES:
        sys.exit(f"Unknown cell type: {args.cell_type}")

    process_cell_type(args.cell_type, args, project_dir, out_dir)


if __name__ == "__main__":
    main()
