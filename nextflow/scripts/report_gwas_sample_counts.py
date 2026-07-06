#!/usr/bin/env python3
"""
Reconstruct per-cohort sample counts for the GWAS pipeline from existing outputs.

Stages covered:
  - input phenotype data
  - phenotype-retained samples
  - clinical/final covariates
  - input genotype / post-QC genotype checkpoints
  - final GWAS samples used by REGENIE step 2

Outputs:
  1. cohort_gwas_sample_counts.tsv
  2. cohort_gwas_sample_count_notes.md
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


CANONICAL_COHORTS = [
    "ROSMAP",
    "Mayo",
    "MSBB",
    "CMC_MSSM",
    "CMC_PENN",
    "CMC_PITT",
    "GTEx",
    "NABEC",
    "NIMH_HBCC_1M",
    "NIMH_HBCC_h650",
    "NIMH_HBCC_Omni5M",
    "GVEX",
]

ALIASES = {
    "CMC_MSSM": ["CMC_MSSM", "MSSM"],
    "GTEx": ["GTEx", "gtex"],
}

RE_FINAL_LOG_LOADED = re.compile(r"^\s*(\d+)\s+samples .* loaded from", re.MULTILINE)
RE_FINAL_LOG_KEEP = re.compile(r"--keep:\s+(\d+)\s+samples remaining\.", re.MULTILINE)
RE_REGENIE_PSAM = re.compile(r"psam\s*:\s*\[[^\]]+\]\s+n_samples\s*=\s*(\d+)")
RE_REGENIE_PHENO = re.compile(r"phenotyped individuals with no missing data\s*=\s*(\d+)")
RE_REGENIE_COVAR = re.compile(r"individuals with covariate data\s*=\s*(\d+)")
RE_REGENIE_USED = re.compile(r"number of individuals used in analysis\s*=\s*(\d+)")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--results-dir",
        default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results",
        help="Results directory containing one subdirectory per cohort",
    )
    p.add_argument(
        "--output-dir",
        default="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/gwas_sample_counts",
        help="Directory for the generated reports",
    )
    return p.parse_args()


def count_delimited_rows(path: Path | None, delimiter: str, has_header: bool = True) -> int | None:
    if path is None or not path.exists():
        return None
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        rows = sum(1 for row in reader if row)
    if has_header:
        rows -= 1
    return max(rows, 0)


def count_whitespace_rows(path: Path | None, has_header: bool = False) -> int | None:
    if path is None or not path.exists():
        return None
    rows = 0
    with path.open() as handle:
        for line in handle:
            if line.strip():
                rows += 1
    if has_header:
        rows -= 1
    return max(rows, 0)


def read_text(path: Path) -> str | None:
    if not path.exists():
        return None
    return path.read_text(errors="ignore")


def first_existing(base: Path, candidates: list[str]) -> Path | None:
    for name in candidates:
        path = base / name
        if path.exists():
            return path
    return None


def resolve_cohort_dir(results_dir: Path, cohort: str) -> tuple[Path | None, list[str]]:
    tried = ALIASES.get(cohort, [cohort])
    for name in tried:
        path = results_dir / name
        if path.exists():
            return path, tried
    return None, tried


def find_named_file(cohort_dir: Path, cohort: str, suffix: str) -> Path | None:
    names = [cohort]
    if cohort == "CMC_MSSM":
        names.append("MSSM")
    if cohort.startswith("NIMH_HBCC_"):
        names.append("CMC_HBCC")
    if cohort == "GTEx":
        names.append("gtex")
    for stem in names:
        candidate = cohort_dir / f"{stem}{suffix}"
        if candidate.exists():
            return candidate
    hits = list(cohort_dir.glob(f"*{suffix}"))
    return hits[0] if hits else None


def parse_final_log(path: Path | None) -> tuple[int | None, int | None]:
    if path is None:
        return None, None
    text = read_text(path)
    if not text:
        return None, None
    loaded = RE_FINAL_LOG_LOADED.search(text)
    kept = RE_FINAL_LOG_KEEP.search(text)
    return (
        int(loaded.group(1)) if loaded else None,
        int(kept.group(1)) if kept else None,
    )


def parse_regenie_step2_logs(regenie_dir: Path) -> dict[str, int | str | None]:
    result = {
        "regenie_step2_log_count": 0,
        "regenie_step2_psam_n": None,
        "regenie_step2_phenotyped_n": None,
        "regenie_step2_covariate_n": None,
        "final_gwas_n": None,
        "regenie_step2_unique_final_n": None,
    }
    if not regenie_dir.exists():
        return result

    psam_vals = []
    pheno_vals = []
    covar_vals = []
    used_vals = []

    for log in sorted(regenie_dir.glob("*.log")):
        text = read_text(log)
        if not text:
            continue
        result["regenie_step2_log_count"] += 1
        if m := RE_REGENIE_PSAM.search(text):
            psam_vals.append(int(m.group(1)))
        if m := RE_REGENIE_PHENO.search(text):
            pheno_vals.append(int(m.group(1)))
        if m := RE_REGENIE_COVAR.search(text):
            covar_vals.append(int(m.group(1)))
        if m := RE_REGENIE_USED.search(text):
            used_vals.append(int(m.group(1)))

    def stable_value(values: list[int]) -> int | None:
        return sorted(set(values))[0] if values else None

    result["regenie_step2_psam_n"] = stable_value(psam_vals)
    result["regenie_step2_phenotyped_n"] = stable_value(pheno_vals)
    result["regenie_step2_covariate_n"] = stable_value(covar_vals)
    result["final_gwas_n"] = stable_value(used_vals)
    if used_vals:
        result["regenie_step2_unique_final_n"] = ",".join(str(v) for v in sorted(set(used_vals)))

    return result


def collect_row_counts(cohort_dir: Path, cohort: str) -> dict[str, int | None]:
    counts = {}
    counts["input_metadata_n"] = count_delimited_rows(cohort_dir / "metadata_cleaned.csv", ",", True)
    counts["input_cell_proportions_n"] = count_delimited_rows(cohort_dir / "cell_proportions.csv", ",", True)
    counts["input_cell_proportions_scaled_n"] = count_delimited_rows(
        cohort_dir / "cell_proportions_scaled.csv", ",", True
    )
    counts["phenotype_retained_n"] = count_whitespace_rows(cohort_dir / "samples_with_phenotypes.txt", has_header=True)
    counts["phenotype_table_n"] = count_delimited_rows(cohort_dir / "phenotypes_RINT.txt", "\t", True)
    counts["clinical_covariates_n"] = count_delimited_rows(cohort_dir / "clinical_covariates.txt", "\t", True)
    counts["final_covariates_n"] = count_delimited_rows(cohort_dir / "covariates.txt", "\t", True)
    counts["post_standard_qc_genotype_n"] = count_whitespace_rows(
        find_named_file(cohort_dir, cohort, ".QC.fam"), has_header=False
    )
    counts["post_het_qc_genotype_n"] = count_whitespace_rows(
        find_named_file(cohort_dir, cohort, ".QC.valid"), has_header=False
    )
    counts["final_genotype_n"] = count_delimited_rows(
        find_named_file(cohort_dir, cohort, ".QC.final.psam"), "\t", True
    )
    counts["pca_n"] = count_delimited_rows(cohort_dir / "pca.csv", ",", True)
    return counts


def build_notes(row: dict[str, object]) -> str:
    notes = []
    if row["cohort_dir_name"] != row["cohort"]:
        notes.append(f"used directory alias {row['cohort_dir_name']}")
    if row["input_metadata_n"] is None:
        notes.append("missing metadata_cleaned.csv")
    if row["phenotype_retained_n"] is None:
        notes.append("missing samples_with_phenotypes.txt")
    if row["post_standard_qc_genotype_n"] is None:
        notes.append("missing *.QC.fam")
    if row["post_het_qc_genotype_n"] is None:
        notes.append("missing *.QC.valid")
    if row["final_genotype_n"] is None and row["regenie_step2_psam_n"] is not None:
        notes.append("final_genotype_n filled from REGENIE step2 psam_n")
    if row["input_genotype_n"] is None and row["regenie_step2_psam_n"] is not None:
        notes.append("input_genotype_n unavailable from QC log")
    if row.get("qc_final_log_inconsistent"):
        notes.append("QC.final.log sample counts inconsistent with cohort-specific QC/REGENIE outputs")
    if row.get("final_psam_inconsistent"):
        notes.append("QC.final.psam count inconsistent with REGENIE step2 psam_n; used REGENIE value")
    if row.get("pca_inconsistent"):
        notes.append("pca.csv count inconsistent with REGENIE step2 psam_n")
    if row.get("covariates_inconsistent"):
        notes.append("covariates.txt count inconsistent with REGENIE step2 covariate_n")
    if row["regenie_step2_log_count"] != 19:
        notes.append(f"expected 19 step2 logs, found {row['regenie_step2_log_count']}")
    if row["regenie_step2_unique_final_n"] and "," in str(row["regenie_step2_unique_final_n"]):
        notes.append(f"step2 final N varies by cell type: {row['regenie_step2_unique_final_n']}")
    return "; ".join(notes) if notes else "complete"


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for cohort in CANONICAL_COHORTS:
        cohort_dir, aliases = resolve_cohort_dir(results_dir, cohort)
        if cohort_dir is None:
            rows.append({
                "cohort": cohort,
                "cohort_dir_name": "MISSING",
                "input_metadata_n": None,
                "input_cell_proportions_n": None,
                "input_cell_proportions_scaled_n": None,
                "phenotype_retained_n": None,
                "phenotype_table_n": None,
                "clinical_covariates_n": None,
                "final_covariates_n": None,
                "input_genotype_n": None,
                "post_standard_qc_genotype_n": None,
                "post_het_qc_genotype_n": None,
                "final_genotype_n": None,
                "pca_n": None,
                "regenie_step2_psam_n": None,
                "regenie_step2_phenotyped_n": None,
                "regenie_step2_covariate_n": None,
                "final_gwas_n": None,
                "regenie_step2_log_count": 0,
                "regenie_step2_unique_final_n": None,
                "count_source_notes": f"missing cohort directory; tried aliases: {', '.join(aliases)}",
            })
            continue

        counts = collect_row_counts(cohort_dir, cohort)
        final_log = find_named_file(cohort_dir, cohort, ".QC.final.log")
        input_genotype_n, final_log_kept_n = parse_final_log(final_log)
        regenie = parse_regenie_step2_logs(cohort_dir / "regenie_step2")

        raw_final_genotype_n = counts["final_genotype_n"]
        final_genotype_n = raw_final_genotype_n
        final_psam_inconsistent = False
        if (
            raw_final_genotype_n is not None
            and regenie["regenie_step2_psam_n"] is not None
            and abs(raw_final_genotype_n - regenie["regenie_step2_psam_n"]) > 5
        ):
            final_genotype_n = regenie["regenie_step2_psam_n"]
            final_psam_inconsistent = True
        elif final_genotype_n is None:
            final_genotype_n = regenie["regenie_step2_psam_n"]

        qc_final_log_inconsistent = False
        if (
            final_log_kept_n is not None
            and counts["post_het_qc_genotype_n"] is not None
            and abs(final_log_kept_n - counts["post_het_qc_genotype_n"]) > 5
        ):
            qc_final_log_inconsistent = True
            input_genotype_n = None

        pca_inconsistent = False
        if (
            counts["pca_n"] is not None
            and regenie["regenie_step2_psam_n"] is not None
            and abs(counts["pca_n"] - regenie["regenie_step2_psam_n"]) > 5
        ):
            pca_inconsistent = True

        covariates_inconsistent = False
        if (
            counts["final_covariates_n"] is not None
            and regenie["regenie_step2_covariate_n"] is not None
            and abs(counts["final_covariates_n"] - regenie["regenie_step2_covariate_n"]) > 5
        ):
            covariates_inconsistent = True

        row = {
            "cohort": cohort,
            "cohort_dir_name": cohort_dir.name,
            "input_metadata_n": counts["input_metadata_n"],
            "input_cell_proportions_n": counts["input_cell_proportions_n"],
            "input_cell_proportions_scaled_n": counts["input_cell_proportions_scaled_n"],
            "phenotype_retained_n": counts["phenotype_retained_n"],
            "phenotype_table_n": counts["phenotype_table_n"],
            "clinical_covariates_n": counts["clinical_covariates_n"],
            "final_covariates_n": counts["final_covariates_n"],
            "input_genotype_n": input_genotype_n,
            "post_standard_qc_genotype_n": counts["post_standard_qc_genotype_n"],
            "post_het_qc_genotype_n": counts["post_het_qc_genotype_n"],
            "final_genotype_n": final_genotype_n,
            "pca_n": counts["pca_n"],
            "regenie_step2_psam_n": regenie["regenie_step2_psam_n"],
            "regenie_step2_phenotyped_n": regenie["regenie_step2_phenotyped_n"],
            "regenie_step2_covariate_n": regenie["regenie_step2_covariate_n"],
            "final_gwas_n": regenie["final_gwas_n"],
            "regenie_step2_log_count": regenie["regenie_step2_log_count"],
            "regenie_step2_unique_final_n": regenie["regenie_step2_unique_final_n"],
            "final_qc_log_keep_n": final_log_kept_n,
            "qc_final_log_inconsistent": qc_final_log_inconsistent,
            "final_psam_inconsistent": final_psam_inconsistent,
            "pca_inconsistent": pca_inconsistent,
            "covariates_inconsistent": covariates_inconsistent,
        }
        row["count_source_notes"] = build_notes(row)
        rows.append(row)

    fieldnames = [
        "cohort",
        "cohort_dir_name",
        "input_metadata_n",
        "input_cell_proportions_n",
        "input_cell_proportions_scaled_n",
        "phenotype_retained_n",
        "phenotype_table_n",
        "clinical_covariates_n",
        "final_covariates_n",
        "input_genotype_n",
        "post_standard_qc_genotype_n",
        "post_het_qc_genotype_n",
        "final_genotype_n",
        "pca_n",
        "regenie_step2_psam_n",
        "regenie_step2_phenotyped_n",
        "regenie_step2_covariate_n",
        "final_gwas_n",
        "regenie_step2_log_count",
        "regenie_step2_unique_final_n",
        "final_qc_log_keep_n",
        "count_source_notes",
    ]

    out_tsv = output_dir / "cohort_gwas_sample_counts.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    out_md = output_dir / "cohort_gwas_sample_count_notes.md"
    with out_md.open("w") as handle:
        handle.write("# GWAS Sample Count Provenance\n\n")
        handle.write("## Stage Definitions\n")
        handle.write("- `input_metadata_n`: rows in `metadata_cleaned.csv`\n")
        handle.write("- `input_cell_proportions_n`: rows in `cell_proportions.csv`\n")
        handle.write("- `phenotype_retained_n`: rows in `samples_with_phenotypes.txt` (excluding header)\n")
        handle.write("- `clinical_covariates_n`: rows in `clinical_covariates.txt`\n")
        handle.write("- `final_covariates_n`: rows in `covariates.txt`\n")
        handle.write("- `input_genotype_n`: `samples loaded` parsed from `*.QC.final.log`\n")
        handle.write("- `post_standard_qc_genotype_n`: rows in `*.QC.fam`\n")
        handle.write("- `post_het_qc_genotype_n`: rows in `*.QC.valid`\n")
        handle.write("- `final_genotype_n`: rows in `*.QC.final.psam`, with fallback to `regenie_step2 psam n_samples`\n")
        handle.write("- `final_gwas_n`: `number of individuals used in analysis` from REGENIE step 2 logs\n\n")
        handle.write("## Cohort Notes\n")
        for row in rows:
            handle.write(f"- `{row['cohort']}`: {row['count_source_notes']}\n")

    print(f"Wrote sample count table: {out_tsv}")
    print(f"Wrote provenance notes:   {out_md}")


if __name__ == "__main__":
    main()
