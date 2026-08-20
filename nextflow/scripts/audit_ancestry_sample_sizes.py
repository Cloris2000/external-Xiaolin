#!/usr/bin/env python3
"""
Audit individual-level ancestry labels for the 15-cohort meta set, intersect with
GWAS-eligible samples (covariates.txt IIDs), write N tables, and build FID/IID
keep-lists for (cohort, ancestry) pairs with N >= min_n.

Hispanic-first mapping: isHispanic/HISP -> AMR; then race/ethnicity -> EUR/AFR/EAS.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
OUT_DIR = PROJECT / "docs" / "ancestry_specific"
KEEP_DIR = OUT_DIR / "keep_lists"

MIN_N_DEFAULT = 50

COHORTS_15 = [
    "ROSMAP",
    "ROSMAP_array",
    "Mayo",
    "MSBB",
    "CMC_MSSM",
    "CMC_PENN",
    "CMC_PITT",
    "GTEx_v10",
    "NABEC",
    "NIMH_HBCC_1M",
    "NIMH_HBCC_h650",
    "NIMH_HBCC_Omni5M",
    "GVEX",
    "AMP_AD_Rush",
    "AMP_AD_Mayo",
]

# Harmonized ancestry codes
EUR, AFR, AMR, EAS, OTHER, UNK = "EUR", "AFR", "AMR", "EAS", "OTHER", "UNKNOWN"

EUR_LABELS = {
    "cauc",
    "white",
    "w",
    "european",
    "us caucasian",
    "1",
    "1.0",
    "nhw",
    "non-hispanic white",
}
AFR_LABELS = {
    "aa",
    "black",
    "b",
    "black or african american",
    "african american",
    "african background",
    "2",
    "2.0",
}
AMR_LABELS = {
    "hisp",
    "h",
    "hispanic",
    "latino",
    "latina",
    "latin american",
    "hispanic or latino",
}
EAS_LABELS = {"as", "a", "asian", "east asian", "3", "3.0"}


def _norm(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().lower()


def harmonize_row(
    race=None,
    ethnicity=None,
    is_hispanic=None,
    *,
    hispanic_first: bool = True,
) -> str:
    """Map self-report fields to EUR/AFR/AMR/EAS/OTHER/UNKNOWN."""
    # Explicit Hispanic / Latino
    hisp_vals = {_norm(is_hispanic), _norm(ethnicity), _norm(race)}
    if hispanic_first:
        if _norm(is_hispanic) in {"true", "1", "1.0", "yes", "y"}:
            return AMR
        if _norm(ethnicity) in AMR_LABELS or _norm(race) in AMR_LABELS:
            return AMR

    for val in (_norm(ethnicity), _norm(race)):
        if not val:
            continue
        if val in AMR_LABELS:
            return AMR
        if val in EUR_LABELS:
            return EUR
        if val in AFR_LABELS:
            return AFR
        if val in EAS_LABELS:
            return EAS
        if val in {"multiracial", "mixed", "non-white", "other", "7", "7.0", "u", "unknown"}:
            return OTHER
    return UNK


def load_gwas_ids(cohort: str) -> pd.DataFrame:
    """Return FID/IID table for samples that entered GWAS (covariates.txt)."""
    cov = RESULTS / cohort / "covariates.txt"
    if not cov.exists():
        raise FileNotFoundError(f"Missing covariates for {cohort}: {cov}")
    df = pd.read_csv(cov, sep="\t")
    # Standardize
    if "FID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "FID"})
    if "IID" not in df.columns:
        df = df.rename(columns={df.columns[1]: "IID"})
    out = df[["FID", "IID"]].copy()
    out["FID"] = out["FID"].astype(str)
    out["IID"] = out["IID"].astype(str)
    return out.drop_duplicates("IID")


def _dedupe_map(df: pd.DataFrame, id_col: str, anc_col: str = "ancestry") -> Dict[str, str]:
    sub = df.dropna(subset=[id_col]).copy()
    sub[id_col] = sub[id_col].astype(str)
    # Prefer non-UNKNOWN if duplicates
    sub["_rank"] = sub[anc_col].map(
        {EUR: 0, AFR: 1, AMR: 2, EAS: 3, OTHER: 4, UNK: 5}
    ).fillna(9)
    sub = sub.sort_values("_rank").drop_duplicates(id_col, keep="first")
    return dict(zip(sub[id_col], sub[anc_col]))


def ancestry_for_rosmap(cohort: str, gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    bio = pd.read_csv(
        "/nethome/kcni/xzhou/GWAS_tut/AMP-AD/ROSMAP_biospecimen_metadata.csv",
        low_memory=False,
    )
    meta = pd.read_csv(
        "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_combined_metrics.csv",
        low_memory=False,
    )
    race_by_ind = (
        meta.dropna(subset=["individualID"])
        .assign(
            ancestry=lambda d: d["race"].map(
                lambda r: harmonize_row(race=r)
            )
        )[["individualID", "ancestry"]]
    )
    # Array: MAP + zero-padded projid
    if cohort == "ROSMAP_array":
        meta2 = meta.dropna(subset=["projid"]).copy()
        meta2["map_id"] = (
            "MAP"
            + meta2["projid"].astype(float).astype(int).astype(str).str.zfill(8)
        )
        meta2["ancestry"] = meta2["race"].map(lambda r: harmonize_row(race=r))
        amap = _dedupe_map(meta2, "map_id")
        anc = gwas["IID"].map(amap)
        src = {
            "label_source": "ROSMAP_combined_metrics.race via projid->MAP*",
            "id_join": "IID == MAP + zfill8(projid)",
            "assignment_mode": "metadata_race",
        }
        return anc, src

    link = bio[["specimenID", "individualID"]].dropna()
    link["specimenID"] = link["specimenID"].astype(str)
    link = link.merge(race_by_ind, on="individualID", how="left")
    amap = _dedupe_map(link, "specimenID")
    anc = gwas["IID"].map(amap)
    src = {
        "label_source": "ROSMAP_combined_metrics.race via biospecimen specimenID",
        "id_join": "covariates.IID == biospecimen.specimenID -> individualID",
        "assignment_mode": "metadata_race",
    }
    return anc, src


def ancestry_for_mayo(gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    meta = pd.read_csv(
        "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/Mayo_meta_tissue_counts_Nov_18.csv",
        low_memory=False,
    )
    meta = meta.dropna(subset=["individualID"]).copy()
    meta["ancestry"] = meta["race"].map(lambda r: harmonize_row(race=r))
    amap = _dedupe_map(meta, "individualID")
    return gwas["IID"].map(amap), {
        "label_source": "Mayo_meta_tissue_counts_Nov_18.race",
        "id_join": "covariates.IID == individualID",
        "assignment_mode": "metadata_race",
    }


def ancestry_for_msbb(gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    bio = pd.read_csv(
        "/nethome/kcni/xzhou/GWAS_tut/AMP-AD/MSBB_biospecimen_metadata.csv",
        low_memory=False,
    )
    meta = pd.read_csv(
        "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/msbb_meta.csv",
        low_memory=False,
    )
    race = meta.dropna(subset=["individualID"]).copy()
    race["ancestry"] = [
        harmonize_row(race=r, ethnicity=e)
        for r, e in zip(race["race"], race["ethnicity"])
    ]
    link = bio[["specimenID", "individualID"]].dropna()
    link["specimenID"] = link["specimenID"].astype(str)
    link = link.merge(race[["individualID", "ancestry"]], on="individualID", how="left")
    amap = _dedupe_map(link, "specimenID")
    return gwas["IID"].map(amap), {
        "label_source": "msbb_meta.race/ethnicity via MSBB biospecimen",
        "id_join": "covariates.IID == biospecimen.specimenID (WGS numeric)",
        "assignment_mode": "metadata_race_ethnicity",
    }


def ancestry_for_cmc(cohort: str, gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    snp = pd.read_csv(
        "/external/rprshnas01/netdata_kcni/stlab/CMC_genotypes/SNPs/Release3/Metadata/CMC_Human_SNP_metadata.csv",
        low_memory=False,
    )
    snp["geno_fid"] = "0_" + snp["Genotyping_Sample_ID"].astype(str)
    mc = pd.read_csv(RESULTS / cohort / "metadata_cleaned.csv", low_memory=False)
    eth_col = "ethnicity" if "ethnicity" in mc.columns else None
    race_col = "race" if "race" in mc.columns else None
    mc = mc.dropna(subset=["individualID"]).copy()
    mc["ancestry"] = [
        harmonize_row(
            race=mc.loc[i, race_col] if race_col else None,
            ethnicity=mc.loc[i, eth_col] if eth_col else None,
        )
        for i in mc.index
    ]
    link = snp[["geno_fid", "Individual_ID"]].rename(
        columns={"Individual_ID": "individualID"}
    )
    link = link.merge(mc[["individualID", "ancestry"]], on="individualID", how="left")
    amap = _dedupe_map(link, "geno_fid")
    return gwas["IID"].map(amap), {
        "label_source": f"{cohort}/metadata_cleaned ethnicity/race via CMC SNP metadata",
        "id_join": "covariates.IID == 0_+Genotyping_Sample_ID -> Individual_ID",
        "assignment_mode": "metadata_ethnicity",
    }


def ancestry_for_hbcc(cohort: str, gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    snp = pd.read_csv(
        "/external/rprshnas01/netdata_kcni/stlab/CMC_genotypes/SNPs/Release3/Metadata/CMC_Human_SNP_metadata.csv",
        low_memory=False,
    )
    snp["geno_fid"] = "0_" + snp["Genotyping_Sample_ID"].astype(str)
    hb = pd.read_csv(PROJECT / "data_input/nimh_hbcc/HBCC_metadata.csv", low_memory=False)
    hb = hb.dropna(subset=["individualID"]).copy()
    hb["ancestry"] = hb["ethnicity"].map(lambda e: harmonize_row(ethnicity=e))
    link = snp[["geno_fid", "Individual_ID"]].rename(
        columns={"Individual_ID": "individualID"}
    )
    link = link.merge(hb[["individualID", "ancestry"]], on="individualID", how="left")
    amap = _dedupe_map(link, "geno_fid")
    return gwas["IID"].map(amap), {
        "label_source": "HBCC_metadata.ethnicity via CMC SNP metadata",
        "id_join": "covariates.IID == 0_+Genotyping_Sample_ID -> Individual_ID",
        "assignment_mode": "metadata_ethnicity",
    }


def ancestry_for_gvex(gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    mc = pd.read_csv(RESULTS / "GVEX" / "metadata_cleaned.csv", low_memory=False)
    cap = pd.read_csv(
        "/external/rprshnas01/external_data/psychencode/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv",
        low_memory=False,
    )
    cap = cap.dropna(subset=["individualID"]).copy()
    cap["ancestry"] = [
        harmonize_row(race=r, ethnicity=e)
        for r, e in zip(cap.get("race", pd.Series(index=cap.index)), cap["ethnicity"])
    ]
    # GWAS IID == specimen or individual in metadata_cleaned
    link_spec = mc[["specimenID", "individualID"]].dropna()
    link_spec["specimenID"] = link_spec["specimenID"].astype(str)
    link_spec = link_spec.merge(
        cap[["individualID", "ancestry"]], on="individualID", how="left"
    )
    amap = _dedupe_map(link_spec, "specimenID")
    # also allow IID == individualID
    amap_ind = _dedupe_map(
        cap[["individualID", "ancestry"]].rename(columns={"individualID": "id"}),
        "id",
    )
    anc = gwas["IID"].map(amap)
    missing = anc.isna()
    anc = anc.where(~missing, gwas["IID"].map(amap_ind))
    return anc, {
        "label_source": "CapstoneCollection_Metadata_Clinical.ethnicity/race",
        "id_join": "covariates.IID == GVEX specimenID/individualID",
        "assignment_mode": "metadata_ethnicity",
    }


def ancestry_for_amp_ad(cohort: str, gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    meta_name = {
        "AMP_AD_Rush": "Rush_DLPFC_metadata.csv",
        "AMP_AD_Mayo": "Mayo_DLPFC_metadata.csv",
    }[cohort]
    meta = pd.read_csv(
        PROJECT / "data_input/amp_ad_diverse" / meta_name, low_memory=False
    )
    meta = meta.copy()
    meta["ancestry"] = [
        harmonize_row(race=r, is_hispanic=h)
        for r, h in zip(meta["race"], meta["isHispanic"])
    ]
    # Rush: IID == individualID; Mayo: IID like 11387_DLPFC_WGS -> specimen 11387_DLPFC
    if cohort == "AMP_AD_Rush":
        amap = _dedupe_map(meta.assign(id=meta["individualID"].astype(str)), "id")
        anc = gwas["IID"].map(amap)
        join = "covariates.IID == individualID"
    else:
        meta["specimen_key"] = meta["specimenID"].astype(str)
        amap = _dedupe_map(meta.rename(columns={"specimen_key": "id"}), "id")
        # try raw and stripped _WGS
        anc = gwas["IID"].map(amap)
        stripped = gwas["IID"].str.replace(r"_WGS$", "", regex=True)
        anc = anc.fillna(stripped.map(amap))
        join = "covariates.IID (strip _WGS) == specimenID"
    return anc, {
        "label_source": f"{meta_name} race + isHispanic (Hispanic-first)",
        "id_join": join,
        "assignment_mode": "metadata_race_hispanic",
    }


def ancestry_for_nabec(gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    meta = pd.read_csv(
        "/nethome/kcni/xzhou/GWAS_tut/NABEC/NABEC_metadata_combined.csv",
        low_memory=False,
    )
    meta = meta.dropna(subset=["individualID"]).copy()
    meta["ancestry"] = meta["RACE"].map(lambda r: harmonize_row(race=r))
    amap = _dedupe_map(meta, "individualID")
    return gwas["IID"].map(amap), {
        "label_source": "NABEC_metadata_combined.RACE",
        "id_join": "covariates.IID == individualID",
        "assignment_mode": "metadata_race",
    }


def ancestry_for_gtex(gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    # No ancestry column in BA9 sample metadata used by the pipeline.
    anc = pd.Series([EUR] * len(gwas), index=gwas.index)
    return anc, {
        "label_source": "none (GTEx_v10 BA9 metadata has no race/ethnicity column)",
        "id_join": "cohort_policy_EUR for all GWAS samples",
        "assignment_mode": "cohort_policy_EUR",
    }


def assign_ancestry(cohort: str, gwas: pd.DataFrame) -> Tuple[pd.Series, dict]:
    if cohort in ("ROSMAP", "ROSMAP_array"):
        return ancestry_for_rosmap(cohort, gwas)
    if cohort == "Mayo":
        return ancestry_for_mayo(gwas)
    if cohort == "MSBB":
        return ancestry_for_msbb(gwas)
    if cohort.startswith("CMC_"):
        return ancestry_for_cmc(cohort, gwas)
    if cohort.startswith("NIMH_HBCC"):
        return ancestry_for_hbcc(cohort, gwas)
    if cohort == "GVEX":
        return ancestry_for_gvex(gwas)
    if cohort.startswith("AMP_AD_"):
        return ancestry_for_amp_ad(cohort, gwas)
    if cohort == "NABEC":
        return ancestry_for_nabec(gwas)
    if cohort == "GTEx_v10":
        return ancestry_for_gtex(gwas)
    raise ValueError(f"No ancestry mapper for {cohort}")


def write_keep_list(path: Path, gwas: pd.DataFrame, mask: pd.Series) -> int:
    sub = gwas.loc[mask, ["FID", "IID"]].drop_duplicates()
    path.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(path, sep="\t", index=False)
    return len(sub)


def main(min_n: int = MIN_N_DEFAULT) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    KEEP_DIR.mkdir(parents=True, exist_ok=True)

    source_rows = []
    count_rows = []
    sample_rows = []
    pair_rows = []

    for cohort in COHORTS_15:
        print(f"[audit] {cohort}")
        gwas = load_gwas_ids(cohort)
        anc, src = assign_ancestry(cohort, gwas)
        anc = anc.fillna(UNK)
        gwas = gwas.copy()
        gwas["ancestry"] = anc.values

        source_rows.append(
            {
                "cohort": cohort,
                "n_gwas_samples": len(gwas),
                "n_labeled": int((gwas["ancestry"] != UNK).sum()),
                "n_unlabeled": int((gwas["ancestry"] == UNK).sum()),
                **src,
            }
        )

        for ancestry, n in gwas["ancestry"].value_counts().items():
            count_rows.append(
                {
                    "cohort": cohort,
                    "ancestry": ancestry,
                    "n": int(n),
                    "n_gwas_total": len(gwas),
                    "pct": round(100.0 * n / len(gwas), 2),
                    "passes_min_n": int(n) >= min_n,
                }
            )

        # Per-sample table (for reproducibility)
        for _, row in gwas.iterrows():
            sample_rows.append(
                {
                    "cohort": cohort,
                    "FID": row["FID"],
                    "IID": row["IID"],
                    "ancestry": row["ancestry"],
                }
            )

        for ancestry in [EUR, AFR, AMR, EAS]:
            mask = gwas["ancestry"] == ancestry
            n = int(mask.sum())
            if n < min_n:
                continue
            keep_path = KEEP_DIR / f"{cohort}_{ancestry}_samples_to_keep.txt"
            n_written = write_keep_list(keep_path, gwas, mask)
            pair_rows.append(
                {
                    "cohort": cohort,
                    "ancestry": ancestry,
                    "n": n_written,
                    "study_name": f"{cohort}_{ancestry}",
                    "keep_list": str(keep_path.relative_to(PROJECT)),
                    "run_gwas": 1,
                }
            )

    sources = pd.DataFrame(source_rows)
    counts = pd.DataFrame(count_rows)
    pairs = pd.DataFrame(pair_rows)
    samples = pd.DataFrame(sample_rows)

    sources.to_csv(OUT_DIR / "ancestry_label_sources.tsv", sep="\t", index=False)
    counts.to_csv(OUT_DIR / "ancestry_sample_counts.tsv", sep="\t", index=False)
    pairs.to_csv(OUT_DIR / f"ancestry_pairs_n{min_n}.tsv", sep="\t", index=False)
    samples.to_csv(OUT_DIR / "ancestry_sample_assignments.tsv", sep="\t", index=False)

    # Meta eligibility: ancestries with >=2 cohorts passing N gate
    meta_rows = []
    if len(pairs):
        for ancestry, sub in pairs.groupby("ancestry"):
            meta_rows.append(
                {
                    "ancestry": ancestry,
                    "n_cohorts": len(sub),
                    "total_n": int(sub["n"].sum()),
                    "cohorts": ",".join(sub["cohort"].tolist()),
                    "run_meta": int(len(sub) >= 2),
                    "note": "meta" if len(sub) >= 2 else "single_cohort_only",
                }
            )
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUT_DIR / f"ancestry_meta_eligibility_n{min_n}.tsv", sep="\t", index=False)

    # README
    readme = OUT_DIR / "README.md"
    readme.write_text(
        f"""# Ancestry-specific GWAS sample audit

## Rules
- Ancestry from **individual-level metadata** (self-report race/ethnicity/isHispanic).
- **Hispanic-first:** `isHispanic`/HISP → `AMR`; else map race/ethnicity → EUR/AFR/EAS.
- GWAS-eligible N = samples in `results/<COHORT>/covariates.txt`.
- Subset GWAS only if **N ≥ {min_n}**.
- Ancestry meta only if **≥ 2** cohorts pass the N gate for that ancestry.
- GTEx_v10 has no race column in the pipeline metadata → `cohort_policy_EUR`.

## Outputs
- `ancestry_label_sources.tsv` — per-cohort label source and ID join
- `ancestry_sample_counts.tsv` — cohort × ancestry counts
- `ancestry_pairs_n{min_n}.tsv` — runnable (cohort, ancestry) GWAS pairs
- `ancestry_meta_eligibility_n{min_n}.tsv` — which ancestry metas have ≥2 cohorts
- `keep_lists/<COHORT>_<ANCESTRY>_samples_to_keep.txt` — FID/IID keep lists
- `ancestry_sample_assignments.tsv` — per-sample assignments

Generated by `scripts/audit_ancestry_sample_sizes.py`.
""",
        encoding="utf-8",
    )

    print("\n=== Label sources ===")
    print(sources.to_string(index=False))
    print("\n=== Counts (wide) ===")
    wide = counts.pivot_table(
        index="cohort", columns="ancestry", values="n", fill_value=0, aggfunc="sum"
    )
    print(wide.to_string())
    print(f"\n=== Pairs N>={min_n} ({len(pairs)}) ===")
    print(pairs.to_string(index=False) if len(pairs) else "(none)")
    print("\n=== Meta eligibility ===")
    print(meta.to_string(index=False) if len(meta) else "(none)")
    print(f"\nWrote outputs under {OUT_DIR}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=MIN_N_DEFAULT)
    args = ap.parse_args()
    main(min_n=args.min_n)
