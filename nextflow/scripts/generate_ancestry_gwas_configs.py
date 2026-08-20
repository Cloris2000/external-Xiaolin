#!/usr/bin/env python3
"""
Generate ancestry-subset GWAS Nextflow configs and a SLURM launch list.

Reads docs/ancestry_specific/reuse_vs_regwas.tsv:
  - re_gwas       -> ancestry_configs/nextflow.config.combined.<study>
  - reuse_existing -> documented for meta only (no new GWAS)

Overlays includeConfig the parent cohort config and override study / output_dir /
samples_to_keep / work_dir / cohort_configs. Uses precomputed cell proportions from
the parent results dir to skip RNA/deconv; genotyping QC + REGENIE still re-run.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
OUT_CFG = PROJECT / "ancestry_configs"
DOCS = PROJECT / "docs" / "ancestry_specific"

# study name (COHORT) -> parent nextflow.config.combined.* basename
PARENT_CONFIG = {
    "ROSMAP": "rosmap",
    "ROSMAP_array": "rosmap_array",
    "Mayo": "mayo",
    "MSBB": "msbb",
    "CMC_MSSM": "cmc_mssm",
    "CMC_PENN": "cmc_penn",
    "CMC_PITT": "cmc_pitt",
    "GTEx_v10": "gtex_v10",
    "NABEC": "nabec",
    "NIMH_HBCC_1M": "nimh_hbcc_1M",
    "NIMH_HBCC_h650": "nimh_hbcc_h650",
    "NIMH_HBCC_Omni5M": "nimh_hbcc_Omni5M",
    "GVEX": "gvex",
    "AMP_AD_Rush": "amp_ad_rush",
    "AMP_AD_Mayo": "amp_ad_mayo",
}

# QC / normalized-VCF basename used by genotyping_qc.py
# ({qc_study}.QC.{{chrom}}.normalized.vcf.gz). Must match parent VCF prefixes
# when ancestry study name differs (e.g. NIMH_HBCC_1M_AFR -> CMC_HBCC).
QC_STUDY = {
    "ROSMAP_array": "ROSMAP_array",
    "MSBB": "MSBB",
    "CMC_MSSM": "MSSM",
    "CMC_PENN": "CMC",
    "CMC_PITT": "CMC",
    "NIMH_HBCC_1M": "CMC_HBCC",
    "NIMH_HBCC_h650": "CMC_HBCC",
    "NIMH_HBCC_Omni5M": "CMC_HBCC",
    # AMP-AD uses vcf_pattern (AMP_AD_Diverse.chr{{chrom}}); no study-named normalized files
}


def config_text(cohort: str, ancestry: str, n: int, keep_rel: str) -> str:
    study = f"{cohort}_{ancestry}"
    parent = PARENT_CONFIG[cohort]
    parent_results = f"${{projectDir}}/results/{cohort}"
    qc_study = QC_STUDY.get(cohort)
    if qc_study:
        geno_line = f'genotyping_study = "{qc_study}"  // VCF/QC basename (not ancestry study name)'
        pgen_base = qc_study
    else:
        # Rely on parent vcf_pattern fallback; name QC outputs after ancestry study
        geno_line = "genotyping_study = null"
        pgen_base = study

    return f"""/*
 * Ancestry-subset GWAS: {study} (N={n})
 *
 * Parent config: nextflow.config.combined.{parent}
 * Keep list: {keep_rel}
 * QC/VCF basename: {pgen_base}
 *
 * Strategy:
 *   - Reuse parent cell proportions / metadata (skip RNA + deconv)
 *   - Restrict samples via samples_to_keep
 *   - Re-run genotyping QC + REGENIE on the ancestry subset
 *
 * Run:
 *   mkdir -p run_dirs/{study} work/{study} logs/{study}
 *   nextflow run combined_pipeline_v2.nf \\
 *       -c ancestry_configs/nextflow.config.combined.{study} \\
 *       -work-dir work/{study} \\
 *       -resume
 */

includeConfig '../nextflow.config.combined.{parent}'

params {{
    cohorts = ["{study}"]
    study   = "{study}"
    {geno_line}

    samples_to_keep = "${{projectDir}}/{keep_rel}"

    precomputed_proportions = "{parent_results}/cell_proportions.csv"
    precomputed_metadata    = "{parent_results}/combined_metrics.csv"

    // Bump so cached CREATE_PGEN from failed pilots is not reused
    qc_version = "v4_ancestry"

    wgs_psam_file = "${{projectDir}}/results/{study}/{pgen_base}.QC.final.psam"
    pca_file      = "${{projectDir}}/results/{study}/pca.csv"

    work_dir   = "${{projectDir}}/work/{study}"
    log_dir    = "${{projectDir}}/logs/{study}"
    output_dir = "${{projectDir}}/results/{study}"

    skip_metal = true

    cohort_configs = [
        "{study}": [
            pgen_file:     "${{projectDir}}/results/{study}/{pgen_base}.QC.final",
            prune_in_file: "${{projectDir}}/results/{study}/{pgen_base}.QC.prune.in",
            pheno_file:    "${{projectDir}}/results/{study}/phenotypes_RINT.txt",
            covar_file:    "${{projectDir}}/results/{study}/covariates.txt"
        ]
    ]
}}
"""


def main() -> None:
    OUT_CFG.mkdir(parents=True, exist_ok=True)
    reuse = pd.read_csv(DOCS / "reuse_vs_regwas.tsv", sep="\t")
    pairs = pd.read_csv(DOCS / "ancestry_pairs_n50.tsv", sep="\t")
    pairs = pairs.merge(
        reuse[["cohort", "ancestry", "recommendation", "reason"]],
        on=["cohort", "ancestry"],
        how="left",
    )

    launch_rows = []
    for _, row in pairs.iterrows():
        cohort, ancestry = row["cohort"], row["ancestry"]
        study = row["study_name"]
        rec = row.get("recommendation", "re_gwas")
        if rec == "reuse_existing":
            launch_rows.append(
                {
                    "study": study,
                    "cohort": cohort,
                    "ancestry": ancestry,
                    "n": int(row["n"]),
                    "action": "reuse_existing",
                    "config": "",
                    "source_results": f"results/{cohort}",
                    "keep_list": row["keep_list"],
                }
            )
            continue

        keep_rel = row["keep_list"]
        cfg_name = f"nextflow.config.combined.{study}"
        cfg_path = OUT_CFG / cfg_name
        cfg_path.write_text(
            config_text(cohort, ancestry, int(row["n"]), keep_rel), encoding="utf-8"
        )
        launch_rows.append(
            {
                "study": study,
                "cohort": cohort,
                "ancestry": ancestry,
                "n": int(row["n"]),
                "action": "re_gwas",
                "config": f"ancestry_configs/{cfg_name}",
                "source_results": f"results/{study}",
                "keep_list": keep_rel,
            }
        )
        print(f"wrote {cfg_path.relative_to(PROJECT)}")

    launch = pd.DataFrame(launch_rows)
    launch_path = DOCS / "ancestry_gwas_launch_list.tsv"
    launch.to_csv(launch_path, sep="\t", index=False)

    # Array task list: only re_gwas rows
    regwas = launch[launch["action"] == "re_gwas"].reset_index(drop=True)
    tasks_path = DOCS / "ancestry_gwas_array_tasks.tsv"
    regwas.to_csv(tasks_path, sep="\t", index=False)

    print(f"\nLaunch list: {launch_path}")
    print(f"Array tasks ({len(regwas)} re_gwas): {tasks_path}")
    print(regwas[["study", "n", "config"]].to_string(index=False))
    print(f"\nreuse_existing ({(launch.action=='reuse_existing').sum()}):")
    print(
        launch.loc[launch.action == "reuse_existing", ["study", "n", "source_results"]].to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
