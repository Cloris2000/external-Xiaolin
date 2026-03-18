# Cell-type GWAS pipeline (Nextflow)

Nextflow pipelines for genome-wide association studies (GWAS) of cell-type proportions across multiple brain cohorts. Phenotypes are estimated cell-type proportions from bulk RNA-seq (deconvolution); genotypes are from WGS or SNP array data. The pipeline runs per-cohort GWAS with REGENIE and meta-analysis with METAL.

## Cohorts

Supported cohorts: **ROSMAP**, **Mayo**, **MSBB**, **CMC_MSSM**, **CMC_PENN**, **CMC_PITT**, **GTEx**, **NABEC**.

- WGS cohorts: ROSMAP, Mayo, MSBB, GTEx, NABEC  
- SNP array (reimputed): CMC_MSSM, CMC_PENN, CMC_PITT  

## Pipeline overview

### Main pipeline (`nextflow/combined_pipeline_v2.nf`)

Single cohort, end-to-end run in this order:

1. **Stage 1 — RNA-seq:** PCA, technical covariate removal, cell-type deconvolution (MGP).
2. **Stage 2 — Phenotype prep:** RINT-transformed cell proportions, sample list, clinical covariates (sex, age). No genotype data used here (avoids circular dependency).
3. **Stage 3 — Genotype QC:** Filter to phenotyped samples, variant QC (MAF, HWE, missingness), LD pruning, heterozygosity filter, PCA on genotypes.
4. **Stage 4 — Covariate prep:** Merge genotype PCA (PC1–PC10) with clinical covariates.
5. **Stage 5 — GWAS:** REGENIE step1 (whole-genome regression), REGENIE step2 (association testing per cell type), LOG10P→P conversion, METAL meta-analysis (if multiple cohorts are configured).

Outputs: per-cohort REGENIE results under `results/<cohort>/regenie_step2/` (e.g. `*_SST_step2.regenie.raw_p`).

### Standalone meta-analysis (`nextflow/standalone_meta_analysis.nf`)

Takes existing `*.regenie.raw_p` files from `results/<cohort>/regenie_step2/`, groups by cell type, and runs METAL. Use this when all cohorts have already been run and you only want to update meta-analysis.

**Cell types (19):** Astrocyte, Endothelial, IT, L4.IT, L5.ET, L5.6.IT.Car3, L5.6.NP, L6.CT, L6b, LAMP5, Microglia, OPC, Oligodendrocyte, PAX6, PVALB, Pericyte, SST, VIP, VLMC.

## Repository layout

```
nextflow/
├── combined_pipeline_v2.nf      # Main pipeline (one cohort per run)
├── genotyping_qc_pipeline.nf    # Genotype QC subworkflow
├── gwas_pipeline.nf             # REGENIE + METAL subworkflow
├── standalone_meta_analysis.nf  # Meta-analysis only
├── modules/                     # Nextflow process definitions
│   ├── pheno_prep.nf
│   ├── cov_prep.nf
│   ├── cell_type_deconv.nf
│   ├── pca_tech_cov.nf
│   ├── remove_tech_covar.nf
│   ├── regenie_step1.nf
│   ├── regenie_step2.nf
│   ├── convert_logp_to_p.nf
│   ├── metal_meta_analysis.nf
│   └── ...
├── scripts/                     # R and Python scripts called by modules
│   ├── pheno_prep.R
│   ├── cov_prep.R
│   ├── cell_type_deconv.R
│   ├── genotyping_qc.py
│   ├── remove_tech_covar.R
│   └── ...
└── nextflow.config.combined.*   # Per-cohort configs (e.g. .nabec, .rosmap, .gtex)
```

## Requirements

- **Nextflow** (tested with 25.x)
- **R** (with packages used in `scripts/*.R`, e.g. for deconvolution, covariates, plotting)
- **Python 3** (for `scripts/genotyping_qc.py`)
- **PLINK2**, **REGENIE**, **METAL** (paths set in config)
- **SLURM** (or change `executor` in config for local/other schedulers)

Conda/env details and exact R/Python package lists are defined per environment; configs reference absolute paths to binaries and data.

## Running the pipeline

**One cohort (e.g. NABEC):**

```bash
cd nextflow
nextflow run combined_pipeline_v2.nf -c nextflow.config.combined.nabec
```

Resume after interruption:

```bash
nextflow run combined_pipeline_v2.nf -c nextflow.config.combined.nabec -resume
```

**Meta-analysis only** (after all cohort runs have produced `*.regenie.raw_p`):

```bash
nextflow run standalone_meta_analysis.nf -c nextflow.config.standalone_meta
```

Config files (`nextflow.config.combined.<cohort>`) define cohort-specific inputs: expression matrices, metadata, VCF paths (or pgen for CMC), reference panels, and output directories. Edit those configs to match your environment and paths.

### Running multi-cohort

The pipeline is one cohort per Nextflow run. To run all 8 cohorts and then meta-analysis, use the wrapper scripts (run from the `nextflow/` directory):

**All 8 cohorts in parallel** (fastest; ~4–6 hours wall time):

```bash
cd nextflow
bash run_8_cohorts_parallel.sh
```

This starts 8 separate Nextflow runs (one per cohort), each with its own config and work directory, then runs `standalone_meta_analysis.nf` after all complete.

**All 8 cohorts sequentially** (~12–16 hours):

```bash
cd nextflow
bash run_8_cohorts_sequential.sh
```

**Summary**

| Goal | Command |
|------|--------|
| One cohort | `nextflow run combined_pipeline_v2.nf -c nextflow.config.combined.<cohort>` |
| All 8 cohorts + meta-analysis | `bash run_8_cohorts_parallel.sh` or `run_8_cohorts_sequential.sh` |
| Meta-analysis only | `nextflow run standalone_meta_analysis.nf -c nextflow.config.standalone_meta` |

## Outputs

- **Per cohort:** `results/<cohort>/phenotypes_RINT.txt`, `clinical_covariates.txt`, `covariates.txt`, `pca.csv`, genotype QC outputs, `regenie_step2/*.regenie.raw_p`.
- **Meta-analysis:** `results/meta_analysis/<CellType>_meta_analysis_<cohorts>.tbl` (METAL summary statistics).

## Benchmarks

**Technical covariate vs. MGP accuracy:** Benchmark comparing MGP cell-type proportion accuracy when different technical covariate sets are regressed out before deconvolution (none vs. shared vs. full), using CMC_MSSM and ROSMAP with snRNA-based ground truth. Code, run instructions, and raw data paths live in the **benchmark_tech_cov** folder:

- [Tech covariate benchmark – README and run guide](nextflow/benchmark_tech_cov/README.md)
- [Tech covariate benchmark – raw data manifest](nextflow/benchmark_tech_cov/DATA.md)

## Citation and data

Cohort data (ROSMAP, Mayo, MSBB, CMC, GTEx, NABEC) are from their respective consortia and data use agreements; this repository contains only the pipeline code and configuration examples, not the data.

## License

See repository license file (if present). Pipeline code is for research use.
