# Technical covariate benchmark

Compare MGP cell-type proportion accuracy when different technical covariate sets are regressed out before deconvolution (none vs. shared vs. full), using CMC_MSSM and ROSMAP with snRNA-based ground truth.

## What you need to run this (e.g. on SCC)

**This folder alone is not enough.** The benchmark script calls the **Nextflow phenotype pipeline**, which lives in the parent **nextflow/** directory. Your PI should:

1. **Clone or download the full repository** (or at least the entire **nextflow/** directory and everything inside it).
2. **Transfer that full nextflow/ tree to SCC** (or clone the repo on SCC).
3. On SCC, install/load Nextflow, R, and any required R packages; then run from **nextflow/** as in “How to run” below.

Required pieces outside this folder (all under **nextflow/**):

- `phenotype_pipeline.nf` — main pipeline
- `nextflow.config.combined`, `nextflow.config.combined.cmc_mssm`, `nextflow.config.combined.rosmap` — base and cohort configs (input data paths, etc.)
- `modules/` — e.g. `pca_tech_cov.nf`, `remove_tech_covar.nf`, `cell_type_deconv.nf`
- `scripts/` — e.g. `pca_tech_cov.R`, `remove_tech_covar.R`, `cell_type_deconv.R` (and any other R/Python scripts used by the pipeline)

So: **download all scripts and the Nextflow pipeline in the repo (the full nextflow/ directory), transfer to SCC, then run the benchmark from nextflow/.**

## Contents of this folder

| Item | Description |
|------|-------------|
| **Option A – Nextflow** | |
| `run_tech_cov_benchmark.sh` | Run all benchmark arms via Nextflow, then per-cohort R summaries. Run from **nextflow/** directory. |
| `benchmark_configs/` | Nextflow configs for each cohort × arm (none/shared/full) and `tech_cov_sets/` (shared covariate lists). |
| **Option B – Standalone (no Nextflow)** | |
| `config_standalone.sh` | Paths for count matrices, metadata, marker file, ground truth. Edit and source before running standalone. |
| `run_standalone_benchmark.sh` | Run full benchmark with R only: PCA → remove tech cov → deconvolution → accuracy comparison. No Nextflow. |
| `scripts/pca_tech_cov.R` | PCA and z-score normalization (same logic as pipeline). |
| `scripts/remove_tech_covar.R` | Remove batch effects and technical covariates (none / shared list / auto). |
| `scripts/cell_type_deconv.R` | MGP cell-type proportion estimation. |
| `scripts/benchmark_tech_cov_accuracy.R` | Compare arm outputs to ground truth, write CSVs and plots. |
| **Docs** | |
| `README.md` | This file. |
| `DATA.md` | Raw data manifest (paths to ground truth, count matrices, metadata). |

**Option A** still needs the parent **nextflow/** directory (pipeline, cohort configs). **Option B** is self-contained: everything needed is in this folder plus your data paths in `config_standalone.sh`.

## How to run

### Option A: With Nextflow (on SCC or wherever Nextflow runs)

1. **From the nextflow directory** (parent of `benchmark_tech_cov`), run:
   ```bash
   cd nextflow
   bash benchmark_tech_cov/run_tech_cov_benchmark.sh all all
   ```
   Or one cohort / one arm:
   ```bash
   bash benchmark_tech_cov/run_tech_cov_benchmark.sh ROSMAP all
   bash benchmark_tech_cov/run_tech_cov_benchmark.sh CMC_MSSM shared
   ```

2. **Summary only** (if pipeline outputs already exist under `nextflow/results/benchmark_tech_cov/`):
   ```bash
   cd nextflow
   Rscript benchmark_tech_cov/scripts/benchmark_tech_cov_accuracy.R --cohort ROSMAP \
     --results_root results/benchmark_tech_cov/ROSMAP
   ```

3. **Override paths** (e.g. for a different clone or machine):
   ```bash
   export NEXTFLOW_BENCHMARK_PROJECT_DIR="/path/to/nextflow"
   bash benchmark_tech_cov/run_tech_cov_benchmark.sh all all
   ```
   For ground truth and results paths, use `--results_root`, `--ground_truth_cmc`, `--ground_truth_rosmap` in the R script.

### Option B: Standalone (no Nextflow – all scripts in this folder)

Runs from start to end using only R and this folder. No Nextflow or parent repo required.

1. **Edit paths** in `config_standalone.sh`: count matrices, metadata, marker file, HGNC mapping, ground truth CSVs (see `DATA.md` for what you need).
2. **Source config and run** (from anywhere; can run on SCC):
   ```bash
   cd /path/to/benchmark_tech_cov   # or cd nextflow/benchmark_tech_cov if inside repo
   source config_standalone.sh     # optional if you export vars yourself
   bash run_standalone_benchmark.sh all all
   ```
   Or one cohort / one arm: `bash run_standalone_benchmark.sh ROSMAP all` or `bash run_standalone_benchmark.sh CMC_MSSM shared`.
3. **Outputs** go under `results_standalone/` (or `RESULTS_ROOT` if you set it in config):  
   `results_standalone/<COHORT>/_pca/`, `results_standalone/<COHORT>/none|shared|full/`, and `results_standalone/<COHORT>/benchmark_summary/`.

R must have the same packages as the Nextflow pipeline (e.g. DESeq2, limma, PCAtools, markerGeneProfile, dplyr, readr, ggplot2). No Java or Nextflow needed.

## Outputs

- **Option A (Nextflow):** `nextflow/results/benchmark_tech_cov/<COHORT>/<arm>/` and `.../benchmark_summary/`.
- **Option B (Standalone):** `results_standalone/<COHORT>/<arm>/` and `results_standalone/<COHORT>/benchmark_summary/` (or `$RESULTS_ROOT/...` if set).

## Raw data

See **DATA.md** in this folder for paths to ground truth CSVs, count matrices, and metadata. Edit cohort configs in **nextflow/** and/or pass `--ground_truth_*` and `--results_root` to the R script to match your environment.
