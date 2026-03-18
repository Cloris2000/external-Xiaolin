# Technical covariate benchmark

Compare MGP cell-type proportion accuracy when different technical covariate sets are regressed out before deconvolution (none vs. shared vs. full), using CMC_MSSM and ROSMAP with snRNA-based ground truth.

## Contents of this folder

| Item | Description |
|------|-------------|
| `run_tech_cov_benchmark.sh` | Run all benchmark arms (Nextflow) then per-cohort R summaries. Run from **nextflow/** directory. |
| `scripts/benchmark_tech_cov_accuracy.R` | R script: compare arm outputs to ground truth, write CSVs and plots. |
| `benchmark_configs/` | Nextflow configs for each cohort × arm (none/shared/full) and `tech_cov_sets/` (shared covariate lists). |
| `README.md` | This file. |
| `DATA.md` | Raw data manifest (paths to ground truth, count matrices, metadata). |

The **Nextflow phenotype pipeline** (`phenotype_pipeline.nf`), **cohort configs** (`nextflow.config.combined.cmc_mssm`, `nextflow.config.combined.rosmap`), and **remove_tech_covar** module/script live in the parent **nextflow/** directory. This folder only contains benchmark-specific configs and the comparison script.

## How to run

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

## Outputs

- **Pipeline:** `nextflow/results/benchmark_tech_cov/<COHORT>/<arm>/` (e.g. `cell_proportions_scaled.csv`).
- **Summaries:** `nextflow/results/benchmark_tech_cov/<COHORT>/benchmark_summary/` (CSVs, bar plot, delta plot, heatmap).

## Raw data

See **DATA.md** in this folder for paths to ground truth CSVs, count matrices, and metadata. Edit cohort configs in **nextflow/** and/or pass `--ground_truth_*` and `--results_root` to the R script to match your environment.
