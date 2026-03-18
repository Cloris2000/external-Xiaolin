# Standalone benchmark config (no Nextflow). Edit paths to match your environment.
# Source this file before running run_standalone_benchmark.sh, or export vars and run.

# === Shared reference data (used by both cohorts) ===
export MARKER_FILE="${MARKER_FILE:-/path/to/new_MTGnCgG_lfct2.5_Publication.csv}"
export HGNC_MAPPING_FILE="${HGNC_MAPPING_FILE:-/path/to/hgnc_complete_set.txt}"

# === CMC_MSSM ===
export CMC_COUNT_MATRIX="${CMC_COUNT_MATRIX:-/path/to/CMC_MSSM_count_matrix.csv}"
export CMC_METADATA="${CMC_METADATA:-/path/to/CMC_MSSM_metadata.csv}"
export CMC_TISSUE_FILTER="${CMC_TISSUE_FILTER:-dorsolateral prefrontal cortex}"
export CMC_COL_SAMPLE_ID="${CMC_COL_SAMPLE_ID:-individualID}"
export CMC_BATCH_COVARIATES="${CMC_BATCH_COVARIATES:-sequencingBatch,libraryPrep}"

# === ROSMAP ===
export ROSMAP_COUNT_MATRIX="${ROSMAP_COUNT_MATRIX:-/path/to/ROSMAP_DLPFC_batch_all.csv}"
export ROSMAP_METADATA="${ROSMAP_METADATA:-/path/to/ROSMAP_combined_metrics.csv}"
export ROSMAP_TISSUE_FILTER="${ROSMAP_TISSUE_FILTER:-dorsolateral prefrontal cortex}"
export ROSMAP_COL_SAMPLE_ID="${ROSMAP_COL_SAMPLE_ID:-}"
export ROSMAP_BATCH_COVARIATES="${ROSMAP_BATCH_COVARIATES:-sequencingBatch,libraryPrep}"

# === Ground truth (for final R comparison step) ===
export GROUND_TRUTH_CMC="${GROUND_TRUTH_CMC:-/path/to/PsychEncode_label_transferred_snCTP.csv}"
export GROUND_TRUTH_ROSMAP="${GROUND_TRUTH_ROSMAP:-/path/to/rosmap_single_nuc_proportions.csv}"

# === Output root (default: results_standalone next to this config) ===
BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-.}")" && pwd)"
export RESULTS_ROOT="${RESULTS_ROOT:-${BENCHMARK_DIR}/results_standalone}"
