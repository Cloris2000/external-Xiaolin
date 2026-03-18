#!/bin/bash
# Run the tech-covariate benchmark from start to end without Nextflow.
# Requires: R with packages (DESeq2, limma, PCAtools, markerGeneProfile, etc.), and config paths set.
# 1. Source config:  source benchmark_tech_cov/config_standalone.sh
# 2. Edit config_standalone.sh with your paths (count matrices, metadata, marker file, ground truth).
# 3. From repo root: bash benchmark_tech_cov/run_standalone_benchmark.sh [all|CMC_MSSM|ROSMAP] [all|none|shared|full]

set -euo pipefail

BENCHMARK_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${BENCHMARK_DIR}/config_standalone.sh"

TARGET_COHORT="${1:-all}"
TARGET_ARM="${2:-all}"

cohorts() {
    case "$TARGET_COHORT" in all) echo "CMC_MSSM ROSMAP" ;; CMC_MSSM|ROSMAP) echo "$TARGET_COHORT" ;; *) echo "Usage: $0 [all|CMC_MSSM|ROSMAP] [all|none|shared|full]" >&2; exit 1 ;; esac
}
arms() {
    case "$TARGET_ARM" in all) echo "none shared full" ;; none|shared|full) echo "$TARGET_ARM" ;; *) echo "Usage: $0 [all|CMC_MSSM|ROSMAP] [all|none|shared|full]" >&2; exit 1 ;; esac
}

run_pca() {
    local cohort="$1"
    local out_pca="$2"
    local count_file="$3"
    local meta_file="$4"
    local tissue="$5"
    local col_id="$6"
    mkdir -p "$out_pca"
    echo "=== PCA + z-score: $cohort ==="
    ( cd "$out_pca" && Rscript "${BENCHMARK_DIR}/scripts/pca_tech_cov.R" \
        --count_matrix_file "$count_file" \
        --metadata_file "$meta_file" \
        --tissue_filter "$tissue" \
        --output_dir . \
        --zscore_output "zscore_data.RData" \
        --metadata_output "metadata_DLPFC.RData" \
        --metrics_output "combined_metrics.csv" \
        $([ -n "$col_id" ] && echo "--col_sample_id_for_matching $col_id") )
}

run_remove_tech_covar() {
    local cohort="$1"
    local arm="$2"
    local zscore="$3"
    local meta_rdata="$4"
    local out_dir="$5"
    local batch_cov="$6"
    local top_n="${7:-10}"
    mkdir -p "$out_dir"
    local mode="auto_top_n"
    local fixed_file=""
    local disable_batch=""
    case "$arm" in
        none)   mode="none"; disable_batch="--disable_batch_correction" ;;
        shared) mode="fixed_list"
               if [ "$cohort" = "CMC_MSSM" ]; then fixed_file="${BENCHMARK_DIR}/benchmark_configs/tech_cov_sets/shared_cmc_mssm.txt"; else fixed_file="${BENCHMARK_DIR}/benchmark_configs/tech_cov_sets/shared_rosmap.txt"; fi ;;
        full)   mode="auto_top_n" ;;
        *) echo "Unknown arm: $arm" >&2; exit 1 ;;
    esac
    echo "=== Remove tech cov: $cohort / $arm ==="
    ( cd "$out_dir" && Rscript "${BENCHMARK_DIR}/scripts/remove_tech_covar.R" \
        --zscore_data "$zscore" \
        --metadata "$meta_rdata" \
        --top_n_tech_cov "$top_n" \
        --tech_cov_mode "$mode" \
        $([ -n "$fixed_file" ] && echo "--tech_covariates_file $fixed_file") \
        --batch_covariates "$batch_cov" \
        $disable_batch \
        --output_dir . \
        --corrected_output "corrected_data.RData" \
        --metadata_output "metadata_cleaned.csv" \
        --tech_cov_output "top_tech_covariates.txt" )
}

run_deconv() {
    local corrected="$1"
    local meta_cleaned="$2"
    local out_dir="$3"
    echo "=== Cell type deconvolution ==="
    ( cd "$out_dir" && Rscript "${BENCHMARK_DIR}/scripts/cell_type_deconv.R" \
        --corrected_data "$corrected" \
        --metadata "$meta_cleaned" \
        --deconv_tool "MGP" \
        --reference_taxonomy "sonny_markers" \
        --marker_file "$MARKER_FILE" \
        --hgnc_mapping_file "$HGNC_MAPPING_FILE" \
        --output_dir . \
        --proportions_output "cell_proportions.csv" \
        --proportions_scaled_output "cell_proportions_scaled.csv" )
}

for cohort in $(cohorts); do
    case "$cohort" in
        CMC_MSSM) count="$CMC_COUNT_MATRIX"; meta="$CMC_METADATA"; tissue="$CMC_TISSUE_FILTER"; col_id="$CMC_COL_SAMPLE_ID"; batch="$CMC_BATCH_COVARIATES" ;;
        ROSMAP)   count="$ROSMAP_COUNT_MATRIX"; meta="$ROSMAP_METADATA"; tissue="$ROSMAP_TISSUE_FILTER"; col_id="$ROSMAP_COL_SAMPLE_ID"; batch="$ROSMAP_BATCH_COVARIATES" ;;
        *) exit 1 ;;
    esac
    out_pca="${RESULTS_ROOT}/${cohort}/_pca"
    run_pca "$cohort" "$out_pca" "$count" "$meta" "$tissue" "$col_id"
    zscore="${out_pca}/zscore_data.RData"
    meta_rdata="${out_pca}/metadata_DLPFC.RData"
    for arm in $(arms); do
        out_arm="${RESULTS_ROOT}/${cohort}/${arm}"
        run_remove_tech_covar "$cohort" "$arm" "$zscore" "$meta_rdata" "$out_arm" "$batch"
        cp "${out_pca}/combined_metrics.csv" "$out_arm/"
        run_deconv "${out_arm}/corrected_data.RData" "${out_arm}/metadata_cleaned.csv" "$out_arm"
    done
done

echo "=== Running accuracy comparison (benchmark summary) ==="
for cohort in $(cohorts); do
    Rscript "${BENCHMARK_DIR}/scripts/benchmark_tech_cov_accuracy.R" \
        --cohort "$cohort" \
        --results_root "${RESULTS_ROOT}/${cohort}" \
        --ground_truth_cmc "$GROUND_TRUTH_CMC" \
        --ground_truth_rosmap "$GROUND_TRUTH_ROSMAP"
done
echo "Done. Results under ${RESULTS_ROOT} and per-cohort benchmark_summary/."
