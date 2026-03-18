#!/bin/bash
set -euo pipefail

# Run from nextflow/ (parent of benchmark_tech_cov). Override via NEXTFLOW_BENCHMARK_PROJECT_DIR or NEXTFLOW_BIN.
BENCHMARK_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="${NEXTFLOW_BENCHMARK_PROJECT_DIR:-$(cd "$BENCHMARK_DIR/.." && pwd)}"
NEXTFLOW_BIN="${NEXTFLOW_BIN:-nextflow}"
WORK_ROOT="${PROJECT_DIR}/work_benchmark_tech_cov"

export JAVA_HOME="${JAVA_HOME:-$HOME/.sdkman/candidates/java/current}"
export PATH="$HOME/.local/bin:${JAVA_HOME}/bin:$PATH"

TARGET_COHORT="${1:-all}"
TARGET_ARM="${2:-all}"

usage() {
    echo "Usage: $0 [all|CMC_MSSM|ROSMAP] [all|none|shared|full]"
    echo "Run from nextflow/ directory: bash benchmark_tech_cov/run_tech_cov_benchmark.sh all all"
}

normalize_cohorts() {
    case "$TARGET_COHORT" in
        all)
            echo "CMC_MSSM ROSMAP"
            ;;
        CMC_MSSM|ROSMAP)
            echo "$TARGET_COHORT"
            ;;
        *)
            usage
            exit 1
            ;;
    esac
}

normalize_arms() {
    case "$TARGET_ARM" in
        all)
            echo "none shared full"
            ;;
        none|shared|full)
            echo "$TARGET_ARM"
            ;;
        *)
            usage
            exit 1
            ;;
    esac
}

config_for() {
    local cohort="$1"
    local arm="$2"
    case "${cohort}:${arm}" in
        CMC_MSSM:none)   echo "${BENCHMARK_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.cmc_mssm.none" ;;
        CMC_MSSM:shared) echo "${BENCHMARK_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.cmc_mssm.shared" ;;
        CMC_MSSM:full)   echo "${BENCHMARK_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.cmc_mssm.full" ;;
        ROSMAP:none)     echo "${BENCHMARK_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.rosmap.none" ;;
        ROSMAP:shared)   echo "${BENCHMARK_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.rosmap.shared" ;;
        ROSMAP:full)     echo "${BENCHMARK_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.rosmap.full" ;;
        *)
            echo "ERROR: Unknown cohort/arm: ${cohort}/${arm}" >&2
            exit 1
            ;;
    esac
}

run_arm() {
    local cohort="$1"
    local arm="$2"
    local config_file
    local work_dir
    local nextflow_home

    config_file="$(config_for "$cohort" "$arm")"
    work_dir="${WORK_ROOT}/${cohort}/${arm}"
    nextflow_home="/tmp/.nextflow_${cohort}_${arm}_$$"

    mkdir -p "$(dirname "$work_dir")"

    echo "================================================"
    echo "Running tech-cov benchmark arm: ${cohort} / ${arm}"
    echo "  config: ${config_file}"
    echo "  work:   ${work_dir}"
    echo "================================================"

    cd "$PROJECT_DIR"
    export NXF_HOME="$nextflow_home"

    "$NEXTFLOW_BIN" run phenotype_pipeline.nf \
        -c "$config_file" \
        -w "$work_dir" \
        -name "techcov_${cohort}_${arm}_$$"

    rm -rf "$nextflow_home"
}

run_summary() {
    local cohort="$1"
    echo "================================================"
    echo "Running benchmark summary for ${cohort}"
    echo "================================================"
    cd "$PROJECT_DIR"
    Rscript "${BENCHMARK_DIR}/scripts/benchmark_tech_cov_accuracy.R" --cohort "$cohort" \
        --results_root "${PROJECT_DIR}/results/benchmark_tech_cov/${cohort}"
}

pids=()
for cohort in $(normalize_cohorts); do
    for arm in $(normalize_arms); do
        run_arm "$cohort" "$arm" &
        pids+=($!)
    done
done
echo "Waiting for all benchmark arms to complete..."
for pid in "${pids[@]}"; do
    wait "$pid" || true
done
echo "All benchmark arms completed. Running summaries..."
for cohort in $(normalize_cohorts); do
    run_summary "$cohort"
done
