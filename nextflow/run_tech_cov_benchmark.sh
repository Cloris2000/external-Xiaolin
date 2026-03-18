#!/bin/bash
set -euo pipefail

PROJECT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
NEXTFLOW_BIN="/nethome/kcni/xzhou/.local/bin/nextflow"
WORK_ROOT="${PROJECT_DIR}/work_benchmark_tech_cov"

export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$HOME/.local/bin:$JAVA_HOME/bin:$PATH"

TARGET_COHORT="${1:-all}"
TARGET_ARM="${2:-all}"

usage() {
    echo "Usage: $0 [all|CMC_MSSM|ROSMAP] [all|none|shared|full]"
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
        CMC_MSSM:none)
            echo "${PROJECT_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.cmc_mssm.none"
            ;;
        CMC_MSSM:shared)
            echo "${PROJECT_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.cmc_mssm.shared"
            ;;
        CMC_MSSM:full)
            echo "${PROJECT_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.cmc_mssm.full"
            ;;
        ROSMAP:none)
            echo "${PROJECT_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.rosmap.none"
            ;;
        ROSMAP:shared)
            echo "${PROJECT_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.rosmap.shared"
            ;;
        ROSMAP:full)
            echo "${PROJECT_DIR}/benchmark_configs/nextflow.config.tech_cov_benchmark.rosmap.full"
            ;;
        *)
            echo "ERROR: Unknown cohort/arm combination: ${cohort}/${arm}" >&2
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
    echo "Running tech-cov benchmark arm"
    echo "  cohort: ${cohort}"
    echo "  arm:    ${arm}"
    echo "  config: ${config_file}"
    echo "  work:   ${work_dir}"
    echo "================================================"

    cd "$PROJECT_DIR"

    # Use isolated Nextflow home to avoid cache conflicts
    export NXF_HOME="$nextflow_home"

    "$NEXTFLOW_BIN" run phenotype_pipeline.nf \
        -c "$config_file" \
        -w "$work_dir" \
        -name "techcov_${cohort}_${arm}_$$"
    
    # Clean up isolated Nextflow home
    rm -rf "$nextflow_home"
}

run_summary() {
    local cohort="$1"

    echo "================================================"
    echo "Running benchmark summary for ${cohort}"
    echo "================================================"

    cd "$PROJECT_DIR"

    Rscript "${PROJECT_DIR}/scripts/benchmark_tech_cov_accuracy.R" \
        --cohort "$cohort"
}

pids=()

# Launch all arms in parallel
for cohort in $(normalize_cohorts); do
    for arm in $(normalize_arms); do
        run_arm "$cohort" "$arm" &
        pids+=($!)
    done
done

# Wait for all arms to complete
echo "Waiting for all benchmark arms to complete..."
for pid in "${pids[@]}"; do
    wait "$pid"
    if [ $? -ne 0 ]; then
        echo "WARNING: Arm process $pid failed"
    fi
done

echo "All benchmark arms completed. Running summaries..."

# Run summaries after all arms are done
for cohort in $(normalize_cohorts); do
    run_summary "$cohort"
done
