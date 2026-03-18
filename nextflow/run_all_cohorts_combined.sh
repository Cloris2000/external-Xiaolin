#!/bin/bash

#
# Run Combined Pipeline for All 8 Cohorts
#
# This script runs the combined pipeline (phenotype + genotyping QC + GWAS)
# for all 8 cohorts. Supports both SEQUENTIAL and PARALLEL execution modes.
#
# Usage:
#   # Sequential (one cohort at a time - safer, less resource intensive)
#   bash run_all_cohorts_combined.sh
#
#   # Parallel (all cohorts simultaneously - faster, requires more resources)
#   bash run_all_cohorts_combined.sh --parallel
#
#   # Run specific cohorts only (sequential)
#   bash run_all_cohorts_combined.sh rosmap mayo msbb
#
#   # Run specific cohorts in parallel
#   bash run_all_cohorts_combined.sh --parallel rosmap mayo msbb
#

set -e  # Exit on error

# Set working directory to the nextflow directory
cd "$(dirname "$0")"
NEXTFLOW_DIR=$(pwd)

# Define all cohorts (names as they appear in config files)
ALL_COHORTS=(
    "rosmap"
    "mayo"
    "msbb"
    "gtex"
    "nabec"
    "cmc_mssm"
    "cmc_penn"
    "cmc_pitt"
)

# Parse arguments
PARALLEL_MODE=false
COHORTS_TO_RUN=()

for arg in "$@"; do
    if [[ "$arg" == "--parallel" ]] || [[ "$arg" == "-p" ]]; then
        PARALLEL_MODE=true
    else
        COHORTS_TO_RUN+=("$arg")
    fi
done

# Use provided cohorts or all cohorts
if [ ${#COHORTS_TO_RUN[@]} -eq 0 ]; then
    COHORTS_TO_RUN=("${ALL_COHORTS[@]}")
fi

# Log file
LOG_DIR="${NEXTFLOW_DIR}/logs"
mkdir -p "${LOG_DIR}"
MASTER_LOG="${LOG_DIR}/combined_pipeline_all_cohorts_$(date +%Y%m%d_%H%M%S).log"

echo "==================================================" | tee -a "$MASTER_LOG"
echo "Combined Pipeline - Multi-Cohort Execution" | tee -a "$MASTER_LOG"
echo "==================================================" | tee -a "$MASTER_LOG"
echo "Cohorts to process: ${COHORTS_TO_RUN[@]}" | tee -a "$MASTER_LOG"
echo "Execution mode: $([ "$PARALLEL_MODE" = true ] && echo "PARALLEL" || echo "SEQUENTIAL")" | tee -a "$MASTER_LOG"
echo "Start time: $(date)" | tee -a "$MASTER_LOG"
echo "==================================================" | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"

if [ "$PARALLEL_MODE" = true ]; then
    echo "⚠️  PARALLEL MODE WARNINGS:" | tee -a "$MASTER_LOG"
    echo "  - Requires significant compute resources" | tee -a "$MASTER_LOG"
    echo "  - Each cohort needs ~24-32 GB RAM" | tee -a "$MASTER_LOG"
    echo "  - Running 8 cohorts = ~192-256 GB RAM total" | tee -a "$MASTER_LOG"
    echo "  - Monitor system resources carefully!" | tee -a "$MASTER_LOG"
    echo "" | tee -a "$MASTER_LOG"
fi

# Counter for tracking progress
TOTAL_COHORTS=${#COHORTS_TO_RUN[@]}
FAILED_COHORTS=()
SUCCESS_COHORTS=()
declare -A COHORT_PIDS

# Function to run pipeline for a single cohort
run_cohort() {
    local cohort=$1
    local cohort_num=$2
    
    if [ "$PARALLEL_MODE" = false ]; then
        echo "" | tee -a "$MASTER_LOG"
        echo "==================================================" | tee -a "$MASTER_LOG"
        echo "[$cohort_num/$TOTAL_COHORTS] Processing: $cohort" | tee -a "$MASTER_LOG"
        echo "==================================================" | tee -a "$MASTER_LOG"
        echo "Start time: $(date)" | tee -a "$MASTER_LOG"
    fi
    
    # Config file name (cohort names are already lowercase)
    local config_file="nextflow.config.combined.${cohort}"
    
    # Check if config file exists
    if [ ! -f "$config_file" ]; then
        echo "❌ WARNING: Config file not found: $config_file" | tee -a "$MASTER_LOG"
        echo "           Please create this file before running the pipeline for $cohort" | tee -a "$MASTER_LOG"
        echo "           You can use nextflow.config.combined.cmc_mssm as a template" | tee -a "$MASTER_LOG"
        FAILED_COHORTS+=("$cohort (no config)")
        return 1
    fi
    
    # Create cohort-specific log file
    local cohort_log="${LOG_DIR}/combined_${cohort}_$(date +%Y%m%d_%H%M%S).log"
    
    # Run the pipeline
    if [ "$PARALLEL_MODE" = false ]; then
        echo "Running: nextflow run combined_pipeline.nf -c $config_file -resume" | tee -a "$MASTER_LOG"
    fi
    
    if nextflow run combined_pipeline.nf \
        -c "$config_file" \
        -resume \
        > "$cohort_log" 2>&1; then
        
        if [ "$PARALLEL_MODE" = false ]; then
            echo "" | tee -a "$MASTER_LOG"
            echo "✅ SUCCESS: $cohort pipeline completed" | tee -a "$MASTER_LOG"
            echo "End time: $(date)" | tee -a "$MASTER_LOG"
        fi
        SUCCESS_COHORTS+=("$cohort")
        return 0
    else
        if [ "$PARALLEL_MODE" = false ]; then
            echo "" | tee -a "$MASTER_LOG"
            echo "❌ ERROR: $cohort pipeline failed" | tee -a "$MASTER_LOG"
            echo "End time: $(date)" | tee -a "$MASTER_LOG"
            echo "Check log file: $cohort_log" | tee -a "$MASTER_LOG"
        fi
        FAILED_COHORTS+=("$cohort")
        return 1
    fi
}

# Run cohorts
if [ "$PARALLEL_MODE" = true ]; then
    # PARALLEL EXECUTION
    echo "Starting all cohorts in parallel..." | tee -a "$MASTER_LOG"
    echo "" | tee -a "$MASTER_LOG"
    
    # Launch all cohorts in background
    cohort_num=0
    for cohort in "${COHORTS_TO_RUN[@]}"; do
        cohort_num=$((cohort_num + 1))
        echo "[$cohort_num/$TOTAL_COHORTS] Launching $cohort in background..." | tee -a "$MASTER_LOG"
        
        # Run in background and capture PID
        run_cohort "$cohort" "$cohort_num" &
        COHORT_PIDS[$cohort]=$!
    done
    
    echo "" | tee -a "$MASTER_LOG"
    echo "All cohorts launched. Waiting for completion..." | tee -a "$MASTER_LOG"
    echo "Monitor logs in: $LOG_DIR/combined_*.log" | tee -a "$MASTER_LOG"
    echo "" | tee -a "$MASTER_LOG"
    
    # Wait for all cohorts to complete
    for cohort in "${COHORTS_TO_RUN[@]}"; do
        pid=${COHORT_PIDS[$cohort]}
        if wait $pid; then
            echo "✅ $cohort completed successfully" | tee -a "$MASTER_LOG"
        else
            echo "❌ $cohort failed" | tee -a "$MASTER_LOG"
        fi
    done
    
else
    # SEQUENTIAL EXECUTION
    cohort_num=0
    for cohort in "${COHORTS_TO_RUN[@]}"; do
        cohort_num=$((cohort_num + 1))
        run_cohort "$cohort" "$cohort_num" | tee -a "$MASTER_LOG"
    done
fi

# Print summary
echo "" | tee -a "$MASTER_LOG"
echo "==================================================" | tee -a "$MASTER_LOG"
echo "PIPELINE EXECUTION SUMMARY" | tee -a "$MASTER_LOG"
echo "==================================================" | tee -a "$MASTER_LOG"
echo "Execution mode: $([ "$PARALLEL_MODE" = true ] && echo "PARALLEL" || echo "SEQUENTIAL")" | tee -a "$MASTER_LOG"
echo "Total cohorts: $TOTAL_COHORTS" | tee -a "$MASTER_LOG"
echo "Successful: ${#SUCCESS_COHORTS[@]}" | tee -a "$MASTER_LOG"
echo "Failed: ${#FAILED_COHORTS[@]}" | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"

if [ ${#SUCCESS_COHORTS[@]} -gt 0 ]; then
    echo "Successful cohorts:" | tee -a "$MASTER_LOG"
    for cohort in "${SUCCESS_COHORTS[@]}"; do
        echo "  ✅ $cohort" | tee -a "$MASTER_LOG"
    done
    echo "" | tee -a "$MASTER_LOG"
fi

if [ ${#FAILED_COHORTS[@]} -gt 0 ]; then
    echo "Failed cohorts:" | tee -a "$MASTER_LOG"
    for cohort in "${FAILED_COHORTS[@]}"; do
        echo "  ❌ $cohort" | tee -a "$MASTER_LOG"
    done
    echo "" | tee -a "$MASTER_LOG"
fi

echo "End time: $(date)" | tee -a "$MASTER_LOG"
echo "Master log: $MASTER_LOG" | tee -a "$MASTER_LOG"
echo "Individual logs: $LOG_DIR/combined_*.log" | tee -a "$MASTER_LOG"
echo "==================================================" | tee -a "$MASTER_LOG"

# Exit with error if any cohort failed
if [ ${#FAILED_COHORTS[@]} -gt 0 ]; then
    exit 1
fi

exit 0
