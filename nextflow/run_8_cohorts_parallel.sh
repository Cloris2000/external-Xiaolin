#!/bin/bash
#
# Run All 8 Cohorts in PARALLEL + Meta-Analysis
#
# This script launches all cohorts simultaneously using background jobs.
# Much faster than sequential processing!
#
# Fixes applied (from past error logs):
#   - Each cohort uses its own NXF_HOME (work_parallel/<cohort>/.nextflow_home)
#     so multiple Nextflow runs do not contend on .nextflow/history lock.
#   - Counter updates use COMPLETED=$((COMPLETED+1)) instead of ((COMPLETED++))
#     so that set -e does not exit when the counter was 0.
#
# Estimated time:
#   - Sequential: 12-16 hours
#   - Parallel: 4-6 hours (depending on resources)
#

set -e
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Cohorts to process
COHORTS=("rosmap" "mayo" "msbb" "cmc_mssm" "cmc_penn" "cmc_pitt" "gtex" "nabec")

# Log directory
LOG_DIR="${SCRIPT_DIR}/logs/parallel_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║        Run All 8 Cohorts in PARALLEL + Meta-Analysis          ║"
echo "╠════════════════════════════════════════════════════════════════╣"
echo "║  Approach: Launch all cohorts simultaneously                   ║"
echo "║                                                                ║"
echo "║  Cohorts (8): ROSMAP, Mayo, MSBB, CMC_MSSM, CMC_PENN,          ║"
echo "║               CMC_PITT, GTEx, NABEC                            ║"
echo "║                                                                ║"
echo "║  Execution: All cohorts run in parallel as background jobs    ║"
echo "║  Meta-analysis: Runs after ALL cohorts complete                ║"
echo "║                                                                ║"
echo "║  Estimated time: 4-6 hours (vs 12-16h sequential)             ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Log directory: ${LOG_DIR}"
echo ""

# Array to track PIDs
declare -a PIDS
declare -A COHORT_PIDS

TOTAL_START=$(date +%s)

echo "=========================================="
echo "Launching All Cohorts in Parallel"
echo "=========================================="

# Launch all cohorts in background
for cohort in "${COHORTS[@]}"; do
    CONFIG="nextflow.config.combined.${cohort}"
    
    if [ ! -f "$CONFIG" ]; then
        echo "⚠ WARNING: Config not found: ${CONFIG}"
        echo "Skipping ${cohort}..."
        continue
    fi
    
    echo "🚀 Launching ${cohort^^}..."
    echo "   Config: ${CONFIG}"
    echo "   Log: ${LOG_DIR}/${cohort}.log"
    
    # Isolate each cohort's Nextflow runtime to avoid .nextflow/history lock conflicts.
    # NXF_WORK = work dir for task outputs; NXF_HOME = Nextflow runtime (history, cache).
    COHORT_WORK_DIR="${SCRIPT_DIR}/work_parallel/${cohort}"
    COHORT_NXF_HOME="${COHORT_WORK_DIR}/.nextflow_home"
    mkdir -p "$COHORT_WORK_DIR" "$COHORT_NXF_HOME"
    
    # Launch in background with cohort-specific NXF_HOME and work dir.
    # NXF_HOME isolates each cohort's .nextflow/history so there is no lock contention.
    (
        export NXF_HOME="$COHORT_NXF_HOME"
        nextflow run combined_pipeline_v2.nf \
            -c "$CONFIG" \
            -work-dir "$COHORT_WORK_DIR" \
            -resume \
            > "${LOG_DIR}/${cohort}.log" 2> "${LOG_DIR}/${cohort}.err"
    ) &
    
    # Store PID
    PID=$!
    PIDS+=($PID)
    COHORT_PIDS[$PID]=$cohort
    
    # Small delay to avoid simultaneous starts
    sleep 2
done

echo ""
echo "All cohorts launched!"
echo "Running: ${#PIDS[@]} cohorts"
echo ""
echo "=========================================="
echo "Monitoring Progress"
echo "=========================================="
echo "You can monitor individual cohorts:"
for cohort in "${COHORTS[@]}"; do
    echo "  tail -f ${LOG_DIR}/${cohort}.log"
done
echo ""
echo "Waiting for all cohorts to complete..."
echo "(This may take 4-6 hours depending on cohort sizes)"
echo ""

# Monitor completion
COMPLETED=0
FAILED=0
SUCCESS_COHORTS=()
FAILED_COHORTS=()

# Wait for each background job and track success/failure
for PID in "${PIDS[@]}"; do
    COHORT="${COHORT_PIDS[$PID]}"
    
    # Wait for this specific PID
    if wait $PID; then
        echo "✓ ${COHORT^^} completed successfully (PID: $PID)"
        COMPLETED=$((COMPLETED + 1))
        SUCCESS_COHORTS+=("$COHORT")
    else
        echo "✗ ${COHORT^^} FAILED (PID: $PID)"
        echo "  Check error log: ${LOG_DIR}/${COHORT}.err"
        FAILED=$((FAILED + 1))
        FAILED_COHORTS+=("$COHORT")
    fi
done

COHORT_END=$(date +%s)
COHORT_DURATION=$((COHORT_END - TOTAL_START))
COHORT_HOURS=$((COHORT_DURATION / 3600))
COHORT_MINUTES=$(((COHORT_DURATION % 3600) / 60))

echo ""
echo "=========================================="
echo "Cohort Processing Complete"
echo "=========================================="
echo "Successful: ${COMPLETED}/${#PIDS[@]}"
echo "Failed: ${FAILED}/${#PIDS[@]}"
echo "Duration: ${COHORT_HOURS}h ${COHORT_MINUTES}m"
echo ""

if [ ${FAILED} -gt 0 ]; then
    echo "Failed cohorts:"
    for cohort in "${FAILED_COHORTS[@]}"; do
        echo "  - ${cohort^^}"
    done
    echo ""
fi

# Only run meta-analysis if at least 2 cohorts succeeded
if [ $COMPLETED -ge 2 ]; then
    echo ""
    echo "=========================================="
    echo "Running Meta-Analysis"
    echo "=========================================="
    echo "Combining results from ${COMPLETED} cohorts..."
    echo ""
    
    META_START=$(date +%s)
    META_MINUTES=0
    
    if nextflow run standalone_meta_analysis.nf \
        -c nextflow.config.standalone_meta \
        -resume \
        > "${LOG_DIR}/meta_analysis.log" 2> "${LOG_DIR}/meta_analysis.err"; then
        
        META_END=$(date +%s)
        META_DURATION=$((META_END - META_START))
        META_MINUTES=$((META_DURATION / 60))
        
        TOTAL_END=$(date +%s)
        TOTAL_DURATION=$((TOTAL_END - TOTAL_START))
        TOTAL_HOURS=$((TOTAL_DURATION / 3600))
        TOTAL_MINUTES=$(((TOTAL_DURATION % 3600) / 60))
        
        echo ""
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║                      ✓ SUCCESS ✓                              ║"
        echo "╠════════════════════════════════════════════════════════════════╣"
        echo "║  Cohorts processed: ${COMPLETED}/${#PIDS[@]} (parallel)                             ║"
        echo "║  Meta-analysis: COMPLETE                                       ║"
        echo "║                                                                ║"
        echo "║  Total duration: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m                                       ║"
        echo "║  Time saved vs sequential: ~$(( (12*60+COHORT_HOURS*60+COHORT_MINUTES) / 60 - TOTAL_HOURS ))h                               ║"
        echo "║                                                                ║"
        echo "║  Results:                                                      ║"
        echo "║    - Individual cohorts: results/{COHORT}/regenie_step2/     ║"
        echo "║    - Meta-analysis: results/meta_analysis/                    ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        echo ""
        
        # List meta-analysis results
        echo "Meta-analysis results (19 cell types):"
        ls -lh results/meta_analysis/*.tbl 2>/dev/null | wc -l | xargs echo "  Files:"
        echo ""
        
    else
        echo ""
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║                ⚠ META-ANALYSIS FAILED ⚠                       ║"
        echo "╠════════════════════════════════════════════════════════════════╣"
        echo "║  Cohorts processed: ${COMPLETED}/${#PIDS[@]}                                       ║"
        echo "║  Meta-analysis: FAILED                                         ║"
        echo "║                                                                ║"
        echo "║  Check error log: ${LOG_DIR}/meta_analysis.err        ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        echo ""
        exit 1
    fi
else
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║              ⚠ INSUFFICIENT COHORTS ⚠                         ║"
    echo "╠════════════════════════════════════════════════════════════════╣"
    echo "║  Only ${COMPLETED} cohort(s) succeeded.                                ║"
    echo "║  Need at least 2 cohorts for meta-analysis.                   ║"
    echo "║  Meta-analysis SKIPPED.                                        ║"
    echo "║                                                                ║"
    echo "║  Fix failed cohorts and re-run with:                          ║"
    echo "║    bash run_8_cohorts_parallel.sh                             ║"
    echo "║                                                                ║"
    echo "║  (Use -resume flag will skip successful cohorts)              ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    exit 1
fi

# Save summary
cat > "${LOG_DIR}/SUMMARY.txt" << EOF
PARALLEL MULTI-COHORT EXECUTION SUMMARY
========================================

Date: $(date)
Execution Mode: PARALLEL (all cohorts simultaneously)

Total Duration: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m
  - Cohort processing: ${COHORT_HOURS}h ${COHORT_MINUTES}m
  - Meta-analysis: ${META_MINUTES}m

Cohorts Processed: ${COMPLETED}/${#PIDS[@]}
Successful: $(IFS=, ; echo "${SUCCESS_COHORTS[*]}")
Failed: ${FAILED}

Pipeline: combined_pipeline_v2.nf (per cohort, in parallel)
Meta-analysis: standalone_meta_analysis.nf

Results:
  - Individual cohorts: results/{COHORT}/regenie_step2/
  - Meta-analysis: results/meta_analysis/

Logs Directory: ${LOG_DIR}
  - Per-cohort logs: ${LOG_DIR}/{cohort}.log
  - Meta-analysis log: ${LOG_DIR}/meta_analysis.log

Performance:
  Estimated sequential time: ~12-16 hours
  Actual parallel time: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m
  Time saved: ~$(( 14 - TOTAL_HOURS ))h

To check results:
  ls -lh results/*/regenie_step2/*.regenie
  ls -lh results/meta_analysis/*.tbl
EOF

echo "Summary saved: ${LOG_DIR}/SUMMARY.txt"
echo ""

exit 0

