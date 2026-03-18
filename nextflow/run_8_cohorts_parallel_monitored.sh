#!/bin/bash
#
# Run All 8 Cohorts in PARALLEL with Real-Time Monitoring
#
# Features:
#   - Live status updates every 30 seconds
#   - Clear indication of running/completed/failed cohorts
#   - Detailed progress tracking
#   - Logs errors immediately when detected
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

# Status file for real-time updates
STATUS_FILE="${LOG_DIR}/status.txt"
touch "$STATUS_FILE"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║     Run All 8 Cohorts in PARALLEL with Live Monitoring        ║"
echo "╠════════════════════════════════════════════════════════════════╣"
echo "║  Cohorts (8): ROSMAP, Mayo, MSBB, CMC_MSSM, CMC_PENN,          ║"
echo "║               CMC_PITT, GTEx, NABEC                            ║"
echo "║                                                                ║"
echo "║  Features:                                                     ║"
echo "║    ✓ Real-time status updates every 30 seconds                ║"
echo "║    ✓ Live progress monitoring                                 ║"
echo "║    ✓ Immediate error detection                                ║"
echo "║    ✓ Detailed logs per cohort                                 ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Log directory: ${LOG_DIR}"
echo "Status file: ${STATUS_FILE}"
echo ""

# Arrays to track everything
declare -a PIDS
declare -A COHORT_PIDS
declare -A COHORT_STATUS  # running, completed, failed
declare -A COHORT_START_TIME

TOTAL_START=$(date +%s)

# Function to update status display
update_status() {
    clear
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║              Multi-Cohort Pipeline Status                      ║"
    echo "╠════════════════════════════════════════════════════════════════╣"
    echo "  Started: $(date -d @$TOTAL_START '+%Y-%m-%d %H:%M:%S')"
    echo "  Elapsed: $(( ($(date +%s) - TOTAL_START) / 60 )) minutes"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    
    local running=0
    local completed=0
    local failed=0
    
    for cohort in "${COHORTS[@]}"; do
        local status="${COHORT_STATUS[$cohort]:-unknown}"
        local symbol="⏳"
        local color="$YELLOW"
        
        case "$status" in
            running)
                symbol="🔄"
                color="$BLUE"
                ((running++))
                ;;
            completed)
                symbol="✓"
                color="$GREEN"
                ((completed++))
                ;;
            failed)
                symbol="✗"
                color="$RED"
                ((failed++))
                ;;
        esac
        
        # Get elapsed time for this cohort
        local elapsed=""
        if [ -n "${COHORT_START_TIME[$cohort]:-}" ]; then
            local cohort_elapsed=$(( ($(date +%s) - ${COHORT_START_TIME[$cohort]}) / 60 ))
            elapsed=" (${cohort_elapsed}m)"
        fi
        
        # Check log for latest activity
        local latest=""
        if [ -f "${LOG_DIR}/${cohort}.log" ]; then
            latest=$(tail -1 "${LOG_DIR}/${cohort}.log" 2>/dev/null | cut -c1-50)
        fi
        
        printf "${color}%-8s${NC} ${symbol} %-15s %s\n" "${cohort^^}" "${status}${elapsed}" "${latest}"
    done
    
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    printf "Running: ${BLUE}%d${NC}  |  Completed: ${GREEN}%d${NC}  |  Failed: ${RED}%d${NC}  |  Total: %d\n" \
        $running $completed $failed ${#COHORTS[@]}
    echo "════════════════════════════════════════════════════════════════"
    echo ""
    echo "Press Ctrl+C to stop monitoring (cohorts will continue running)"
    echo "Logs: tail -f ${LOG_DIR}/{cohort}.log"
    echo ""
}

# Function to check cohort status
check_cohort_status() {
    local cohort=$1
    local pid=$2
    
    # Check if process is still running
    if kill -0 $pid 2>/dev/null; then
        # Still running - check for errors in log
        if grep -q "ERROR\|FAILED\|Exception" "${LOG_DIR}/${cohort}.err" 2>/dev/null; then
            COHORT_STATUS[$cohort]="failed"
            return 1
        else
            COHORT_STATUS[$cohort]="running"
            return 0
        fi
    else
        # Process finished - check exit code
        if wait $pid 2>/dev/null; then
            COHORT_STATUS[$cohort]="completed"
            return 0
        else
            COHORT_STATUS[$cohort]="failed"
            return 1
        fi
    fi
}

# Launch all cohorts
echo "Launching cohorts in parallel..."
echo ""

for cohort in "${COHORTS[@]}"; do
    CONFIG="nextflow.config.combined.${cohort}"
    
    if [ ! -f "$CONFIG" ]; then
        echo "⚠ WARNING: Config not found: ${CONFIG}"
        COHORT_STATUS[$cohort]="failed"
        continue
    fi
    
    echo "🚀 Launching ${cohort^^}..."
    
    # Create separate work directory for this cohort to avoid .nextflow/history lock conflicts
    COHORT_WORK_DIR="${SCRIPT_DIR}/work_parallel/${cohort}"
    mkdir -p "$COHORT_WORK_DIR"
    
    # Launch in background with separate work directory
    NXF_WORK="$COHORT_WORK_DIR" nextflow run combined_pipeline_v2.nf \
        -c "$CONFIG" \
        -work-dir "$COHORT_WORK_DIR" \
        -resume \
        > "${LOG_DIR}/${cohort}.log" 2> "${LOG_DIR}/${cohort}.err" &
    
    # Store info
    PID=$!
    PIDS+=($PID)
    COHORT_PIDS[$PID]=$cohort
    COHORT_STATUS[$cohort]="running"
    COHORT_START_TIME[$cohort]=$(date +%s)
    
    sleep 1
done

echo ""
echo "All cohorts launched! Starting monitoring..."
sleep 3

# Monitor loop
MONITORING=true
trap 'MONITORING=false; echo ""; echo "Stopping monitoring... (cohorts will continue)"; echo ""' INT

while $MONITORING; do
    # Update status for all cohorts
    ALL_DONE=true
    for PID in "${PIDS[@]}"; do
        cohort="${COHORT_PIDS[$PID]}"
        if check_cohort_status "$cohort" "$PID"; then
            if [ "${COHORT_STATUS[$cohort]}" == "running" ]; then
                ALL_DONE=false
            fi
        else
            # Failed - show error snippet
            if [ "${COHORT_STATUS[$cohort]}" == "failed" ]; then
                echo ""
                echo -e "${RED}✗ ${cohort^^} FAILED${NC}"
                if [ -f "${LOG_DIR}/${cohort}.err" ]; then
                    echo "Last error:"
                    tail -5 "${LOG_DIR}/${cohort}.err" | sed 's/^/  /'
                fi
                echo ""
            fi
        fi
    done
    
    # Display current status
    update_status
    
    # Break if all done
    if $ALL_DONE; then
        break
    fi
    
    # Wait before next update
    sleep 30
done

# Final status check
echo ""
echo "All cohorts finished processing!"
echo ""

# Count final results
COMPLETED=0
FAILED=0
SUCCESS_COHORTS=()
FAILED_COHORTS=()

for cohort in "${COHORTS[@]}"; do
    if [ "${COHORT_STATUS[$cohort]}" == "completed" ]; then
        ((COMPLETED++))
        SUCCESS_COHORTS+=("$cohort")
    else
        ((FAILED++))
        FAILED_COHORTS+=("$cohort")
    fi
done

COHORT_END=$(date +%s)
COHORT_DURATION=$((COHORT_END - TOTAL_START))
COHORT_HOURS=$((COHORT_DURATION / 3600))
COHORT_MINUTES=$(((COHORT_DURATION % 3600) / 60))

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║              Cohort Processing Complete                        ║"
echo "╠════════════════════════════════════════════════════════════════╣"
printf "║  Successful: ${GREEN}%-2d${NC}/8                                                ║\n" $COMPLETED
printf "║  Failed:     ${RED}%-2d${NC}/8                                                ║\n" $FAILED
echo "║  Duration:   ${COHORT_HOURS}h ${COHORT_MINUTES}m                                                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

if [ ${FAILED} -gt 0 ]; then
    echo -e "${RED}Failed cohorts:${NC}"
    for cohort in "${FAILED_COHORTS[@]}"; do
        echo "  ✗ ${cohort^^}"
        echo "    Log: ${LOG_DIR}/${cohort}.err"
    done
    echo ""
fi

# Only run meta-analysis if at least 2 cohorts succeeded
if [ $COMPLETED -ge 2 ]; then
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                   Running Meta-Analysis                        ║"
    echo "╠════════════════════════════════════════════════════════════════╣"
    echo "║  Combining results from ${COMPLETED} successful cohorts                   ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    
    META_START=$(date +%s)
    
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
        echo "║  Cohorts processed: ${COMPLETED}/8 (parallel)                             ║"
        echo "║  Meta-analysis: COMPLETE                                       ║"
        echo "║                                                                ║"
        echo "║  Total time: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m                                           ║"
        echo "║  Meta-analysis: ${META_MINUTES}m                                          ║"
        echo "║                                                                ║"
        echo "║  Results:                                                      ║"
        echo "║    Individual: results/{COHORT}/regenie_step2/                ║"
        echo "║    Meta-analysis: results/meta_analysis/                      ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        echo ""
        
        # Show meta-analysis results
        echo "Meta-analysis results (19 cell types):"
        ls -lh results/meta_analysis/*.tbl 2>/dev/null | head -5
        echo "  ... and more"
        echo ""
        
    else
        echo ""
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║                ⚠ META-ANALYSIS FAILED ⚠                       ║"
        echo "╠════════════════════════════════════════════════════════════════╣"
        echo "║  Check: ${LOG_DIR}/meta_analysis.err                  ║"
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
    echo "║  Need at least 2 for meta-analysis.                           ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    exit 1
fi

# Save final summary
cat > "${LOG_DIR}/SUMMARY.txt" << EOF
PARALLEL EXECUTION SUMMARY
==========================

Date: $(date)
Mode: Parallel with real-time monitoring

Duration: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m
  - Cohort processing: ${COHORT_HOURS}h ${COHORT_MINUTES}m
  - Meta-analysis: ${META_MINUTES}m

Results:
  Successful: ${COMPLETED}/8
  Failed: ${FAILED}/8

Successful cohorts: ${SUCCESS_COHORTS[@]:-none}
Failed cohorts: ${FAILED_COHORTS[@]:-none}

Logs: ${LOG_DIR}
  - Per-cohort: {cohort}.log / {cohort}.err
  - Meta-analysis: meta_analysis.log / meta_analysis.err

Results:
  - Individual: results/{COHORT}/regenie_step2/
  - Meta-analysis: results/meta_analysis/
EOF

echo "Summary: ${LOG_DIR}/SUMMARY.txt"
echo ""

exit 0

