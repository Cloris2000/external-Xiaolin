#!/bin/bash
#
# Run All 8 Cohorts in TRULY PARALLEL with Complete Isolation
#
# Each cohort gets its own:
#   - Nextflow cache directory (NXF_HOME)
#   - Work directory
#   - Temp directory
# This prevents ALL conflicts!
#

set -e
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Cohorts to process
COHORTS=("rosmap" "mayo" "msbb" "cmc_mssm" "cmc_penn" "cmc_pitt" "gtex" "nabec")
#COHORTS=("gtex")
# Log directory
LOG_DIR="${SCRIPT_DIR}/logs/truly_parallel_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"

# Isolated environments directory
ISOLATED_DIR="${SCRIPT_DIR}/isolated_runs"
mkdir -p "$ISOLATED_DIR"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║     Run All 8 Cohorts in TRULY PARALLEL (Isolated Envs)       ║"
echo "╠════════════════════════════════════════════════════════════════╣"
echo "║  Each cohort runs in completely isolated Nextflow environment ║"
echo "║    - Separate cache (NXF_HOME)                                 ║"
echo "║    - Separate work directory                                   ║"
echo "║    - Separate temp directory                                   ║"
echo "║    - NO conflicts possible!                                    ║"
echo "║                                                                ║"
echo "║  Cohorts (8): ROSMAP, Mayo, MSBB, CMC_MSSM, CMC_PENN,          ║"
echo "║               CMC_PITT, GTEx, NABEC                            ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Log directory: ${LOG_DIR}"
echo "Isolated environments: ${ISOLATED_DIR}"
echo ""

# Arrays to track PIDs
declare -a PIDS
declare -A COHORT_PIDS
declare -A COHORT_STATUS
declare -A COHORT_START_TIME

TOTAL_START=$(date +%s)

# Function to update status
update_status() {
    local running=0
    local completed=0
    local failed=0
    
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  Elapsed: $(( ($(date +%s) - TOTAL_START) / 60 )) minutes"
    echo "════════════════════════════════════════════════════════════════"
    
    for cohort in "${COHORTS[@]}"; do
        local status="${COHORT_STATUS[$cohort]:-unknown}"
        local symbol="⏳"
        local color="$YELLOW"
        
        case "$status" in
            running)
                symbol="🔄"
                color="$BLUE"
                ((running++))
                # Show last activity
                if [ -f "${LOG_DIR}/${cohort}.log" ]; then
                    local last_line=$(tail -1 "${LOG_DIR}/${cohort}.log" 2>/dev/null | cut -c1-60)
                    printf "${color}%-8s${NC} ${symbol} running    %s\n" "${cohort^^}" "$last_line"
                else
                    printf "${color}%-8s${NC} ${symbol} running\n" "${cohort^^}"
                fi
                ;;
            completed)
                symbol="✓"
                color="$GREEN"
                ((completed++))
                printf "${color}%-8s${NC} ${symbol} completed\n" "${cohort^^}"
                ;;
            failed)
                symbol="✗"
                color="$RED"
                ((failed++))
                printf "${color}%-8s${NC} ${symbol} FAILED\n" "${cohort^^}"
                ;;
            *)
                printf "%-8s %s starting...\n" "${cohort^^}" "$symbol"
                ;;
        esac
    done
    
    echo "════════════════════════════════════════════════════════════════"
    printf "Running: ${BLUE}%d${NC}  |  Completed: ${GREEN}%d${NC}  |  Failed: ${RED}%d${NC}\n" \
        $running $completed $failed
    echo "════════════════════════════════════════════════════════════════"
    echo ""
}

# Launch all cohorts with complete isolation
echo "Launching all cohorts in parallel with isolated environments..."
echo ""

for cohort in "${COHORTS[@]}"; do
    CONFIG="nextflow.config.combined.${cohort}"
    
    if [ ! -f "$CONFIG" ]; then
        echo "⚠ WARNING: Config not found: ${CONFIG}"
        COHORT_STATUS[$cohort]="failed"
        continue
    fi
    
    echo "🚀 Launching ${cohort^^}..."
    
    # Create completely isolated environment for this cohort
    COHORT_ENV="${ISOLATED_DIR}/${cohort}"
    mkdir -p "${COHORT_ENV}/cache"
    mkdir -p "${COHORT_ENV}/work"
    mkdir -p "${COHORT_ENV}/temp"
    mkdir -p "${COHORT_ENV}/launch"  # Separate launch directory
    
    # Link necessary files/directories to launch directory (symlinks for efficiency)
    ln -sf "${SCRIPT_DIR}/combined_pipeline_v2.nf" "${COHORT_ENV}/launch/"
    ln -sf "${SCRIPT_DIR}/genotyping_qc_pipeline.nf" "${COHORT_ENV}/launch/"  # Included by combined_pipeline_v2
    ln -sf "${SCRIPT_DIR}/gwas_pipeline.nf" "${COHORT_ENV}/launch/"
    ln -sf "${SCRIPT_DIR}/nextflow.config.combined" "${COHORT_ENV}/launch/"  # Base config
    ln -sf "${SCRIPT_DIR}/$CONFIG" "${COHORT_ENV}/launch/$(basename $CONFIG)"
    ln -sf "${SCRIPT_DIR}/modules" "${COHORT_ENV}/launch/"
    ln -sf "${SCRIPT_DIR}/scripts" "${COHORT_ENV}/launch/"
    ln -sf "${SCRIPT_DIR}/bin" "${COHORT_ENV}/launch/"
    
    # Launch with isolated environment variables from isolated directory
    (
        export NXF_HOME="${COHORT_ENV}/cache"
        export NXF_WORK="${COHORT_ENV}/work"
        export NXF_TEMP="${COHORT_ENV}/temp"
        export TMPDIR="${COHORT_ENV}/temp"
        
        cd "${COHORT_ENV}/launch"  # Launch from isolated directory!
        
        CONFIG_BASENAME="$(basename $CONFIG)"
        
        nextflow run combined_pipeline_v2.nf \
            -c "$CONFIG_BASENAME" \
            -work-dir "${COHORT_ENV}/work" \
            -resume \
            > "${LOG_DIR}/${cohort}.log" 2> "${LOG_DIR}/${cohort}.err"
    ) &
    
    # Store PID
    PID=$!
    PIDS+=($PID)
    COHORT_PIDS[$PID]=$cohort
    COHORT_STATUS[$cohort]="running"
    COHORT_START_TIME[$cohort]=$(date +%s)
    
    echo "  PID: $PID"
    echo "  Environment: ${COHORT_ENV}"
    echo "  Log: ${LOG_DIR}/${cohort}.log"
    
    sleep 2
done

echo ""
echo "All cohorts launched!"
update_status

# Monitor progress
echo ""
echo "Monitoring progress (updates every 60 seconds)..."
echo "Press Ctrl+C to stop monitoring (pipelines will continue)"
echo ""

MONITORING=true
trap 'MONITORING=false' INT

while $MONITORING; do
    sleep 60
    
    # Check status of all cohorts
    ALL_DONE=true
    for PID in "${PIDS[@]}"; do
        cohort="${COHORT_PIDS[$PID]}"
        
        if kill -0 $PID 2>/dev/null; then
            # Still running
            COHORT_STATUS[$cohort]="running"
            ALL_DONE=false
            
            # Check for errors
            if grep -q "ERROR\|FAILED" "${LOG_DIR}/${cohort}.err" 2>/dev/null; then
                echo ""
                echo -e "${RED}⚠ ${cohort^^} has errors (but still running)${NC}"
                tail -3 "${LOG_DIR}/${cohort}.err" | sed 's/^/  /'
            fi
        else
            # Process finished
            if wait $PID 2>/dev/null; then
                COHORT_STATUS[$cohort]="completed"
            else
                COHORT_STATUS[$cohort]="failed"
                echo ""
                echo -e "${RED}✗ ${cohort^^} FAILED${NC}"
                if [ -f "${LOG_DIR}/${cohort}.err" ]; then
                    echo "Last errors:"
                    tail -5 "${LOG_DIR}/${cohort}.err" | sed 's/^/  /'
                fi
            fi
        fi
    done
    
    update_status
    
    if $ALL_DONE; then
        break
    fi
done

# Final status
echo ""
echo "All cohorts finished!"
echo ""

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

# Run meta-analysis if enough cohorts succeeded
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
        echo "║  Cohorts: ${COMPLETED}/8 (truly parallel, isolated)                       ║"
        echo "║  Meta-analysis: COMPLETE                                       ║"
        echo "║                                                                ║"
        echo "║  Total time: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m                                           ║"
        echo "║                                                                ║"
        echo "║  Results:                                                      ║"
        echo "║    Individual: results/{COHORT}/regenie_step2/                ║"
        echo "║    Meta-analysis: results/meta_analysis/                      ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        echo ""
        
        echo "Meta-analysis results:"
        ls -lh results/meta_analysis/*.tbl 2>/dev/null | wc -l | xargs echo "  Files:"
        echo ""
        
    else
        echo ""
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║                ⚠ META-ANALYSIS FAILED ⚠                       ║"
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

cat > "${LOG_DIR}/SUMMARY.txt" << EOF
TRULY PARALLEL EXECUTION SUMMARY
=================================

Date: $(date)
Mode: Truly Parallel with Complete Isolation

Duration: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m

Results:
  Successful: ${COMPLETED}/8
  Failed: ${FAILED}/8

Successful: ${SUCCESS_COHORTS[@]:-none}
Failed: ${FAILED_COHORTS[@]:-none}

Isolation Strategy:
  Each cohort ran with:
    - Separate NXF_HOME (cache): ${ISOLATED_DIR}/{cohort}/cache
    - Separate work dir: ${ISOLATED_DIR}/{cohort}/work
    - Separate temp dir: ${ISOLATED_DIR}/{cohort}/temp
  
  This prevents ALL Nextflow conflicts!

Logs: ${LOG_DIR}
Results: results/meta_analysis/
EOF

echo "Summary: ${LOG_DIR}/SUMMARY.txt"
echo ""

exit 0

