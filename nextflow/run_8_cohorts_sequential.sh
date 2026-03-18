#!/bin/bash
#
# Run All 8 Cohorts Sequentially, Then Meta-Analysis
#
# This script:
# 1. Runs combined_pipeline_v2.nf for EACH cohort (one at a time or in parallel)
# 2. After ALL cohorts complete, runs standalone meta-analysis
#
# Note: combined_pipeline_v2.nf is designed for ONE cohort at a time.
#       For multiple cohorts, we run it 8 times.
#

set -e
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Cohorts to process
COHORTS=("rosmap" "mayo" "msbb" "cmc_mssm" "cmc_penn" "cmc_pitt" "gtex" "nabec")

# Log directory
LOG_DIR="${SCRIPT_DIR}/logs/multi_cohort_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║        Run All 8 Cohorts + Meta-Analysis                       ║"
echo "╠════════════════════════════════════════════════════════════════╣"
echo "║  Approach: Run each cohort separately, then meta-analysis      ║"
echo "║                                                                ║"
echo "║  Cohorts (8): ROSMAP, Mayo, MSBB, CMC_MSSM, CMC_PENN,          ║"
echo "║               CMC_PITT, GTEx, NABEC                            ║"
echo "║                                                                ║"
echo "║  Pipeline: combined_pipeline_v2.nf (per cohort)                ║"
echo "║  Meta-analysis: standalone_meta_analysis.nf (at the end)       ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

TOTAL_START=$(date +%s)
SUCCESS_COUNT=0
FAILED_COUNT=0

# Run each cohort
for cohort in "${COHORTS[@]}"; do
    CONFIG="nextflow.config.combined.${cohort}"
    
    echo ""
    echo "=========================================="
    echo "Processing Cohort: ${cohort^^}"
    echo "=========================================="
    echo "Config: ${CONFIG}"
    echo "Log: ${LOG_DIR}/${cohort}.log"
    
    if [ ! -f "$CONFIG" ]; then
        echo "⚠ WARNING: Config not found: ${CONFIG}"
        echo "Skipping ${cohort}..."
        ((FAILED_COUNT++))
        continue
    fi
    
    START_TIME=$(date +%s)
    
    # Run pipeline for this cohort
    if nextflow run combined_pipeline_v2.nf \
        -c "$CONFIG" \
        -resume \
        > "${LOG_DIR}/${cohort}.log" 2> "${LOG_DIR}/${cohort}.err"; then
        
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))
        MINUTES=$((DURATION / 60))
        
        echo "✓ ${cohort^^} completed in ${MINUTES} minutes"
        ((SUCCESS_COUNT++))
    else
        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))
        MINUTES=$((DURATION / 60))
        
        echo "✗ ${cohort^^} FAILED after ${MINUTES} minutes"
        echo "  Check error log: ${LOG_DIR}/${cohort}.err"
        ((FAILED_COUNT++))
    fi
done

echo ""
echo "=========================================="
echo "Cohort Processing Summary"
echo "=========================================="
echo "Successful: ${SUCCESS_COUNT}/8"
echo "Failed: ${FAILED_COUNT}/8"
echo "=========================================="
echo ""

# Only run meta-analysis if at least 2 cohorts succeeded
if [ $SUCCESS_COUNT -ge 2 ]; then
    echo ""
    echo "=========================================="
    echo "Running Meta-Analysis"
    echo "=========================================="
    echo "Combining results from ${SUCCESS_COUNT} cohorts..."
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
        echo "║  Cohorts processed: ${SUCCESS_COUNT}/8                                       ║"
        echo "║  Meta-analysis: COMPLETE                                       ║"
        echo "║                                                                ║"
        echo "║  Total duration: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m                                       ║"
        echo "║  Meta-analysis: ${META_MINUTES}m                                           ║"
        echo "║                                                                ║"
        echo "║  Results:                                                      ║"
        echo "║    - Individual cohorts: results/{COHORT}/regenie_step2/     ║"
        echo "║    - Meta-analysis: results/meta_analysis/                    ║"
        echo "╚════════════════════════════════════════════════════════════════╝"
        echo ""
        
        # List meta-analysis results
        echo "Meta-analysis results:"
        ls -lh results/meta_analysis/*.tbl 2>/dev/null | head -10
        echo ""
        
    else
        echo ""
        echo "╔════════════════════════════════════════════════════════════════╗"
        echo "║                ⚠ META-ANALYSIS FAILED ⚠                       ║"
        echo "╠════════════════════════════════════════════════════════════════╣"
        echo "║  Cohorts processed: ${SUCCESS_COUNT}/8                                       ║"
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
    echo "║  Only ${SUCCESS_COUNT} cohort(s) succeeded.                                ║"
    echo "║  Need at least 2 cohorts for meta-analysis.                   ║"
    echo "║  Meta-analysis SKIPPED.                                        ║"
    echo "║                                                                ║"
    echo "║  Fix failed cohorts and re-run with:                          ║"
    echo "║    bash run_8_cohorts_sequential.sh                           ║"
    echo "║                                                                ║"
    echo "║  (Failed cohorts will be retried, successful ones skipped     ║"
    echo "║   thanks to -resume flag)                                     ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    exit 1
fi

# Save summary
cat > "${LOG_DIR}/SUMMARY.txt" << EOF
MULTI-COHORT EXECUTION SUMMARY
==============================

Date: $(date)
Total Duration: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m

Cohorts Processed: ${SUCCESS_COUNT}/8
Failed Cohorts: ${FAILED_COUNT}/8

Pipeline: combined_pipeline_v2.nf (per cohort)
Meta-analysis: standalone_meta_analysis.nf

Results:
  - Individual cohorts: results/{COHORT}/regenie_step2/
  - Meta-analysis: results/meta_analysis/

Logs Directory: ${LOG_DIR}
  - Per-cohort logs: ${LOG_DIR}/{cohort}.log
  - Meta-analysis log: ${LOG_DIR}/meta_analysis.log

To check individual cohort results:
  ls -lh results/*/regenie_step2/*.regenie

To check meta-analysis results:
  ls -lh results/meta_analysis/*.tbl
EOF

echo "Summary saved: ${LOG_DIR}/SUMMARY.txt"
echo ""

exit 0


