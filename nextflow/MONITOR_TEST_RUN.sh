#!/bin/bash
# Monitor script for SynapseID test run
# Run ID: test_synapseid_20260130_163359

LOG_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/test_synapseid_20260130_163359"

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║         SynapseID Test Run - Monitoring Dashboard              ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Function to check sample counts
check_samples() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Sample Count Validation (KEY METRIC)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    for cohort in rosmap cmc_mssm cmc_penn cmc_pitt; do
        echo ""
        echo ">>> $cohort <<<"
        grep -E "Samples with both phenotype|Final sample count:" "${LOG_DIR}/${cohort}.log" 2>/dev/null | tail -4 || echo "  Waiting for PHENO_PREP..."
    done
    echo ""
}

# Function to check pipeline status
check_status() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Pipeline Status"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    for cohort in rosmap cmc_mssm cmc_penn cmc_pitt; do
        echo ""
        echo ">>> $cohort <<<"
        tail -3 "${LOG_DIR}/${cohort}.log" 2>/dev/null || echo "  Log not available yet"
    done
    echo ""
}

# Function to check for errors
check_errors() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Error Check"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    for cohort in rosmap cmc_mssm cmc_penn cmc_pitt; do
        errors=$(grep -i "ERROR\|Failed\|terminated with" "${LOG_DIR}/${cohort}.log" 2>/dev/null | grep -v "errorStrategy" | tail -2)
        if [ -n "$errors" ]; then
            echo ""
            echo "⚠️  $cohort has errors:"
            echo "$errors"
        fi
    done
    echo ""
}

# Function to check completion
check_completion() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Completion Status"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    for cohort in rosmap cmc_mssm cmc_penn cmc_pitt; do
        completion=$(grep -E "Completed at:|Duration.*:|Succeeded.*:" "${LOG_DIR}/${cohort}.log" 2>/dev/null | tail -3)
        if [ -n "$completion" ]; then
            echo ""
            echo "✅ $cohort COMPLETED:"
            echo "$completion"
        fi
    done
    echo ""
}

# Main monitoring loop
case "${1:-status}" in
    samples)
        check_samples
        ;;
    errors)
        check_errors
        ;;
    completion)
        check_completion
        ;;
    all)
        check_samples
        check_status
        check_errors
        check_completion
        ;;
    *)
        check_status
        ;;
esac

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║ Usage:                                                         ║"
echo "║   bash MONITOR_TEST_RUN.sh [samples|errors|completion|all]    ║"
echo "║                                                                ║"
echo "║ Or watch in real-time:                                        ║"
echo "║   watch -n 30 'bash MONITOR_TEST_RUN.sh samples'              ║"
echo "╚════════════════════════════════════════════════════════════════╝"
