#!/bin/bash
set -euo pipefail

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# TEST RUN: Only ROSMAP and CMC cohorts to validate synapseID fix
COHORTS=("rosmap" "cmc_mssm" "cmc_penn" "cmc_pitt")

# Log directory
LOG_DIR="${SCRIPT_DIR}/logs/test_synapseid_$(date +%Y%m%d_%H%M%S)"
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
echo "║         TEST RUN: SynapseID Fix Validation                     ║"
echo "╠════════════════════════════════════════════════════════════════╣"
echo "║  Testing synapseID support for Synapse cohorts:               ║"
echo "║    - ROSMAP                                                    ║"
echo "║    - CMC_MSSM                                                  ║"
echo "║    - CMC_PENN                                                  ║"
echo "║    - CMC_PITT                                                  ║"
echo "║                                                                ║"
echo "║  Each cohort runs in isolated Nextflow environment            ║"
echo "║  Validating: Sample ID matching via synapseID column          ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Log directory: $LOG_DIR"
echo "Isolated environments: $ISOLATED_DIR"
echo ""
echo "Launching test cohorts in parallel with isolated environments..."
echo ""

# Track PIDs
PIDS=()

# Function to launch a cohort
launch_cohort() {
    local cohort=$1
    local cohort_upper=$(echo "$cohort" | tr '[:lower:]' '[:upper:]')
    
    # Create isolated environment for this cohort
    local cohort_env="${ISOLATED_DIR}/${cohort}"
    mkdir -p "$cohort_env/launch"
    
    # Launch in background with isolated environment
    (
        cd "$cohort_env/launch"
        
        # Set isolated Nextflow environment
        export NXF_HOME="$cohort_env/.nextflow"
        export NXF_WORK="$cohort_env/work"
        export NXF_TEMP="$cohort_env/temp"
        export NXF_CACHE_DIR="$cohort_env/cache"
        
        # Create marker file
        touch .pipeline_start
        
        # Run pipeline
        nextflow run "$SCRIPT_DIR/combined_pipeline_v2.nf" \
            -c "$SCRIPT_DIR/nextflow.config.combined.${cohort}" \
            -resume \
            > "${LOG_DIR}/${cohort}.log" 2>&1
    ) &
    
    local pid=$!
    PIDS+=($pid)
    
    echo "🚀 Launching ${cohort_upper}..."
    echo "  PID: $pid"
    echo "  Environment: $cohort_env"
    echo "  Log: ${LOG_DIR}/${cohort}.log"
}

# Launch all cohorts
for cohort in "${COHORTS[@]}"; do
    launch_cohort "$cohort"
done

echo ""
echo "All test cohorts launched!"
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Monitoring tips:"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "Watch all logs:"
echo "  tail -f ${LOG_DIR}/*.log"
echo ""
echo "Check sample counts (after ~5-10 min):"
echo "  grep 'Samples with both phenotype' ${LOG_DIR}/*.log"
echo ""
echo "Monitor ROSMAP specifically:"
echo "  tail -f ${LOG_DIR}/rosmap.log"
echo ""
echo "Quick status check:"
echo "  for c in rosmap cmc_mssm cmc_penn cmc_pitt; do"
echo "    echo \"=== \$c ===\";"
echo "    tail -5 ${LOG_DIR}/\${c}.log;"
echo "  done"
echo ""
echo "════════════════════════════════════════════════════════════════"
