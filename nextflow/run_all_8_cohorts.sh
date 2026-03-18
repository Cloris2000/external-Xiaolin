#!/bin/bash
#
# Run All 8 Cohorts + Automatic Meta-Analysis
#
# This script runs the combined pipeline for all 8 cohorts.
# Meta-analysis automatically runs AFTER all cohorts complete!
#
# How it works:
#   - gwas_pipeline.nf uses groupTuple(by: 0) which WAITS for all cohorts
#   - Meta-analysis only starts after ALL cohorts finish GWAS
#   - No manual synchronization needed!
# * Combined Pipeline for Multiple Cohorts with Final Meta-Analysis
# /* 
# * This pipeline:
# * 1. Processes ALL cohorts in parallel through the complete pipeline
# * 2. Collects GWAS results from ALL cohorts
# * 3. Runs meta-analysis ONCE after ALL cohorts complete
# * 
# * Execution Flow:
# * 
# * FOR EACH COHORT (in parallel):
# *   ├─ RNA-seq Processing
# *   ├─ Phenotype Prep
# *   ├─ Genotype QC  
# *   ├─ Covariate Prep
# *   └─ GWAS (Regenie Step 1 & 2)
# * 
# * THEN (after all cohorts complete):
# *   └─ Meta-Analysis (across all cohorts)
# * 
# * Cohorts: ROSMAP, Mayo, MSBB, CMC_MSSM, CMC_PENN, CMC_PITT, GTEx, NABEC
# */
# See META_ANALYSIS_TIMING.md for technical details.
#

set -e
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Configuration
CONFIG_FILE="nextflow.config.combined.all_cohorts"
LOG_DIR="${SCRIPT_DIR}/logs/run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║            Run All 8 Cohorts + Meta-Analysis                   ║"
echo "╠════════════════════════════════════════════════════════════════╣"
echo "║  Pipeline: combined_pipeline_v2.nf                             ║"
echo "║  Config: ${CONFIG_FILE}                    ║"
echo "║                                                                ║"
echo "║  Cohorts (8):                                                  ║"
echo "║    ROSMAP, Mayo, MSBB, CMC_MSSM, CMC_PENN,                     ║"
echo "║    CMC_PITT, GTEx, NABEC                                       ║"
echo "║                                                                ║"
echo "║  Pipeline Flow:                                                ║"
echo "║    1. Process all cohorts (in parallel)                        ║"
echo "║    2. Run GWAS for each cohort (19 cell types)                 ║"
echo "║    3. Meta-analysis AUTOMATICALLY waits for ALL cohorts        ║"
echo "║                                                                ║"
echo "║  Technical Detail:                                             ║"
echo "║    groupTuple() in gwas_pipeline.nf ensures meta-analysis      ║"
echo "║    ONLY starts after ALL cohorts complete!                     ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Check config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: ${CONFIG_FILE}"
    echo ""
    echo "Please ensure the file exists and contains all 8 cohorts."
    exit 1
fi

# Display cohort list from config
echo "Cohorts configured:"
grep "cohorts =" "$CONFIG_FILE" | head -1
echo ""

echo "Log directory: ${LOG_DIR}"
echo "Main log: ${LOG_DIR}/pipeline.log"
echo "Errors: ${LOG_DIR}/pipeline.err"
echo ""
echo "Starting pipeline at $(date)..."
echo "This will take several hours (4-16h depending on parallelization)."
echo ""
echo "Monitor progress:"
echo "  tail -f ${LOG_DIR}/pipeline.log"
echo "  tail -f .nextflow.log"
echo ""

START_TIME=$(date +%s)

# Run the pipeline
if nextflow run combined_pipeline_v2.nf \
    -c "$CONFIG_FILE" \
    -resume \
    > "${LOG_DIR}/pipeline.log" 2> "${LOG_DIR}/pipeline.err"; then
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                        ✓ SUCCESS ✓                            ║"
    echo "╠════════════════════════════════════════════════════════════════╣"
    echo "║  All 8 cohorts processed successfully!                        ║"
    echo "║  Meta-analysis completed for 19 cell types!                   ║"
    echo "║                                                                ║"
    echo "║  Duration: ${HOURS}h ${MINUTES}m                                             ║"
    echo "║                                                                ║"
    echo "║  Results:                                                      ║"
    echo "║    Individual cohorts: results/{COHORT}/regenie_step2/       ║"
    echo "║    Meta-analysis:      results/meta_analysis/                 ║"
    echo "║                                                                ║"
    echo "║  Meta-analysis files (19 cell types):                         ║"
    echo "║    results/meta_analysis/Astrocyte_meta1.tbl                  ║"
    echo "║    results/meta_analysis/Microglia_meta1.tbl                  ║"
    echo "║    results/meta_analysis/SST_meta1.tbl                        ║"
    echo "║    ... (and 16 more)                                           ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    
    # List meta-analysis results
    echo "Meta-analysis results:"
    ls -lh results/meta_analysis/*.tbl 2>/dev/null | awk '{print "  " $9, "(" $5 ")"}'
    echo ""
    
    # Save summary
    cat > "${LOG_DIR}/summary.txt" << EOF
EXECUTION SUMMARY
=================

Status: SUCCESS ✓
Start: $(date -d @$START_TIME '+%Y-%m-%d %H:%M:%S')
End: $(date -d @$END_TIME '+%Y-%m-%d %H:%M:%S')
Duration: ${HOURS}h ${MINUTES}m

Pipeline: combined_pipeline_v2.nf
Config: ${CONFIG_FILE}

Cohorts Processed: 8
  - ROSMAP
  - Mayo
  - MSBB
  - CMC_MSSM
  - CMC_PENN
  - CMC_PITT
  - GTEx
  - NABEC

Cell Types: 19
Meta-analysis: Completed

Results:
  - Individual GWAS: results/{COHORT}/regenie_step2/
  - Meta-analysis: results/meta_analysis/

Technical Note:
  Meta-analysis automatically waited for ALL cohorts to complete
  before running (via groupTuple in gwas_pipeline.nf).
  
Logs:
  - Pipeline: ${LOG_DIR}/pipeline.log
  - Errors: ${LOG_DIR}/pipeline.err
  - Nextflow: .nextflow.log
EOF
    
    echo "Summary saved: ${LOG_DIR}/summary.txt"
    echo ""
    
else
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                        ✗ FAILED ✗                             ║"
    echo "╠════════════════════════════════════════════════════════════════╣"
    echo "║  Pipeline execution failed!                                    ║"
    echo "║                                                                ║"
    echo "║  Duration before failure: ${HOURS}h ${MINUTES}m                           ║"
    echo "║                                                                ║"
    echo "║  Check error logs:                                             ║"
    echo "║    ${LOG_DIR}/pipeline.err              ║"
    echo "║    .nextflow.log                                               ║"
    echo "║                                                                ║"
    echo "║  To resume after fixing:                                       ║"
    echo "║    nextflow run combined_pipeline_v2.nf \\                     ║"
    echo "║      -c ${CONFIG_FILE} \\                   ║"
    echo "║      -resume                                                   ║"
    echo "║                                                                ║"
    echo "║  The -resume flag will:                                        ║"
    echo "║    - Skip cohorts that already completed                       ║"
    echo "║    - Re-run only failed parts                                  ║"
    echo "║    - Still wait for ALL cohorts before meta-analysis           ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Show last few error lines
    if [ -s "${LOG_DIR}/pipeline.err" ]; then
        echo "Last 20 error lines:"
        tail -20 "${LOG_DIR}/pipeline.err"
        echo ""
    fi
    
    exit 1
fi

exit 0

