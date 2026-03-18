#!/bin/bash
# Force COMPLETE FRESH rerun of ALL 19 cell types
# This script removes ALL cached work and results to ensure a completely fresh run
# Use this if you want to guarantee no caching at all

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

echo "================================================"
echo "COMPLETE FRESH RERUN - All 19 Cell Types"
echo "================================================"
echo "Started at: $(date)"
echo ""
echo "WARNING: This will remove ALL cached Nextflow work"
echo "Press Ctrl+C within 5 seconds to cancel..."
sleep 5

# Backup the old results
BACKUP_DIR="results/meta_analysis_backup_FULL_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"

echo ""
echo "Creating backup of ALL results in: $BACKUP_DIR"
cp results/meta_analysis/*_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP1.tbl* "$BACKUP_DIR/" 2>/dev/null

# Remove the work directory to force complete rerun
echo "Removing Nextflow work directory..."
rm -rf work/

# Remove .nextflow cache
echo "Removing .nextflow cache..."
rm -rf .nextflow/
rm -f .nextflow.log*

# Remove all final output files
echo "Removing all final meta-analysis outputs..."
rm -f results/meta_analysis/*_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP1.tbl*
rm -f results/meta_analysis/*.done

echo ""
echo "================================================"
echo "Starting FRESH meta-analysis run..."
echo "================================================"
echo ""

# Run WITHOUT -resume to ensure fresh start
nohup nextflow run standalone_meta_analysis.nf \
    -c nextflow.config.standalone_meta \
    > pipeline_fresh_run_all.log 2>&1 &

PIPELINE_PID=$!

echo "Pipeline started with PID: $PIPELINE_PID"
echo "Log file: pipeline_fresh_run_all.log"
echo ""
echo "To monitor progress:"
echo "  tail -f pipeline_fresh_run_all.log"
echo ""
echo "To check running jobs:"
echo "  squeue -u \$USER"
echo ""
echo "Backup of old results: $BACKUP_DIR"
echo "================================================"
