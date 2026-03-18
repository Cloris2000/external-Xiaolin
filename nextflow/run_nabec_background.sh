#!/bin/bash
# Run NABEC Combined Pipeline in background
# This script ensures NABEC runs completely in the background

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

echo "================================================"
echo "Starting NABEC Combined Pipeline"
echo "================================================"
echo "Started at: $(date)"
echo "Config: nextflow.config.combined.nabec"
echo "Work directory: work/NABEC"
echo "Output directory: results/NABEC/"
echo "Log file: pipeline_run_nabec.log"
echo "================================================"

# Run with nohup to ensure it continues after logout
nohup nextflow run combined_pipeline.nf \
    -c nextflow.config.combined.nabec \
    -resume \
    > pipeline_run_nabec.log 2>&1 &

# Get the process ID
NABEC_PID=$!

echo "NABEC pipeline started with PID: $NABEC_PID"
echo "To monitor progress, run:"
echo "  tail -f pipeline_run_nabec.log"
echo ""
echo "To check status:"
echo "  ps aux | grep $NABEC_PID"
echo ""
echo "================================================"
