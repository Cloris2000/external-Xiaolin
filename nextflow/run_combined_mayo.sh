#!/bin/bash
# Launch Mayo combined pipeline with proper settings to avoid NFS issues

# Load environment (PATH and Java)
export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$HOME/.local/bin:$JAVA_HOME/bin:$PATH"

# Use shared NFS directory for work (accessible from all compute nodes)
# With concurrency limits (maxForks=10, queueSize=20) to avoid NFS overload
WORK_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/work_combined_Mayo"

# Clean up old work directory if it exists
# rm -rf "$WORK_DIR" 2>/dev/null

echo "================================================"
echo "Running Mayo Combined Pipeline"
echo "================================================"
echo "Work directory (temporary): $WORK_DIR"
echo "Output directory (permanent): results/Mayo/"
echo "================================================"
echo ""

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

/nethome/kcni/xzhou/.local/bin/nextflow run combined_pipeline.nf \
    -c nextflow.config.combined.mayo \
    -w "$WORK_DIR" \
    -resume \
    "$@"

