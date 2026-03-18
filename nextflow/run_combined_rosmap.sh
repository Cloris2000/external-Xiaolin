#!/bin/bash
# Launch ROSMAP combined pipeline with proper settings to avoid NFS issues

# Load environment (PATH and Java)
export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$HOME/.local/bin:$JAVA_HOME/bin:$PATH"

# Use shared NFS directory for work (accessible from all compute nodes)
# With concurrency limits (maxForks=10, queueSize=20) to avoid NFS overload
WORK_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/work_combined_ROSMAP"

# Clean up old work directory if it exists
# rm -rf "$WORK_DIR" 2>/dev/null

echo "================================================"
echo "Running ROSMAP Combined Pipeline"
echo "================================================"
echo "Work directory (temporary): $WORK_DIR"
echo "Output directory (permanent): results/ROSMAP/"
echo "================================================"
echo ""

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

/nethome/kcni/xzhou/.local/bin/nextflow run combined_pipeline.nf \
    -c nextflow.config.combined.rosmap \
    -w "$WORK_DIR" \
    -resume \
    "$@"

exit_code=$?

echo ""
echo "================================================"
echo "Pipeline finished with exit code: $exit_code"
echo "================================================"

# Optional: Clean up work directory after successful completion
# if [ $exit_code -eq 0 ]; then
#     echo "Cleaning up work directory..."
#     rm -rf "$WORK_DIR"
# fi

exit $exit_code
