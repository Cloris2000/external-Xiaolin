#!/bin/bash
# Launch ROSMAP_array combined pipeline
# Genotype source: merged_overlap_rs BED split to per-chr VCFs (170 TOPmed array samples)

export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$HOME/.local/bin:$JAVA_HOME/bin:$PATH"

WORK_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/work_combined_ROSMAP_array"

echo "================================================"
echo "Running ROSMAP_array Combined Pipeline"
echo "================================================"
echo "Work directory (temporary): $WORK_DIR"
echo "Output directory (permanent): results/ROSMAP_array/"
echo "================================================"
echo ""

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

/nethome/kcni/xzhou/.local/bin/nextflow run combined_pipeline_v2.nf \
    -c nextflow.config.combined.rosmap_array \
    -w "$WORK_DIR" \
    -resume \
    "$@"

exit_code=$?

echo ""
echo "================================================"
echo "Pipeline finished with exit code: $exit_code"
echo "================================================"

exit $exit_code
