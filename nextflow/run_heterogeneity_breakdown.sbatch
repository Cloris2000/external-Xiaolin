#!/bin/bash
# Submit parallel heterogeneity breakdown (19 cell types) + merge + R plots.
#
# Usage:
#   bash run_heterogeneity_breakdown.sbatch
#
# Prerequisites:
#   results/meta_analysis_15cohorts/*.tbl
#   results/{cohort}/regenie_step2/*.regenie.raw_p

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/het_breakdown"
cd "${NF_DIR}"

META_DIR="${NF_DIR}/results/meta_analysis_15cohorts"
count=$(ls "${META_DIR}"/*.tbl 2>/dev/null | grep -cv '1\.tbl$' || true)
if [ "${count}" -eq 0 ]; then
  echo "ERROR: No METAL .tbl files in ${META_DIR}/"
  exit 1
fi
echo "Found ${count} METAL .tbl files."

echo "Submitting 19 parallel cell-type array jobs..."
ARRAY_JOB=$(sbatch --parsable "${NF_DIR}/run_heterogeneity_breakdown_array.sbatch")
echo "Array job: ${ARRAY_JOB}"

echo "Submitting merge job..."
MERGE_JOB=$(sbatch --parsable --dependency=afterok:${ARRAY_JOB} "${NF_DIR}/run_heterogeneity_breakdown_merge.sbatch")
echo "Merge job: ${MERGE_JOB} (depends on ${ARRAY_JOB})"

echo "Submitting R plots job..."
PLOT_JOB=$(sbatch --parsable --dependency=afterok:${MERGE_JOB} "${NF_DIR}/run_heterogeneity_breakdown_plots.sbatch")
echo "Plots job: ${PLOT_JOB} (depends on ${MERGE_JOB})"

echo ""
echo "Monitor: squeue -u \$USER | grep het"
echo "Outputs: results/meta_analysis_15cohorts/heterogeneity_breakdown/"
