#!/bin/bash
#
# Monitor standalone meta-analysis progress.
# Verifies that all 8 cohorts are used for each cell type.
#
# Usage: ./monitor_meta_analysis.sh
#

BASE_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
META_DIR="${BASE_DIR}/results/meta_analysis"
WORK_DIR="${BASE_DIR}/work"
SUFFIX="CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP"

echo "=============================================="
echo "Standalone Meta-Analysis Monitor"
echo "=============================================="
echo ""

# 1. Raw_p input check: all 8 cohorts?
echo "--- Input: raw_p files per cohort ---"
for c in ROSMAP Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx NABEC; do
  n=$(ls ${BASE_DIR}/results/$c/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
  echo "  $c: $n files"
done
echo "  Total cohorts with data: $(for c in ROSMAP Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx NABEC; do [ $(ls ${BASE_DIR}/results/$c/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l) -gt 0 ] && echo $c; done | wc -l)/8"
echo ""

# 2. Output progress: how many cell types done?
echo "--- Output: meta-analysis results ---"
if [[ -d "$META_DIR" ]]; then
  done_count=$(ls ${META_DIR}/*_meta_analysis_${SUFFIX}.tbl 2>/dev/null | wc -l)
  echo "  Completed: $done_count / 19 cell types"
  if [[ $done_count -gt 0 ]]; then
    echo "  Done: $(ls ${META_DIR}/*_meta_analysis_${SUFFIX}.tbl 2>/dev/null | xargs -I{} basename {} .tbl | sed 's/_meta_analysis.*//' | tr '\n' ' ')"
  fi
else
  echo "  (meta_analysis dir not found yet)"
fi
echo ""

# 3. Verify a running METAL task uses 8 cohorts (if work dir exists)
echo "--- Verify 8 cohorts per cell type ---"
recent_meta=$(find "$WORK_DIR" -name "*_metal_script.txt" -mmin -120 2>/dev/null | head -1)
if [[ -n "$recent_meta" ]]; then
  proc_count=$(grep -c "^PROCESS " "$recent_meta" 2>/dev/null || echo 0)
  echo "  Sample METAL script: $(dirname $recent_meta | xargs basename)"
  echo "  PROCESS lines (cohort count): $proc_count"
  if [[ $proc_count -eq 8 ]]; then
    echo "  OK: All 8 cohorts included"
  else
    echo "  Note: Expected 8, found $proc_count (may vary per task)"
  fi
else
  echo "  (No recent METAL script in work/ - run may have finished or not started)"
fi
echo ""

# 4. Nextflow / slurm status
echo "--- Run status ---"
if pgrep -f "standalone_meta_analysis.nf" >/dev/null 2>&1; then
  echo "  Nextflow: RUNNING"
  squeue -u $(whoami) 2>/dev/null | head -15 || echo "  (squeue not available)"
else
  echo "  Nextflow: not running (finished or not started)"
fi
echo ""
echo "To follow live: tail -f .nextflow.log"
echo "=============================================="
