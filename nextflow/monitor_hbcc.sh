#!/bin/bash
# Monitor HBCC h650 pipeline status
# Usage: ./monitor_hbcc.sh

set -e

STUDY="NIMH_HBCC_h650"
WORK_DIR="${1:-work/NIMH_HBCC_h650}"
RESULT_DIR="results/${STUDY}"
LOG_DIR="logs"

echo "=============================================="
echo "  HBCC h650 Pipeline Monitor"
echo "  $(date)"
echo "=============================================="
echo ""

# 1. SLURM jobs
echo "--- SLURM jobs (user: $(whoami)) ---"
squeue -u $(whoami) 2>/dev/null || echo "(squeue not available)"
echo ""

# 2. Nextflow processes (if nextflow run in background)
echo "--- Nextflow/Java processes ---"
ps aux | grep -E "nextflow|Nextflow" | grep -v grep || true
echo ""

# 3. Work dir status
echo "--- Work dir: ${WORK_DIR} ---"
if [ -d "$WORK_DIR" ]; then
  echo "  Genotyping QC artifacts:"
  ls -la "$WORK_DIR"/CMC_HBCC.QC.*.pgen 2>/dev/null | wc -l | xargs echo "    .pgen files:"
  ls -la "$WORK_DIR"/CMC_HBCC.QC.merged.pgen 2>/dev/null || echo "    merged.pgen: not yet"
  echo "  Task dirs:"
  find "$WORK_DIR" -maxdepth 2 -type d -name "[0-9a-f]*" 2>/dev/null | wc -l | xargs echo "    count:"
else
  echo "  (not found)"
fi
echo ""

# 4. Results dir - key outputs
echo "--- Results: ${RESULT_DIR} ---"
if [ -d "$RESULT_DIR" ]; then
  for f in phenotypes_RINT.txt samples_with_phenotypes.txt clinical_covariates.txt covariates.txt CMC_HBCC.QC.final.pgen pca.csv; do
    if [ -f "$RESULT_DIR/$f" ]; then
      echo "  [OK] $f"
    else
      echo "  [--] $f (missing)"
    fi
  done
else
  echo "  (not found)"
fi
echo ""

# 5. Latest log
echo "--- Latest log ---"
LATEST=$(ls -t "$LOG_DIR"/hbcc_h650*.log 2>/dev/null | head -1)
if [ -n "$LATEST" ]; then
  echo "  $LATEST"
  tail -20 "$LATEST"
else
  echo "  No hbcc_h650 logs found"
fi
