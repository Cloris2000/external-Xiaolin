#!/bin/bash
#
# Monitor GTEx pipeline progress.
# Usage: ./monitor_gtex_run.sh [LOG_DIR]
#   LOG_DIR defaults to the most recent logs/truly_parallel_* directory.
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

if [[ -n "${1:-}" ]]; then
  LOG_DIR="$1"
else
  LOG_DIR=$(ls -td logs/truly_parallel_* 2>/dev/null | head -1)
fi

if [[ -z "$LOG_DIR" || ! -d "$LOG_DIR" ]]; then
  echo "No log directory found. Specify: $0 <path-to-logs/truly_parallel_YYYYMMDD_HHMMSS>"
  exit 1
fi

GTEX_LOG="${LOG_DIR}/gtex.log"
GTEX_ERR="${LOG_DIR}/gtex.err"

echo "=============================================="
echo "GTEx pipeline monitor"
echo "Log dir: $LOG_DIR"
echo "Main log: $GTEX_LOG"
echo "=============================================="
echo ""

# Is Nextflow / run still active?
if pgrep -f "run_8_cohorts_truly_parallel" >/dev/null 2>&1 || pgrep -f "nextflow.*combined_pipeline_v2" >/dev/null 2>&1; then
  echo "Status: RUNNING"
else
  echo "Status: run script / Nextflow not detected (may have finished or been run elsewhere)"
fi
echo ""

# Last N lines of main log
echo "--- Last 40 lines of gtex.log ---"
tail -40 "$GTEX_LOG" 2>/dev/null || echo "(no log yet)"
echo ""
echo "--- Last 20 lines of gtex.err ---"
tail -20 "$GTEX_ERR" 2>/dev/null || echo "(no err yet)"
echo ""
echo "To follow live: tail -f $GTEX_LOG"
