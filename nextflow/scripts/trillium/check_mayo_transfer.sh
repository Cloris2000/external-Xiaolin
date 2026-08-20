#!/usr/bin/env bash
# Check whether the Mayo transfer (run inside tmux) has finished.
#
# Run on Trillium (tri-dm2), same node where you started the transfer:
#   bash check_mayo_transfer.sh
#
# It reports:
#   1. Whether the tmux session is still alive
#   2. Whether an rsync process is still running
#   3. The last lines of the newest transfer log
#   4. Whether the "transfer finished" marker was printed
#   5. Destination size + file counts

set -uo pipefail

SESSION="${SESSION:-mayo_xfer}"
DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/Mayo}"
LOG_DIR="${LOG_DIR:-${DEST_ROOT}/transfer_logs}"

echo "=================================================="
echo "Mayo transfer status check  ($(date '+%F %T'))"
echo "  SESSION  = ${SESSION}"
echo "  DEST     = ${DEST_ROOT}"
echo "=================================================="
echo

# 1) tmux session
echo "--- tmux session ---"
if command -v tmux >/dev/null 2>&1 && tmux has-session -t "${SESSION}" 2>/dev/null; then
  echo "STATUS: session '${SESSION}' is STILL ALIVE (transfer may be running or you left the shell open)"
  SESSION_ALIVE=1
else
  echo "STATUS: session '${SESSION}' not found (it ended, or was never named that)"
  SESSION_ALIVE=0
fi
echo "Other tmux sessions:"
tmux ls 2>/dev/null || echo "  (none)"
echo

# 2) running rsync processes pulling from CAMH
echo "--- rsync processes ---"
RSYNC_PROCS=$(pgrep -af 'rsync.*192\.197\.205\.74' 2>/dev/null)
if [[ -n "${RSYNC_PROCS}" ]]; then
  echo "STATUS: rsync STILL RUNNING:"
  echo "${RSYNC_PROCS}"
  RSYNC_RUNNING=1
else
  echo "STATUS: no active rsync from 192.197.205.74"
  RSYNC_RUNNING=0
fi
echo

# 3) newest log tail
echo "--- newest transfer log ---"
if [[ -d "${LOG_DIR}" ]]; then
  LATEST_LOG=$(ls -1t "${LOG_DIR}"/xfer_mayo_*.log 2>/dev/null | head -1)
  if [[ -n "${LATEST_LOG}" ]]; then
    echo "LOG: ${LATEST_LOG}"
    tail -n 15 "${LATEST_LOG}"
  else
    echo "No xfer_mayo_*.log yet in ${LOG_DIR}"
  fi
else
  echo "Log dir ${LOG_DIR} does not exist yet"
fi
echo

# 4) completion marker across all mayo logs
echo "--- completion marker ---"
if [[ -d "${LOG_DIR}" ]] && grep -rqs "Mayo transfer finished" "${LOG_DIR}" 2>/dev/null; then
  echo "FOUND: 'Mayo transfer finished' in logs"
  FINISHED_MARKER=1
else
  # rsync per-step logs may not contain the script's final echo; also check the
  # 4 expected steps each show a completed file list without errors.
  echo "No explicit finish marker in ${LOG_DIR}"
  echo "(The final 'Mayo transfer finished' line prints to the tmux pane, not the per-step logs.)"
  FINISHED_MARKER=0
fi
echo

# 5) destination contents
echo "--- destination contents ---"
if [[ -d "${DEST_ROOT}" ]]; then
  for d in \
    "RNA/Mayo_raw_counts_Nov_18_ensembl.csv" \
    "Metadata/Mayo_meta_tissue_counts_Nov_18.csv" \
    "Genotype/Mayo_joint_WGS_vcf" \
    "Genotype/Mayo_joint_WGS_vcf_normalized"; do
    p="${DEST_ROOT}/${d}"
    if [[ -e "${p}" ]]; then
      if [[ -d "${p}" ]]; then
        n=$(find "${p}" -type f 2>/dev/null | wc -l)
        sz=$(du -sh "${p}" 2>/dev/null | cut -f1)
        echo "  OK   ${d}  (${n} files, ${sz})"
      else
        sz=$(du -sh "${p}" 2>/dev/null | cut -f1)
        echo "  OK   ${d}  (${sz})"
      fi
    else
      echo "  MISS ${d}"
    fi
  done
  echo
  echo "Total ${DEST_ROOT}: $(du -sh "${DEST_ROOT}" 2>/dev/null | cut -f1)"
  # Expected: raw ~107G + normalized ~110G => ~217G genotypes
  echo "Expected genotypes ~217G (raw ~107G + normalized ~110G)."
else
  echo "Destination ${DEST_ROOT} does not exist yet"
fi
echo

# Verdict
echo "=================================================="
if [[ "${RSYNC_RUNNING}" -eq 1 ]]; then
  echo "VERDICT: IN PROGRESS — rsync is still running."
elif [[ "${SESSION_ALIVE}" -eq 1 ]]; then
  echo "VERDICT: LIKELY DONE or IDLE — session open but no rsync running."
  echo "         Attach to confirm:  tmux attach -t ${SESSION}"
else
  echo "VERDICT: FINISHED (no session, no rsync). Verify sizes above match expectations."
fi
echo "=================================================="
