#!/usr/bin/env bash
# Shared helpers for Trillium <- CAMH rsync transfers.
# Source from cohort transfer scripts (do not run directly).
#
# CAMH @ 192.197.205.74 is rssh-only (rsync/scp/sftp). Gateway sees
# netdata_kcni; paths under external_data / nethome must be staged first
# via stage_remaining_cohorts_on_camh.sh.

: "${SRC_HOST:=xzhou@192.197.205.74}"
: "${SSH_KEY:=${HOME}/.ssh/id_ed25519_camh}"
: "${STAGE_ROOT:=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data_input/trillium_staging}"

RSYNC_OPTS=(-avP)
DRY_RUN=0
CTRL_DIR="${HOME}/.ssh/sockets"
CTRL_PATH="${CTRL_DIR}/camh-xfer-%C"
SSH_OPTS=(
  -o "ControlMaster=auto"
  -o "ControlPath=${CTRL_PATH}"
  -o "ControlPersist=8h"
  -o "ServerAliveInterval=60"
  -o "ServerAliveCountMax=3"
)
if [[ -f "${SSH_KEY}" ]]; then
  SSH_OPTS+=(-i "${SSH_KEY}" -o "IdentitiesOnly=yes")
fi

trillium_parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --dry-run) DRY_RUN=1; shift ;;
      --dest) DEST_ROOT="$2"; shift 2 ;;
      --src) SRC_HOST="$2"; shift 2 ;;
      -h|--help) trillium_usage; exit 0 ;;
      *) echo "Unknown option: $1" >&2; trillium_usage; exit 1 ;;
    esac
  done
  if [[ "$DRY_RUN" -eq 1 ]]; then
    RSYNC_OPTS+=(-n)
  fi
  RSYNC_OPTS+=(-e "ssh ${SSH_OPTS[*]}")
}

trillium_init() {
  local dest_root="$1"
  LOG_DIR="${LOG_DIR:-${dest_root}/transfer_logs}"
  mkdir -p "${dest_root}" "${LOG_DIR}" "${CTRL_DIR}"
  chmod 700 "${CTRL_DIR}"
  trap trillium_cleanup_ssh EXIT
  echo "Transfer starting (dry_run=${DRY_RUN})"
  echo "DEST_ROOT=${dest_root}"
  echo "SRC_HOST=${SRC_HOST}"
  if [[ -f "${SSH_KEY}" ]]; then
    echo "SSH_KEY=${SSH_KEY}"
  else
    echo "WARNING: ${SSH_KEY} not found — password may be required"
  fi
  echo
}

trillium_cleanup_ssh() {
  ssh -O exit "${SSH_OPTS[@]}" "${SRC_HOST}" 2>/dev/null || true
}

trillium_timestamp() { date '+%Y-%m-%d_%H%M%S'; }

trillium_rsync() {
  local label="$1"
  local src_path="$2"
  local dest_path="$3"
  local log="${LOG_DIR}/xfer_${label}_$(trillium_timestamp).log"

  mkdir -p "${dest_path}"
  echo "========================================"
  echo "[$(date '+%F %T')] ${label}"
  echo "  SRC : ${SRC_HOST}:${src_path}"
  echo "  DEST: ${dest_path}"
  echo "  LOG : ${log}"
  echo "========================================"

  rsync "${RSYNC_OPTS[@]}" \
    "${SRC_HOST}:${src_path}" \
    "${dest_path}" \
    2>&1 | tee "${log}"

  echo "[$(date '+%F %T')] Done: ${label}"
  echo
}
