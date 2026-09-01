#!/usr/bin/env bash
# Shared helpers for Trillium <- CAMH rsync transfers.
# Source from cohort transfer scripts (do not run directly).
#
# CAMH @ 192.197.205.74 is rssh-only (rsync/scp/sftp). Gateway sees
# netdata_kcni; paths under external_data / nethome must be staged first
# via stage_remaining_cohorts_on_camh.sh.

: "${SRC_HOST:=xzhou@192.197.205.74}"
: "${SSH_KEY:=${HOME}/.ssh/id_ed25519_camh}"
: "${STAGE_ROOT:=/external/rprshnas01/netdata_kcni/stlab/DELETE_ME/Xiaolin/nextflow/data_input/trillium_staging}"

# Location of this library / the SSH_ASKPASS helper alongside it.
LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${CAMH_ASKPASS:=${LIB_DIR}/camh_askpass.sh}"

RSYNC_OPTS=(-avP)
DRY_RUN=0
# IMPORTANT: The CAMH data-mover (192.197.205.74) is extremely picky — the ONLY
# invocation observed to work is a BARE `ssh` with no custom options (matching a
# plain `rsync -avP host:path dest`). Adding -o options (ControlMaster, key,
# ServerAlive, PreferredAuthentications, etc.) makes it close the session right
# after the password ("Connection closed by ... port 22"). So we pass NO ssh
# options by default. The non-default key filename is never auto-offered and
# ControlMaster is off by default, so bare ssh already avoids those pitfalls.
# You will be prompted per rsync step, UNLESS you provide a password file for
# sshpass (unattended):
#   printf '%s' 'YOUR_CAMH_PASSWORD' > ~/.camh_pass && chmod 600 ~/.camh_pass
#   export CAMH_PASSFILE=~/.camh_pass
# Extra ssh options can be injected via EXTRA_SSH_OPTS if ever needed.
SSH_OPTS=()
if [[ -n "${EXTRA_SSH_OPTS:-}" ]]; then
  # shellcheck disable=SC2206
  SSH_OPTS=(${EXTRA_SSH_OPTS})
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
  # Optional unattended auth using CAMH_PASSFILE:
  #   1) sshpass if available; else
  #   2) SSH_ASKPASS (OpenSSH >= 8.4) — feeds the password without a pty,
  #      keeping rsync's stdin/stdout clean (safer than an expect pty wrapper).
  local rsh="ssh ${SSH_OPTS[*]}"
  if [[ -n "${CAMH_PASSFILE:-}" ]]; then
    if command -v sshpass >/dev/null 2>&1; then
      rsh="sshpass -f ${CAMH_PASSFILE} ssh ${SSH_OPTS[*]}"
    elif [[ -x "${CAMH_ASKPASS}" ]]; then
      export CAMH_PASSFILE
      export SSH_ASKPASS="${CAMH_ASKPASS}"
      export SSH_ASKPASS_REQUIRE="force"
      export DISPLAY="${DISPLAY:-:0}"
    else
      echo "WARNING: CAMH_PASSFILE set but neither sshpass nor ${CAMH_ASKPASS} usable — interactive prompts" >&2
    fi
  fi
  RSYNC_OPTS+=(-e "${rsh}")
}

trillium_init() {
  local dest_root="$1"
  LOG_DIR="${LOG_DIR:-${dest_root}/transfer_logs}"
  mkdir -p "${dest_root}" "${LOG_DIR}"
  echo "Transfer starting (dry_run=${DRY_RUN})"
  echo "DEST_ROOT=${dest_root}"
  echo "SRC_HOST=${SRC_HOST}"
  if [[ -n "${CAMH_PASSFILE:-}" ]] && command -v sshpass >/dev/null 2>&1; then
    echo "AUTH=password via sshpass (CAMH_PASSFILE=${CAMH_PASSFILE}) — unattended"
  elif [[ -n "${CAMH_PASSFILE:-}" && -n "${SSH_ASKPASS:-}" ]]; then
    echo "AUTH=password via SSH_ASKPASS (CAMH_PASSFILE=${CAMH_PASSFILE}) — unattended"
  else
    echo "AUTH=password (key & ControlMaster disabled — CAMH mover rejects both)"
    echo "You will be prompted for the CAMH password once per step."
  fi
  echo
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
