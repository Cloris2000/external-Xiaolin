#!/usr/bin/env bash
# Transfer ROSMAP + ROSMAP_array pipeline inputs from CAMH to Trillium.
#
# Run on Trillium (e.g. tri-dm2 or tri-login):
#   bash transfer_rosmap.sh
#   bash transfer_rosmap.sh --dry-run
#
# Skips RNA already present under /project/rrg-shreejoy/ROSMAP/
# (ROSMAP_Raw_Counts_Bulk, ROSMAP_RNAseq_provenance).
#
# Transfers:
#   1. Shared ROSMAP biospecimen metadata (once; pipeline copy from nethome)
#   2. ROSmaster.rds + assay metadata (snpArray, wholeGenomeSeq)
#   3. ROSMAP joint WGS raw VCFs           (~426G)
#   4. ROSMAP joint WGS normalized VCFs    (~437G)
#   5. ROSMAP_array raw TOPmed imputed VCFs (~370G)
#        n1686 + n381 dose VCFs + merged/merged_overlap_rs BED
#   6. ROSMAP_array normalized VCFs        (~439M)
#   7. ROSMAP_array sample ID maps (to regenerate normalized from BED)
#
# Approx. total genotype volume: ~1.2 TB — run in tmux/screen.
#
# Auth notes:
#   CAMH @ 192.197.205.74 is rssh-only (rsync/scp/sftp — no shell).
#   Prefer SSH key setup via setup_camh_ssh_key.sh (install pubkey on CAMH).
#   ControlMaster reuses the first rsync SSH session (password once if no key).

set -euo pipefail

SRC_HOST="${SRC_HOST:-xzhou@192.197.205.74}"
DEST_ROOT="${DEST_ROOT:-/project/rrg-shreejoy/ROSMAP}"
LOG_DIR="${LOG_DIR:-${DEST_ROOT}/transfer_logs}"
# NOTE: Trillium rsync gateway sees netdata_kcni, but NOT external_data / nethome.
# Small files + TOPmed are staged on CAMH into STAGE (see stage_rosmap_on_camh.sh).
# WGS dirs already live under netdata and can be pulled directly.
WGS_BASE="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS"
STAGE="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/data_input/trillium_staging/ROSMAP"
STAGE_META="${STAGE}/Metadata"
STAGE_TOPMED="${STAGE}/Genotype/ROSMAP_array_TOPmed_imputed_vcf"
STAGE_SAMPLE_MAPS="${STAGE}/Genotype/ROSMAP_array_sample_maps"
SSH_KEY="${SSH_KEY:-${HOME}/.ssh/id_ed25519_camh}"
SKIP_TOPMED="${SKIP_TOPMED:-0}"

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

usage() {
  cat <<EOF
Usage: $(basename "$0") [--dry-run] [--dest DIR] [--src HOST]

  --dry-run     Show what would be transferred (rsync -n)
  --dest DIR    Destination root (default: ${DEST_ROOT})
  --src HOST    Source SSH host (default: ${SRC_HOST})

Environment overrides: SRC_HOST, DEST_ROOT, LOG_DIR

Tip: CAMH is rssh-only — use setup_camh_ssh_key.sh, then install the
pubkey on CAMH with a normal shell login (see that script's instructions).
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift ;;
    --dest) DEST_ROOT="$2"; shift 2 ;;
    --src) SRC_HOST="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ "$DRY_RUN" -eq 1 ]]; then
  RSYNC_OPTS+=(-n)
fi

# Reuse one SSH connection across rsyncs (first rsync may ask for password
# if the key is not yet installed on CAMH; later ones reuse ControlMaster).
# Do NOT use bare `ssh -fN` here — rssh rejects non-rsync/scp/sftp commands.
RSYNC_OPTS+=(-e "ssh ${SSH_OPTS[*]}")

mkdir -p "${DEST_ROOT}" "${LOG_DIR}" "${CTRL_DIR}" \
  "${DEST_ROOT}/Metadata" \
  "${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf" \
  "${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf_normalized" \
  "${DEST_ROOT}/Genotype/ROSMAP_array_TOPmed_imputed_vcf" \
  "${DEST_ROOT}/Genotype/ROSMAP_array_vcf_normalized" \
  "${DEST_ROOT}/Genotype/ROSMAP_array_sample_maps"
chmod 700 "${CTRL_DIR}"

cleanup_ssh() {
  # Closes local ControlMaster socket only (no remote shell needed)
  ssh -O exit "${SSH_OPTS[@]}" "${SRC_HOST}" 2>/dev/null || true
}
trap cleanup_ssh EXIT

timestamp() { date '+%Y-%m-%d_%H%M%S'; }

rsync_one() {
  local label="$1"
  local src_path="$2"
  local dest_path="$3"
  local log="${LOG_DIR}/rosmap_${label}_$(timestamp).log"

  echo "========================================"
  echo "[$(date '+%F %T')] ${label}"
  echo "  SRC : ${SRC_HOST}:${src_path}"
  echo "  DEST: ${dest_path}"
  echo "  LOG : ${log}"
  echo "========================================"

  # Trailing slash on src dir => copy contents into dest
  rsync "${RSYNC_OPTS[@]}" \
    "${SRC_HOST}:${src_path}" \
    "${dest_path}" \
    2>&1 | tee "${log}"

  echo "[$(date '+%F %T')] Done: ${label}"
  echo
}

echo "ROSMAP transfer starting (dry_run=${DRY_RUN})"
echo "DEST_ROOT=${DEST_ROOT}"
if [[ -f "${SSH_KEY}" ]]; then
  echo "SSH_KEY=${SSH_KEY}"
else
  echo "WARNING: ${SSH_KEY} not found — will prompt for password on first rsync"
fi
echo

# 1) Metadata (staged on CAMH from external_data + nethome → netdata)
rsync_one "metadata" \
  "${STAGE_META}/" \
  "${DEST_ROOT}/Metadata/"

# 2) ROSMAP WGS genotypes — raw (already on netdata)
rsync_one "wgs_raw" \
  "${WGS_BASE}/ROSMAP_joint_WGS_vcf/" \
  "${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf/"

# 3) ROSMAP WGS genotypes — normalized (already on netdata)
rsync_one "wgs_normalized" \
  "${WGS_BASE}/ROSMAP_joint_WGS_vcf_normalized/" \
  "${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf_normalized/"

# 4) ROSMAP_array genotypes — raw TOPmed (must be staged onto netdata on CAMH first)
#    Includes n1686/, n381/, merged/ (~433G). Stage with stage_rosmap_on_camh.sh.
#    Note: cannot locally test STAGE_TOPMED from Trillium — it lives on CAMH.
if [[ "${SKIP_TOPMED}" -eq 1 ]]; then
  echo "SKIP_TOPMED=1 — skipping TOPmed raw transfer"
else
  rsync_one "array_topmed_raw" \
    "${STAGE_TOPMED}/" \
    "${DEST_ROOT}/Genotype/ROSMAP_array_TOPmed_imputed_vcf/"
fi

# 5) ROSMAP_array genotypes — normalized (already on netdata)
rsync_one "array_normalized" \
  "${WGS_BASE}/ROSMAP_array_vcf_normalized/" \
  "${DEST_ROOT}/Genotype/ROSMAP_array_vcf_normalized/"

# 6) Sample ID maps (staged)
rsync_one "array_sample_maps" \
  "${STAGE_SAMPLE_MAPS}/" \
  "${DEST_ROOT}/Genotype/ROSMAP_array_sample_maps/"

echo "All ROSMAP transfers finished."
echo "Layout:"
echo "  ${DEST_ROOT}/Metadata/ROSMAP_biospecimen_metadata.csv"
echo "  ${DEST_ROOT}/Metadata/ROSmaster.rds"
echo "  ${DEST_ROOT}/Metadata/ROSMAP_assay_wholeGenomeSeq_metadata.csv"
echo "  ${DEST_ROOT}/Metadata/ROSMAP_assay_snpArray_metadata.csv"
echo "  ${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf/                 # WGS raw"
echo "  ${DEST_ROOT}/Genotype/ROSMAP_joint_WGS_vcf_normalized/      # WGS normalized"
echo "  ${DEST_ROOT}/Genotype/ROSMAP_array_TOPmed_imputed_vcf/      # array raw (TOPmed)"
echo "  ${DEST_ROOT}/Genotype/ROSMAP_array_vcf_normalized/          # array normalized"
echo "  ${DEST_ROOT}/Genotype/ROSMAP_array_sample_maps/             # ID maps for prep"
echo "Existing RNA (not transferred):"
echo "  ${DEST_ROOT}/ROSMAP_Raw_Counts_Bulk"
echo "  ${DEST_ROOT}/ROSMAP_RNAseq_provenance"
