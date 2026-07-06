#!/usr/bin/env bash
#SBATCH --job-name=sn_bulk_ldsc_rg
#SBATCH --partition=mediumtmp
#SBATCH --time=8:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-15%15
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity/ldsc/slurm_logs/sn_bulk_ldsc_rg_%A_%a.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity/ldsc/slurm_logs/sn_bulk_ldsc_rg_%A_%a.err

set -euo pipefail

NF_DIR="${NF_DIR:-/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow}"
OUTDIR="${OUTDIR:-${NF_DIR}/results/sn_bulk_meta_similarity}"
MANIFEST="${MANIFEST:-${OUTDIR}/ldsc/ldsc_input_manifest.tsv}"
SUMSTATS_DIR="${OUTDIR}/ldsc/sumstats"
SNP_INPUT_DIR="${OUTDIR}/ldsc/snp_only_inputs"
RG_DIR="${OUTDIR}/ldsc/rg_results"
LOG_DIR="${OUTDIR}/ldsc/logs"
SLURM_LOG_DIR="${OUTDIR}/ldsc/slurm_logs"

PYTHON27="${PYTHON27:-/nethome/kcni/xzhou/.anaconda3/envs/ldsc_py2/bin/python}"
MUNGE_PY="${MUNGE_PY:-/external/rprshnas01/kcni/mwainberg/ldsc/munge_sumstats.py}"
LDSC_PY="${LDSC_PY:-/external/rprshnas01/kcni/mwainberg/ldsc/ldsc.py}"
REF_LD_CHR="${REF_LD_CHR:-/external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr/}"
W_LD_CHR="${W_LD_CHR:-/external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr/}"
MERGE_ALLELES="${MERGE_ALLELES:-/external/rprshnas01/kcni/mwainberg/ldsc/w_hm3.snplist}"

mkdir -p "${SUMSTATS_DIR}" "${SNP_INPUT_DIR}" "${RG_DIR}" "${LOG_DIR}" "${SLURM_LOG_DIR}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
ROW="$(awk -v task_id="${TASK_ID}" 'BEGIN { FS = "\t" } NR == task_id + 1 { print }' "${MANIFEST}")"
if [[ -z "${ROW}" ]]; then
  echo "ERROR: No manifest row for task ${TASK_ID}" >&2
  exit 1
fi

IFS=$'\t' read -r SN_CT BULK_CT BULK_INPUT BULK_FOUND SN_INPUT SN_ROWS <<< "${ROW}"

echo "======================================================"
echo "Task ${TASK_ID}: ${SN_CT} vs ${BULK_CT}"
echo "Started: $(date)"
echo "SN input: ${SN_INPUT}"
echo "Bulk input: ${BULK_INPUT}"
echo "SN LDSC rows: ${SN_ROWS}"
echo "======================================================"

if [[ "${BULK_FOUND}" != "TRUE" && "${BULK_FOUND}" != "True" && "${BULK_FOUND}" != "true" ]]; then
  echo "ERROR: Bulk LDSC input missing according to manifest: ${BULK_INPUT}" >&2
  exit 1
fi

filter_snp_input() {
  local input_tsv="$1"
  local label="$2"
  local output_tsv="${SNP_INPUT_DIR}/${label}.snp_only.ldsc_input.tsv"

  if [[ -s "${output_tsv}" ]]; then
    echo "[skip] SNP-only input exists: ${output_tsv}"
    printf '%s\n' "${output_tsv}"
    return
  fi

  echo "[filter] SNP-only LDSC input for ${label}" >&2
  awk -F '\t' 'BEGIN { OFS = "\t" }
    NR == 1 { print; next }
    length($2) == 1 && length($3) == 1 && $2 ~ /^[ACGT]$/ && $3 ~ /^[ACGT]$/ { print }
  ' "${input_tsv}" > "${output_tsv}"

  printf '%s\n' "${output_tsv}"
}

munge_one() {
  local input_tsv="$1"
  local label="$2"
  local output_prefix="${SUMSTATS_DIR}/${label}"
  local output_gz="${output_prefix}.sumstats.gz"
  local lock_dir="${output_prefix}.munge.lock"

  if [[ -s "${output_gz}" ]]; then
    if gzip -t "${output_gz}" 2>/dev/null; then
      echo "[skip] Munged sumstats exist and passed gzip check: ${output_gz}"
      return
    fi
    echo "[warn] Existing munged sumstats failed gzip check; rebuilding: ${output_gz}"
    rm -f "${output_gz}" "${output_prefix}.log"
  fi
  if [[ ! -s "${input_tsv}" ]]; then
    echo "ERROR: LDSC input missing or empty for ${label}: ${input_tsv}" >&2
    exit 1
  fi

  while ! mkdir "${lock_dir}" 2>/dev/null; do
    echo "[wait] Another process is munging ${label}; waiting 60s"
    sleep 60
    if [[ -s "${output_gz}" ]]; then
      if gzip -t "${output_gz}" 2>/dev/null; then
        echo "[skip] Munged sumstats appeared and passed gzip check: ${output_gz}"
        return
      fi
      echo "[warn] Munged sumstats appeared but failed gzip check; waiting lock owner"
    fi
  done
  trap 'rm -rf "${lock_dir}"' RETURN

  if [[ -s "${output_gz}" ]]; then
    if gzip -t "${output_gz}" 2>/dev/null; then
      echo "[skip] Munged sumstats exist after lock and passed gzip check: ${output_gz}"
      rm -rf "${lock_dir}"
      trap - RETURN
      return
    fi
    echo "[warn] Existing munged sumstats failed gzip check after lock; rebuilding: ${output_gz}"
    rm -f "${output_gz}" "${output_prefix}.log"
  fi

  echo "[munge] ${label}"
  "${PYTHON27}" "${MUNGE_PY}" \
    --sumstats "${input_tsv}" \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --N-col N \
    --signed-sumstats BETA,0 \
    --chunksize 500000 \
    --merge-alleles "${MERGE_ALLELES}" \
    --out "${output_prefix}" \
    > "${LOG_DIR}/${label}_munge.log" 2>&1

  rm -rf "${lock_dir}"
  trap - RETURN
}

SN_LABEL="sn_${SN_CT}_hm3"
BULK_LABEL="bulk_${BULK_CT}_hm3"
OUT_PREFIX="${RG_DIR}/${SN_CT}_vs_${BULK_CT}"

SN_SNP_INPUT="$(filter_snp_input "${SN_INPUT}" "${SN_LABEL}")"
BULK_SNP_INPUT="$(filter_snp_input "${BULK_INPUT}" "${BULK_LABEL}")"

munge_one "${SN_SNP_INPUT}" "${SN_LABEL}"
munge_one "${BULK_SNP_INPUT}" "${BULK_LABEL}"

if [[ -s "${OUT_PREFIX}.log" ]]; then
  echo "[skip] LDSC rg exists: ${OUT_PREFIX}.log"
else
  echo "[rg] ${SN_CT} vs ${BULK_CT}"
  "${PYTHON27}" "${LDSC_PY}" \
    --rg "${SUMSTATS_DIR}/${SN_LABEL}.sumstats.gz,${SUMSTATS_DIR}/${BULK_LABEL}.sumstats.gz" \
    --ref-ld-chr "${REF_LD_CHR}" \
    --w-ld-chr "${W_LD_CHR}" \
    --out "${OUT_PREFIX}" \
    > "${LOG_DIR}/${SN_CT}_vs_${BULK_CT}_rg.stdout.log" 2>&1
fi

echo "Finished: $(date)"
