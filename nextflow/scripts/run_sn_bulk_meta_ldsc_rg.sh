#!/usr/bin/env bash
# Run LDSC genetic correlation for snRNAseq-derived vs bulk-derived meta-analysis pairs.

set -euo pipefail

NF_DIR="${NF_DIR:-/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow}"
OUTDIR="${OUTDIR:-${NF_DIR}/results/sn_bulk_meta_similarity}"
MANIFEST="${MANIFEST:-${OUTDIR}/ldsc/ldsc_input_manifest.tsv}"
SUMSTATS_DIR="${OUTDIR}/ldsc/sumstats"
RG_DIR="${OUTDIR}/ldsc/rg_results"
LOG_DIR="${OUTDIR}/ldsc/logs"

PYTHON27="${PYTHON27:-/nethome/kcni/xzhou/.anaconda3/envs/ldsc_py2/bin/python}"
MUNGE_PY="${MUNGE_PY:-/external/rprshnas01/kcni/mwainberg/ldsc/munge_sumstats.py}"
LDSC_PY="${LDSC_PY:-/external/rprshnas01/kcni/mwainberg/ldsc/ldsc.py}"
REF_LD_CHR="${REF_LD_CHR:-/external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr/}"
W_LD_CHR="${W_LD_CHR:-/external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr/}"

mkdir -p "${SUMSTATS_DIR}" "${RG_DIR}" "${LOG_DIR}"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: LDSC input manifest not found: ${MANIFEST}" >&2
  echo "Run scripts/compare_sn_bulk_meta_similarity.R first." >&2
  exit 1
fi

munge_one() {
  local input_tsv="$1"
  local label="$2"
  local output_prefix="${SUMSTATS_DIR}/${label}"
  local output_gz="${output_prefix}.sumstats.gz"

  if [[ -s "${output_gz}" ]]; then
    echo "[skip] Munged sumstats exist: ${output_gz}"
    return
  fi
  if [[ ! -s "${input_tsv}" ]]; then
    echo "ERROR: LDSC input missing or empty for ${label}: ${input_tsv}" >&2
    return 1
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
    --out "${output_prefix}" \
    > "${LOG_DIR}/${label}_munge.log" 2>&1
}

run_rg() {
  local sn_ct="$1"
  local bulk_ct="$2"
  local sn_gz="${SUMSTATS_DIR}/sn_${sn_ct}.sumstats.gz"
  local bulk_gz="${SUMSTATS_DIR}/bulk_${bulk_ct}.sumstats.gz"
  local out_prefix="${RG_DIR}/${sn_ct}_vs_${bulk_ct}"

  if [[ -s "${out_prefix}.log" ]]; then
    echo "[skip] LDSC rg exists: ${out_prefix}.log"
    return
  fi
  if [[ ! -s "${sn_gz}" || ! -s "${bulk_gz}" ]]; then
    echo "ERROR: Missing munged sumstats for ${sn_ct} vs ${bulk_ct}" >&2
    return 1
  fi

  echo "[rg] ${sn_ct} vs ${bulk_ct}"
  "${PYTHON27}" "${LDSC_PY}" \
    --rg "${sn_gz},${bulk_gz}" \
    --ref-ld-chr "${REF_LD_CHR}" \
    --w-ld-chr "${W_LD_CHR}" \
    --out "${out_prefix}" \
    > "${LOG_DIR}/${sn_ct}_vs_${bulk_ct}_rg.stdout.log" 2>&1
}

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sn_ct bulk_ct bulk_input bulk_found sn_input sn_rows; do
  if [[ "${bulk_found}" != "TRUE" && "${bulk_found}" != "True" && "${bulk_found}" != "true" ]]; then
    echo "ERROR: Bulk LDSC input missing for ${sn_ct} vs ${bulk_ct}: ${bulk_input}" >&2
    continue
  fi
  echo "=== ${sn_ct} vs ${bulk_ct} ==="
  munge_one "${sn_input}" "sn_${sn_ct}"
  munge_one "${bulk_input}" "bulk_${bulk_ct}"
  run_rg "${sn_ct}" "${bulk_ct}"
done

python3 "${NF_DIR}/scripts/summarize_ldsc_rg_results.py" \
  --manifest "${MANIFEST}" \
  --rg-dir "${RG_DIR}" \
  --output-tsv "${RG_DIR}/ldsc_rg_summary.tsv"

echo "LDSC rg summary written to ${RG_DIR}/ldsc_rg_summary.tsv"
