#!/bin/bash
#SBATCH --job-name=coloc_all_loci
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/zhoux156/logs/coloc/coloc_full_%j.out
#SBATCH --error=/scratch/zhoux156/logs/coloc/coloc_full_%j.err

# Colocalize ALL suggestive loci (p < 1e-5) using full 500 kb windows.
# All 19 cell types run as parallel workers on one whole node.
#
#   META_DIR=... RESULTS_DIR=... sbatch run_coloc_all_loci.sbatch
#
# After the job finishes, summarize from a login node:
#   Rscript scripts/coloc/04_summarize_coloc.R \
#     --results_dir $RESULTS_ROOT/coloc/coloc_results_full \
#     --loci_dir    $RESULTS_ROOT/coloc/loci_full \
#     --output_dir  $RESULTS_ROOT/coloc/full

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

META_DIR="${META_DIR:-${RESULTS_ROOT}/meta_analysis_15cohorts}"
COLOC_LOCI_DIR="${COLOC_LOCI_DIR:-${RESULTS_ROOT}/coloc/loci_full}"
COLOC_RES_DIR="${COLOC_RES_DIR:-${RESULTS_ROOT}/coloc/coloc_results_full}"
DISEASE_DIR="${DISEASE_DIR:-${SCC_DIR}/results/coloc/disease_gwas}"
SCRIPTS_DIR="${NF_DIR}/scripts/coloc"
CELL_TYPE_N="${CELL_TYPE_N:-2812}"
MIN_SNPS="${MIN_SNPS:-10}"
# Optional comma-separated disease labels (basename without _hg19.tsv) passed to
# 03_run_coloc.R --disease_include.  Unset = every *_hg19.tsv under DISEASE_DIR,
# INCLUDING the superseded sets kept in DISEASE_DIR/archive/ (MDD_Wray2018,
# BD_Mullins2021), which the manuscript run did not use.
COLOC_DISEASES="${COLOC_DISEASES:-}"
WORKER_LOGS="${LOG_ROOT}/coloc/workers_full_${SLURM_JOB_ID:-manual}"

mkdir -p "${COLOC_LOCI_DIR}" "${COLOC_RES_DIR}" "${WORKER_LOGS}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

echo "Meta tables : ${META_DIR}"
echo "Loci        : ${COLOC_LOCI_DIR}"
echo "Results     : ${COLOC_RES_DIR}"
echo "Disease GWAS: ${DISEASE_DIR}"

run_one() {
    local ct="$1"
    local tbl
    tbl=$(ls "${META_DIR}/${ct}"_meta_analysis_*.tbl 2>/dev/null | grep -v '[0-9]\.tbl$' | head -1)
    if [ -z "${tbl}" ]; then
        echo "  FAIL ${ct}: no .tbl in ${META_DIR}" >&2
        return 1
    fi

    local loci="${COLOC_LOCI_DIR}/${ct}_loci.tsv"
    local ldata="${COLOC_LOCI_DIR}/${ct}_locus_data.tsv.gz"
    if [ -s "${loci}" ] && [ -s "${ldata}" ]; then
        echo "  skip extract ${ct}" >> "${WORKER_LOGS}/${ct}.log"
    else
        "${RSCRIPT}" "${SCRIPTS_DIR}/01_extract_loci.R" \
            --cell_type        "${ct}" \
            --tbl_file         "${tbl}" \
            --output_dir       "${COLOC_LOCI_DIR}" \
            --p_threshold      1e-5 \
            --window_kb        500 \
            --full_window_snps TRUE \
            > "${WORKER_LOGS}/${ct}.log" 2>&1
    fi

    "${RSCRIPT}" "${SCRIPTS_DIR}/03_run_coloc.R" \
        --cell_type    "${ct}" \
        --loci_file    "${COLOC_LOCI_DIR}/${ct}_loci.tsv" \
        --locus_data   "${COLOC_LOCI_DIR}/${ct}_locus_data.tsv.gz" \
        --disease_dir  "${DISEASE_DIR}" \
        --output_dir   "${COLOC_RES_DIR}" \
        --cell_type_N  "${CELL_TYPE_N}" \
        --run_susie    FALSE \
        --min_snps     "${MIN_SNPS}" \
        ${COLOC_DISEASES:+--disease_include "${COLOC_DISEASES}"} \
        >> "${WORKER_LOGS}/${ct}.log" 2>&1 \
        && echo "  OK   ${ct}" \
        || { echo "  FAIL ${ct} (see ${WORKER_LOGS}/${ct}.log)"; return 1; }
}
export -f run_one
export META_DIR COLOC_LOCI_DIR COLOC_RES_DIR DISEASE_DIR SCRIPTS_DIR
export CELL_TYPE_N MIN_SNPS RSCRIPT WORKER_LOGS COLOC_DISEASES

printf '%s\n' "${CELL_TYPES[@]}" \
  | xargs -P "${#CELL_TYPES[@]}" -I{} bash -c 'run_one "$@"' _ {} \
  && rc=0 || rc=$?

echo ""
if [ "$rc" -ne 0 ]; then
    echo "ERROR: at least one cell type failed. Logs: ${WORKER_LOGS}"
    exit 1
fi
echo "Coloc complete: ${COLOC_RES_DIR}"
