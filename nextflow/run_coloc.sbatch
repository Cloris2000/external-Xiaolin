#!/bin/bash
#SBATCH --job-name=coloc_15cohorts
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=8:00:00
#SBATCH --output=/scratch/zhoux156/logs/coloc/coloc_%j.out
#SBATCH --error=/scratch/zhoux156/logs/coloc/coloc_%j.err

# Original coloc driver (suggestive-SNP windows, min_snps=50, susie on).
# All 19 cell types run as parallel workers on one whole node.
# Prefer run_coloc_all_loci.sbatch for the manuscript full-window analysis.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

META_DIR="${META_DIR:-${RESULTS_ROOT}/meta_analysis_15cohorts}"
COLOC_LOCI_DIR="${COLOC_LOCI_DIR:-${RESULTS_ROOT}/coloc/loci}"
COLOC_RES_DIR="${COLOC_RES_DIR:-${RESULTS_ROOT}/coloc/coloc_results}"
DISEASE_DIR="${DISEASE_DIR:-${SCC_DIR}/results/coloc/disease_gwas}"
SCRIPTS_DIR="${NF_DIR}/scripts/coloc"
CELL_TYPE_N="${CELL_TYPE_N:-2812}"
WORKER_LOGS="${LOG_ROOT}/coloc/workers_${SLURM_JOB_ID:-manual}"

mkdir -p "${COLOC_LOCI_DIR}" "${COLOC_RES_DIR}" "${WORKER_LOGS}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

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
    if [ ! -f "${loci}" ] || [ ! -f "${ldata}" ]; then
        "${RSCRIPT}" "${SCRIPTS_DIR}/01_extract_loci.R" \
            --cell_type   "${ct}" \
            --tbl_file    "${tbl}" \
            --output_dir  "${COLOC_LOCI_DIR}" \
            --p_threshold 1e-5 \
            --window_kb   500 \
            > "${WORKER_LOGS}/${ct}.log" 2>&1
    fi

    local out="${COLOC_RES_DIR}/${ct}_coloc_results.tsv"
    if [ -f "${out}" ]; then
        echo "  SKIP ${ct} (exists: ${out})"
        return 0
    fi

    "${RSCRIPT}" "${SCRIPTS_DIR}/03_run_coloc.R" \
        --cell_type    "${ct}" \
        --loci_file    "${loci}" \
        --locus_data   "${ldata}" \
        --disease_dir  "${DISEASE_DIR}" \
        --output_dir   "${COLOC_RES_DIR}" \
        --cell_type_N  "${CELL_TYPE_N}" \
        --run_susie    TRUE \
        --pp_h4_susie  0.5 \
        --min_snps     50 \
        >> "${WORKER_LOGS}/${ct}.log" 2>&1 \
        && echo "  OK   ${ct}" \
        || { echo "  FAIL ${ct} (see ${WORKER_LOGS}/${ct}.log)"; return 1; }
}
export -f run_one
export META_DIR COLOC_LOCI_DIR COLOC_RES_DIR DISEASE_DIR SCRIPTS_DIR
export CELL_TYPE_N RSCRIPT WORKER_LOGS

printf '%s\n' "${CELL_TYPES[@]}" \
  | xargs -P "${#CELL_TYPES[@]}" -I{} bash -c 'run_one "$@"' _ {} \
  && rc=0 || rc=$?

echo ""
if [ "$rc" -ne 0 ]; then
    echo "ERROR: at least one cell type failed. Logs: ${WORKER_LOGS}"
    exit 1
fi
echo "Coloc complete: ${COLOC_RES_DIR}"
