#!/bin/bash
#SBATCH --job-name=coloc_bd_mdd_upd
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=8:00:00
#SBATCH --output=/scratch/zhoux156/logs/coloc/coloc_upd_%j.out
#SBATCH --error=/scratch/zhoux156/logs/coloc/coloc_upd_%j.err

# Re-run coloc.abf for updated BD (bip2024) and MDD (MDD2025) only,
# then merge into coloc_results_full replacing Mullins/Wray rows.
# All 19 cell types run as parallel workers on one whole node.

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

COLOC_LOCI_DIR="${COLOC_LOCI_DIR:-${RESULTS_ROOT}/coloc/loci_full}"
COLOC_EXISTING="${COLOC_EXISTING:-${RESULTS_ROOT}/coloc/coloc_results_full}"
COLOC_NEW="${COLOC_NEW:-${RESULTS_ROOT}/coloc/coloc_results_updated}"
DISEASE_DIR="${DISEASE_DIR:-${SCC_DIR}/results/coloc/disease_gwas}"
SCRIPTS_DIR="${NF_DIR}/scripts/coloc"
CELL_TYPE_N="${CELL_TYPE_N:-2812}"
MIN_SNPS="${MIN_SNPS:-10}"
DISEASE_INCLUDE="BD_bip2024,MDD_MDD2025"
REPLACE_DISEASES="BD_Mullins2021,MDD_Wray2018"
WORKER_LOGS="${LOG_ROOT}/coloc/workers_upd_${SLURM_JOB_ID:-manual}"

mkdir -p "${COLOC_NEW}" "${WORKER_LOGS}"

CELL_TYPES=(
  Astrocyte Endothelial IT L4.IT L5.6.IT.Car3 L5.6.NP
  L5.ET L6.CT L6b LAMP5 Microglia OPC Oligodendrocyte
  PAX6 PVALB Pericyte SST VIP VLMC
)

for f in "${DISEASE_DIR}/BD_bip2024/BD_bip2024_hg19.tsv" \
         "${DISEASE_DIR}/MDD_MDD2025/MDD_MDD2025_hg19.tsv"; do
    if [ ! -s "${f}" ]; then
        echo "ERROR: Missing standardized disease file: ${f}" >&2
        exit 1
    fi
done

run_one() {
    local ct="$1"
    local loci="${COLOC_LOCI_DIR}/${ct}_loci.tsv"
    local ldata="${COLOC_LOCI_DIR}/${ct}_locus_data.tsv.gz"
    if [ ! -f "${loci}" ] || [ ! -f "${ldata}" ]; then
        echo "  FAIL ${ct}: missing loci in ${COLOC_LOCI_DIR}" >&2
        return 1
    fi

    "${RSCRIPT}" "${SCRIPTS_DIR}/03_run_coloc.R" \
        --cell_type       "${ct}" \
        --loci_file       "${loci}" \
        --locus_data      "${ldata}" \
        --disease_dir     "${DISEASE_DIR}" \
        --output_dir      "${COLOC_NEW}" \
        --cell_type_N     "${CELL_TYPE_N}" \
        --run_susie       FALSE \
        --min_snps        "${MIN_SNPS}" \
        --disease_include "${DISEASE_INCLUDE}" \
        > "${WORKER_LOGS}/${ct}.log" 2>&1

    "${RSCRIPT}" "${SCRIPTS_DIR}/merge_coloc_results.R" \
        --existing_dir     "${COLOC_EXISTING}" \
        --new_file         "${COLOC_NEW}/${ct}_coloc_results.tsv" \
        --replace_diseases "${REPLACE_DISEASES}" \
        --output_dir       "${COLOC_EXISTING}" \
        >> "${WORKER_LOGS}/${ct}.log" 2>&1 \
        && echo "  OK   ${ct}" \
        || { echo "  FAIL ${ct} (see ${WORKER_LOGS}/${ct}.log)"; return 1; }
}
export -f run_one
export COLOC_LOCI_DIR COLOC_EXISTING COLOC_NEW DISEASE_DIR SCRIPTS_DIR
export CELL_TYPE_N MIN_SNPS DISEASE_INCLUDE REPLACE_DISEASES RSCRIPT WORKER_LOGS

printf '%s\n' "${CELL_TYPES[@]}" \
  | xargs -P "${#CELL_TYPES[@]}" -I{} bash -c 'run_one "$@"' _ {} \
  && rc=0 || rc=$?

echo ""
if [ "$rc" -ne 0 ]; then
    echo "ERROR: at least one cell type failed. Logs: ${WORKER_LOGS}"
    exit 1
fi
echo "Updated BD/MDD coloc merged into ${COLOC_EXISTING}"
