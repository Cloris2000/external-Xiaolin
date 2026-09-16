#!/bin/bash
#SBATCH --job-name=meta_15cohorts_hg19_v2
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_hg19_v2_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_hg19_v2_%j.err

# Corrected 15-cohort meta-analysis (v2) over the cohort_gwas_v2 GWAS.  Output is
# hg19; the four hg38 cohorts (GTEx_v10, AMP_AD_Rush, AMP_AD_Mayo, ROSMAP_array)
# are lifted down to hg19 before harmonization.
#
# Requires all 15 cohorts to have finished run_cohort_gwas_v2.sbatch
# (19 raw_p files each under ${RESULTS_ROOT}/cohort_gwas_v2/results/<cohort>/regenie_step2).
#
# Runs LIFTOVER -> HARMONIZE -> METAL -> allele annotation -> heterogeneity
# summary -> plots, then compares each cell type's .tbl with the SCC
# meta_analysis_15cohorts table (Effect / StdErr / P correlation, lead SNP).
# Expected vs SCC: Freq1 ~ 1 - old Freq1 (bug fix), ROSMAP_array now contributes
# ~5.6M variants instead of 32k, and cohort-level differences propagate.
#
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
source "${NF_DIR}/scripts/cohort_gwas_v2_lib.sh"
module load "${NEXTFLOW_MODULE}"

WORK_DIR="${WORK_ROOT}/meta_15cohorts_hg19_v2"
RUN_DIR="${SCRATCH_ROOT}/run_dirs/meta_15cohorts_hg19_v2"
OUT_DIR="${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2"
CMP_DIR="${OUT_DIR}/compare_vs_scc"
SCC_META="${SCC_RES}/meta_analysis_15cohorts"
mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${OUT_DIR}" "${CMP_DIR}"
cd "${RUN_DIR}"

echo "Checking cohort_gwas_v2 REGENIE output for all 15 cohorts..."
all_good=true
for cohort in "${ALL_COHORTS[@]}"; do
    count=$(ls "${V2_RESULTS}/${cohort}"/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
    if [ "$count" -eq 19 ]; then echo "  OK  ${cohort}: ${count}/19"; else echo "  FAIL ${cohort}: ${count}/19 (INCOMPLETE)"; all_good=false; fi
    if [ -s "${V2_RESULTS}/${cohort}/deconv_orientation.tsv" ] && grep -q NEGATIVE "${V2_RESULTS}/${cohort}/deconv_orientation.tsv"; then
        echo "  WARN ${cohort}: deconv_orientation.tsv has a NEGATIVE cell type"; all_good=false
    fi
done
[ "$all_good" = true ] || { echo "ERROR: cohort GWAS incomplete or orientation problems unresolved. Aborting."; exit 1; }

CHAIN="${REFS_ROOT}/reference_data/hg38ToHg19.over.chain.gz"
[ -s "${CHAIN}" ] || { echo "ERROR: liftover chain file missing: ${CHAIN}" >&2; exit 1; }
echo "Liftover chain present; GTEx_v10, AMP_AD_Rush, AMP_AD_Mayo, ROSMAP_array will be lifted hg38 -> hg19."

nextflow run "${NF_DIR}/standalone_meta_analysis.nf" \
    -c "${NF_DIR}/nextflow.config.standalone_meta_15cohorts_hg19_v2" \
    -work-dir "${WORK_DIR}" \
    -resume

echo ""
echo "Comparing v2 meta tables with the SCC meta_analysis_15cohorts tables..."
rows=()
for tbl in "${OUT_DIR}"/*_meta_analysis_*.tbl; do
    ct=$(basename "${tbl}" | sed 's/_meta_analysis_.*//')
    scc=$(ls "${SCC_META}/${ct}_meta_analysis_"*1.tbl 2>/dev/null | head -1)   # SCC kept METAL's <prefix>1.tbl
    [ -s "${scc}" ] || scc=$(ls "${SCC_META}/${ct}_meta_analysis_"*.tbl 2>/dev/null | head -1)
    if [ ! -s "${scc}" ]; then echo "  ${ct}: no SCC table"; continue; fi
    "${RSCRIPT}" "${NF_DIR}/scripts/compare_raw_p.R" --cohort META15 --cell_type "${ct}" \
        --new "${tbl}" --old "${scc}" --out "${CMP_DIR}/${ct}.tsv" >/dev/null 2>&1 && rows+=("${CMP_DIR}/${ct}.tsv") || echo "  ${ct}: comparison failed"
done
if [ ${#rows[@]} -gt 0 ]; then
    { head -1 "${rows[0]}"; for r in "${rows[@]}"; do tail -n +2 "$r"; done; } > "${CMP_DIR}/summary.tsv"
    column -t "${CMP_DIR}/summary.tsv"
fi
echo ""
echo "Allele annotation summaries:"; cat "${OUT_DIR}"/*.annotated.summary.txt 2>/dev/null | sort | uniq -c | sort -k2 | head -20

echo ""
echo "Corrected 15-cohort meta-analysis complete.  Results: ${OUT_DIR}/"
echo "Sync to ${PROJECT_RESULTS}/ from a login node with:"
echo "  ${NF_DIR}/scripts/sync_results_to_project.sh meta_analysis_15cohorts_hg19_v2"
