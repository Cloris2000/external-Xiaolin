#!/bin/bash
#SBATCH --job-name=cohort_manhattan_v2
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/zhoux156/logs/downstream_v2/cohort_manhattan_v2_%j.out
#SBATCH --error=/scratch/zhoux156/logs/downstream_v2/cohort_manhattan_v2_%j.err

# Per-cell-type Manhattan plots for ONE cohort's v2 REGENIE step-2 results
# (scripts/generate_cohort_manhattan.R: topr, hg19, suggestive SNPs gene-annotated).
#   sbatch run_cohort_manhattan_v2.sbatch ROSMAP          # default ROSMAP (bulk WGS)
# Reads  ${RESULTS_ROOT}/cohort_gwas_v2/results/<cohort>/regenie_step2/<cohort>_<ct>_step2.regenie.raw_p
# Writes ${RESULTS_ROOT}/cohort_gwas_v2/plots/<cohort>/manhattan/<cohort>_<ct>_manhattan_annotated.png
#        + <cohort>_manhattan_summary.tsv (variants, lambda, min P, hits)
# All 19 cell types run in parallel (each reads a ~770 MB table). Whole-node allocation; no --mem.

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
COHORT="${1:-ROSMAP}"
IN="${RESULTS_ROOT}/cohort_gwas_v2/results/${COHORT}/regenie_step2"
OUT="${RESULTS_ROOT}/cohort_gwas_v2/plots/${COHORT}/manhattan"
LOGD="${LOG_ROOT}/downstream_v2/${SLURM_JOB_ID:-manual}_manhattan_${COHORT}"
mkdir -p "${OUT}" "${LOGD}"
echo "cohort_manhattan_v2  ${COHORT}  $(date -Is)"
[ -d "${IN}" ] || { echo "missing ${IN}"; exit 1; }

pids=()
for f in "${IN}"/${COHORT}_*_step2.regenie.raw_p; do
    ct=$(basename "${f}" _step2.regenie.raw_p); ct=${ct#${COHORT}_}
    ( "${RSCRIPT}" "${NF_DIR}/scripts/generate_cohort_manhattan.R" --input "${f}" --cell_type "${ct}" --cohort "${COHORT}" \
          --output_dir "${OUT}" --output_prefix "${COHORT}_${ct}" > "${LOGD}/${ct}.log" 2>&1 \
      && echo "  OK    ${ct}" || echo "  FAIL  ${ct}  (see ${LOGD}/${ct}.log)" ) &
    pids+=($!)
done
wait "${pids[@]}"

# Summary table from the per-cell-type logs
{
  printf "cohort\tcell_type\tn_variants\tlambda\tmin_p\tn_p5e8\tn_p1e5\tplot\n"
  for l in "${LOGD}"/*.log; do
      ct=$(basename "${l}" .log)
      nv=$(grep -m1 "Total variants:"  "${l}" | awk '{print $NF}')
      la=$(grep -m1 "Lambda:"          "${l}" | awk '{print $NF}')
      mp=$(grep -m1 "Minimum p-value:" "${l}" | awk '{print $NF}')
      n8=$(grep -m1 "p < 5e-8" "${l}" | awk '{print $NF}')
      n5=$(grep -m1 "p < 1e-5" "${l}" | awk '{print $NF}')
      png="${OUT}/${COHORT}_${ct}_manhattan_annotated.png"; [ -s "${png}" ] || png="MISSING"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${COHORT}" "${ct}" "${nv:-NA}" "${la:-NA}" "${mp:-NA}" "${n8:-NA}" "${n5:-NA}" "${png}"
  done
} > "${OUT}/${COHORT}_manhattan_summary.tsv"
echo ""; column -t "${OUT}/${COHORT}_manhattan_summary.tsv" | cut -c1-120
echo "Done  $(date -Is)   plots: ${OUT}"
