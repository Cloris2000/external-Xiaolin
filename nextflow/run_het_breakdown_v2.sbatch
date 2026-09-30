#!/bin/bash
#SBATCH --job-name=het_breakdown_v2
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/zhoux156/logs/downstream_v2/het_breakdown_v2_%j.out
#SBATCH --error=/scratch/zhoux156/logs/downstream_v2/het_breakdown_v2_%j.err

# Heterogeneity breakdown (+ figures / forest plots) for the corrected 15-cohort meta (v2),
# standalone.  Same settings as the het_breakdown step of run_downstream_v2.sbatch:
# per-cohort effects come from the harmonized, hg19-lifted files METAL consumed, so
# GTEx_v10 / AMP-AD / ROSMAP_array are included.  Variables are set here rather than
# passed at submit time because sbatch is wrapped with --export=NONE.
# Trillium allocates whole nodes; do not add --mem.

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

V2META="${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2"
HB="${RESULTS_ROOT}/downstream_v2/heterogeneity_breakdown"
mkdir -p "${HB}/figures" "${HB}/forest_plots"

export META_DIR="${V2META}" OUT_DIR="${HB}"
export SUMSTAT_PATTERN="${V2META}/harmonized/{cell_type}/harmonized/{cohort}_{cell_type}_harmonized.raw_p"
echo "pattern: ${SUMSTAT_PATTERN}"

if bash "${NF_DIR}/run_heterogeneity_breakdown_array.sbatch"; then
    echo "breakdown OK; plotting"
    "${RSCRIPT}" "${NF_DIR}/scripts/plot_heterogeneity_breakdown.R" \
        --per-cohort-effects "${HB}/per_cohort_effects.tsv" --group-i2-summary "${HB}/group_i2_summary.tsv" \
        --overall-group-i2 "${HB}/overall_group_i2.tsv" --output-dir "${HB}/figures" --forest-dir "${HB}/forest_plots" \
        --max-forest-plots 100 && echo "plots OK: ${HB}/figures  ${HB}/forest_plots"
else
    echo "breakdown FAILED"; exit 1
fi
