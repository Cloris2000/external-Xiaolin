#!/bin/bash
#SBATCH --job-name=downstream_v2_redo
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/zhoux156/logs/downstream_v2/downstream_v2_redo_%j.out
#SBATCH --error=/scratch/zhoux156/logs/downstream_v2/downstream_v2_redo_%j.err

# Redo of the parts of run_downstream_v2.sbatch (job 2325212) affected by three bugs
# found in its output, on the same v2 inputs:
#
#   coloc          the disease scan included the superseded sets in disease_gwas/archive/
#                  (MDD_Wray2018, BD_Mullins2021); now restricted to the six primary GWAS
#                  the manuscript figures map.  Loci (p<1e-5, +/-500 kb) are reused.
#   figure 3       preprocess_gwas_to_cache.R dropped the 'snp' column figure3 needs
#                  -> cache rebuilt, figure 3 + multi-CT Manhattan re-run
#   heterogeneity  per-cohort reader used fixed column indices, so the four cohorts
#                  with an INFO column (GVEX, NIMH_HBCC x3) were silently dropped
#                  -> breakdown + plots re-run with header-based parsing
#   figures 4 / coloc-rearranged / 6 depend on the coloc results -> re-run
#
# Previous v2 coloc outputs are kept alongside as *_incl_archive, not overwritten.
# Trillium allocates whole nodes; do not add --mem.

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

V2META="${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2"
V2COHORTS="${RESULTS_ROOT}/cohort_gwas_v2/results"
DS="${RESULTS_ROOT}/downstream_v2"
FIG="${DS}/manuscript_figures"
LOGD="${LOG_ROOT}/downstream_v2/${SLURM_JOB_ID:-manual}_redo"
SCC_RES="${SCC_DIR}/results"
DISEASES="MDD_MDD2025,BD_bip2024,SCZ_Trubetskoy2022,AD_Bellenguez2022,LBD_Chia2021,PD_Nalls2019"
mkdir -p "${LOGD}" "${FIG}"
echo "downstream_v2_redo  $(date -Is)   diseases: ${DISEASES}"

step() { local name="$1"; shift; local t0; t0=$(date +%s)
    if "$@" > "${LOGD}/${name}.log" 2>&1; then printf "  OK    %-26s %4ss\n" "${name}" "$(( $(date +%s) - t0 ))"
    else printf "  FAIL  %-26s %4ss  (see %s)\n" "${name}" "$(( $(date +%s) - t0 ))" "${LOGD}/${name}.log"; echo "${name}" >> "${DS}/.failed_steps_redo"; fi; }
rm -f "${DS}/.failed_steps_redo"

# ---- keep the archive-contaminated coloc outputs for reference, then redo -------
for d in coloc_results_full full; do
    if [ -d "${DS}/coloc/${d}" ] && [ ! -d "${DS}/coloc/${d}_incl_archive" ]; then
        mv "${DS}/coloc/${d}" "${DS}/coloc/${d}_incl_archive" && echo "  kept previous coloc/${d} as coloc/${d}_incl_archive"
    fi
done

echo ""; echo "Phase 1  $(date -Is)"
(
  export META_DIR="${V2META}" COLOC_LOCI_DIR="${DS}/coloc/loci_full" COLOC_RES_DIR="${DS}/coloc/coloc_results_full"
  export DISEASE_DIR="${SCC_RES}/coloc/disease_gwas" CELL_TYPE_N=3554 COLOC_DISEASES="${DISEASES}"
  step coloc_all_loci bash "${NF_DIR}/run_coloc_all_loci.sbatch" \
  && step coloc_summarize "${RSCRIPT}" "${NF_DIR}/scripts/coloc/04_summarize_coloc.R" \
        --results_dir "${DS}/coloc/coloc_results_full" --loci_dir "${DS}/coloc/loci_full" --output_dir "${DS}/coloc/full" \
  && step coloc_regional_plots "${RSCRIPT}" "${NF_DIR}/scripts/coloc/05_regional_coloc_plots.R" \
        --sig_file "${DS}/coloc/full/coloc_significant.tsv" --loci_dir "${DS}/coloc/loci_full" \
        --disease_dir "${SCC_RES}/coloc/disease_gwas" --meta_dir "${V2META}" \
        --output_dir "${DS}/coloc/full/plots/regional" --genome_build 37
) &
(
  export META_DIR="${V2META}" OUT_DIR="${DS}/heterogeneity_breakdown"
  export SUMSTAT_PATTERN="${V2META}/harmonized/{cell_type}/harmonized/{cohort}_{cell_type}_harmonized.raw_p"
  HB="${DS}/heterogeneity_breakdown"; rm -rf "${HB}/per_celltype" "${HB}/figures" "${HB}/forest_plots"; mkdir -p "${HB}/figures" "${HB}/forest_plots"
  step het_breakdown bash "${NF_DIR}/run_heterogeneity_breakdown_array.sbatch" \
  && step het_plots "${RSCRIPT}" "${NF_DIR}/scripts/plot_heterogeneity_breakdown.R" \
        --per-cohort-effects "${HB}/per_cohort_effects.tsv" --group-i2-summary "${HB}/group_i2_summary.tsv" \
        --overall-group-i2 "${HB}/overall_group_i2.tsv" --output-dir "${HB}/figures" --forest-dir "${HB}/forest_plots" \
        --max-forest-plots 100
) &
(
  # the converter skips existing cache files, so the snp-less cache must go first
  export MF_META_DIR="${V2META}" MF_CACHE_DIR="${FIG}/gwas_cache" MF_OUT_DIR="${FIG}"
  rm -rf "${FIG}/gwas_cache" "${FIG}/figure3_data_checkpoint.rds"
  cd "${FIG}"
  step fig3_gwas_cache        "${RSCRIPT}" "${NF_DIR}/manuscript_figure/preprocess_gwas_to_cache.R" \
  && step fig3_gwas_discovery "${RSCRIPT}" "${NF_DIR}/manuscript_figure/figure3_gwas_discovery.R" \
  && step fig3_multict_manhattan "${RSCRIPT}" "${NF_DIR}/manuscript_figure/test_multict_gws_manhattan.R"
) &
wait

echo ""; echo "Phase 2 (coloc-dependent figures)  $(date -Is)"
MF="${NF_DIR}/manuscript_figure"
export MF_OUT_DIR="${FIG}" MF_META_DIR="${V2META}" MF_AUTO_REGIONAL=1
export MF_COLOC_DIR="${DS}/coloc/full" MF_LOCI_DIR="${DS}/coloc/loci_full" MF_DISEASE_DIR="${SCC_RES}/coloc/disease_gwas"
cd "${FIG}"
step fig4_coloc_validation  "${RSCRIPT}" "${MF}/figure4_coloc_validation.R"
step fig_coloc_rearranged   "${RSCRIPT}" "${MF}/figure_coloc_rearranged.R"
step fig6_heterogeneity_v2  "${RSCRIPT}" "${MF}/figure6_cohort_heterogeneity_v2.R" \
    --meta_dir "${V2META}" --metadata "${NF_DIR}/docs/meta_cohort_metadata.tsv" --out_dir "${FIG}" \
    --from_coloc "${DS}/coloc/full/coloc_all_results.tsv" --loci_dir "${DS}/coloc/loci_full" --n_loci 3

echo ""; echo "Done  $(date -Is)"
echo "coloc events (PP.H4>=0.5): $(( $(wc -l < "${DS}/coloc/full/coloc_significant.tsv" 2>/dev/null || echo 1) - 1 ))"
echo "het cohorts covered: $(awk -F'\t' 'NR>1{c[$6]=1} END{print length(c)}' "${DS}/heterogeneity_breakdown/per_cohort_effects.tsv" 2>/dev/null)"
echo "figures:"; ls -la --time-style=+%H:%M "${FIG}"/*.png 2>/dev/null | awk '{print "  "$6"  "$7}' | sed 's#.*/##'
if [ -s "${DS}/.failed_steps_redo" ]; then echo "FAILED steps:"; sed 's/^/  /' "${DS}/.failed_steps_redo"; exit 1; fi
