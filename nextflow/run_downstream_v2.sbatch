#!/bin/bash
#SBATCH --job-name=downstream_v2
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/zhoux156/logs/downstream_v2/downstream_v2_%j.out
#SBATCH --error=/scratch/zhoux156/logs/downstream_v2/downstream_v2_%j.err

# Downstream analyses + manuscript figures for the corrected 15-cohort meta (v2).
#
# Phase 1 (parallel on the node)
#   heterogeneity breakdown   run_heterogeneity_breakdown_array.sbatch  (per-cohort effects from the
#                             harmonized hg19 files, so GTEx_v10/AMP-AD/ROSMAP_array are included)
#   bulk vs sn concordance    run_concordance_hodge5u.sbatch  (bulk = v2 meta; sn = SCC hodge5u meta)
#   colocalisation            run_coloc_all_loci.sbatch       (all v2 suggestive loci x 9 disease GWAS)
# Phase 2
#   coloc summary (04) + regional plots (05); heterogeneity figures + forest plots
# Phase 3  manuscript figures (all written to scratch; /project is read-only here)
#   GWAS cache -> Figure 3 (+ multi-CT Manhattan) -> Figure 5 -> Figure 4 + coloc rearranged
#   -> Figure 6 v2 (data-driven forest plots) -> Figure 2 (proportion phenotypes, validation, merged)
#
# The existing drivers are reused by exporting their variables and running them with
# bash (sbatch here is wrapped with --export=NONE, so nothing can be passed from the
# submitting shell; everything is set inside this job).
#
# Optional argument: path to a v2 LDSC h2 summary for Figure 2's h2 panel; falls back to
# the SCC (old) summary and says so.   sbatch run_downstream_v2.sbatch [ldsc_h2_summary.tsv]
#
# Decision recorded here: coloc's --cell_type_N is 3554 (the individuals analysed in v2,
# sum of REGENIE N across the 15 cohorts) instead of the drivers' stale default 2812.
# Trillium allocates whole nodes; do not add --mem.

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

V2META="${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2"
V2COHORTS="${RESULTS_ROOT}/cohort_gwas_v2/results"
DS="${RESULTS_ROOT}/downstream_v2"
FIG="${DS}/manuscript_figures"
LOGD="${LOG_ROOT}/downstream_v2/${SLURM_JOB_ID:-manual}"
SCC_RES="${SCC_DIR}/results"
H2_FILE="${1:-}"
mkdir -p "${DS}" "${FIG}" "${LOGD}"

echo "downstream_v2  $(date -Is)"
echo "  v2 meta    : ${V2META}"
echo "  v2 cohorts : ${V2COHORTS}"
echo "  outputs    : ${DS}"
echo "  figures    : ${FIG}"
echo "  logs       : ${LOGD}"

# ---- preflight ---------------------------------------------------------------
fail=0
[ "$(ls "${V2META}"/*_meta_analysis_*.tbl 2>/dev/null | grep -vc '1\.tbl$')" -eq 19 ] || { echo "FAIL: expected 19 v2 .tbl"; fail=1; }
[ "$(ls "${V2META}"/*.annotated.tsv 2>/dev/null | wc -l)" -eq 19 ] || { echo "FAIL: expected 19 annotated.tsv"; fail=1; }
[ "$(find "${V2META}/harmonized" -name '*_harmonized.raw_p' 2>/dev/null | wc -l)" -eq 285 ] || { echo "FAIL: expected 285 harmonized files"; fail=1; }
[ -d "${SCC_RES}/coloc/disease_gwas" ] || { echo "FAIL: disease GWAS dir missing"; fail=1; }
[ "$(ls "${SCC_RES}/meta_analysis_sn_hodge5u"/*_meta_analysis_*.tbl 2>/dev/null | wc -l)" -gt 0 ] || { echo "FAIL: sn hodge5u meta missing"; fail=1; }
for c in ROSMAP NIMH_HBCC_1M MSBB CMC_MSSM; do [ -s "${V2COHORTS}/${c}/cell_proportions.csv" ] || { echo "FAIL: ${c} cell_proportions.csv"; fail=1; }; done
[ "${fail}" -eq 0 ] || { echo "Preflight failed; aborting."; exit 1; }
echo "  preflight OK"

step() {  # step <name> <cmd...>  -> logs to ${LOGD}/<name>.log, records OK/FAIL, never aborts the job
    local name="$1"; shift
    local t0; t0=$(date +%s)
    if "$@" > "${LOGD}/${name}.log" 2>&1; then
        printf "  OK    %-28s %4ss\n" "${name}" "$(( $(date +%s) - t0 ))"
    else
        printf "  FAIL  %-28s %4ss  (see %s)\n" "${name}" "$(( $(date +%s) - t0 ))" "${LOGD}/${name}.log"; echo "${name}" >> "${DS}/.failed_steps"
    fi
}
rm -f "${DS}/.failed_steps"

# ---- Phase 1: three analyses in parallel --------------------------------------
echo ""; echo "Phase 1  $(date -Is)"
(
  export META_DIR="${V2META}" OUT_DIR="${DS}/heterogeneity_breakdown"
  export SUMSTAT_PATTERN="${V2META}/harmonized/{cell_type}/harmonized/{cohort}_{cell_type}_harmonized.raw_p"
  step het_breakdown bash "${NF_DIR}/run_heterogeneity_breakdown_array.sbatch"
) &
(
  export BULK_DIR="${V2META}" SN_DIR="${SCC_RES}/meta_analysis_sn_hodge5u"
  export OUTDIR="${DS}/sn_bulk_meta_similarity_hodge5u/top_hits"
  step concordance_hodge5u bash "${NF_DIR}/run_concordance_hodge5u.sbatch"
) &
(
  export META_DIR="${V2META}" COLOC_LOCI_DIR="${DS}/coloc/loci_full" COLOC_RES_DIR="${DS}/coloc/coloc_results_full"
  # Disease GWAS on hg19.  The SCC copies of AD_Bellenguez2022 / LBD_Chia2021 / PD_Nalls2019
  # are GRCh38 despite the _hg19 name (GWAS-Catalog harmonised coordinates); the lifted set
  # is produced by run_downstream_v2_redo2.sbatch (scripts/coloc/liftover_disease_gwas.py).
  export DISEASE_DIR="${DS}/disease_gwas_hg19" CELL_TYPE_N=3554
  [ -s "${DISEASE_DIR}/AD_Bellenguez2022/AD_Bellenguez2022_hg19.tsv" ] || \
      echo "  WARNING: ${DISEASE_DIR} missing - run run_downstream_v2_redo2.sbatch (Phase 0) first" >&2
  # Primary, most-recent GWAS per disease only (the six the manuscript figures map);
  # the SCC dir also holds superseded sets under archive/ which must not be scanned.
  export COLOC_DISEASES="MDD_MDD2025,BD_bip2024,SCZ_Trubetskoy2022,AD_Bellenguez2022,LBD_Chia2021,PD_Nalls2019"
  step coloc_all_loci bash "${NF_DIR}/run_coloc_all_loci.sbatch"
) &
wait

# ---- Phase 2: summaries + regional plots + heterogeneity figures ---------------
echo ""; echo "Phase 2  $(date -Is)"
(
  step coloc_summarize "${RSCRIPT}" "${NF_DIR}/scripts/coloc/04_summarize_coloc.R" \
      --results_dir "${DS}/coloc/coloc_results_full" --loci_dir "${DS}/coloc/loci_full" --output_dir "${DS}/coloc/full" \
  && step coloc_regional_plots "${RSCRIPT}" "${NF_DIR}/scripts/coloc/05_regional_coloc_plots.R" \
      --sig_file "${DS}/coloc/full/coloc_significant.tsv" --loci_dir "${DS}/coloc/loci_full" \
      --disease_dir "${SCC_RES}/coloc/disease_gwas" --meta_dir "${V2META}" \
      --output_dir "${DS}/coloc/full/plots/regional" --genome_build 37
) &
(
  HB="${DS}/heterogeneity_breakdown"; mkdir -p "${HB}/figures" "${HB}/forest_plots"
  step het_plots "${RSCRIPT}" "${NF_DIR}/scripts/plot_heterogeneity_breakdown.R" \
      --per-cohort-effects "${HB}/per_cohort_effects.tsv" --group-i2-summary "${HB}/group_i2_summary.tsv" \
      --overall-group-i2 "${HB}/overall_group_i2.tsv" --output-dir "${HB}/figures" --forest-dir "${HB}/forest_plots" \
      --max-forest-plots 100
) &
wait

# ---- Phase 3: manuscript figures ----------------------------------------------
echo ""; echo "Phase 3  $(date -Is)"
MF="${NF_DIR}/manuscript_figure"
export MF_OUT_DIR="${FIG}" FIG2_OUT_DIR="${FIG}"
export MF_META_DIR="${V2META}" MF_CACHE_DIR="${FIG}/gwas_cache"
export MF_RESULTS_DIR="${DS}"                       # figure5: ${DS}/sn_bulk_meta_similarity_hodge5u/top_hits
export MF_COLOC_DIR="${DS}/coloc/full" MF_LOCI_DIR="${DS}/coloc/loci_full" MF_DISEASE_DIR="${SCC_RES}/coloc/disease_gwas"
export MF_AUTO_REGIONAL=1
export MF_BULK_RESULTS_DIR="${V2COHORTS}"
[ -n "${H2_FILE}" ] && export MF_H2_FILE="${H2_FILE}"   # only consumed by the (skipped) h2 figure
cd "${FIG}"

step fig3_gwas_cache        "${RSCRIPT}" "${MF}/preprocess_gwas_to_cache.R"
step fig3_gwas_discovery    "${RSCRIPT}" "${MF}/figure3_gwas_discovery.R"
step fig3_multict_manhattan "${RSCRIPT}" "${MF}/test_multict_gws_manhattan.R"
step fig5_sn_concordance    "${RSCRIPT}" "${MF}/figure5_sn_concordance.R"
step fig4_coloc_validation  "${RSCRIPT}" "${MF}/figure4_coloc_validation.R"
step fig_coloc_rearranged   "${RSCRIPT}" "${MF}/figure_coloc_rearranged.R"
# Figure 6 loci = the three colocalisation events with the highest PP.H4 (one per
# locus); falls back to the manuscript trio if the coloc summary is missing.
step fig6_heterogeneity_v2  "${RSCRIPT}" "${MF}/figure6_cohort_heterogeneity_v2.R" \
    --meta_dir "${V2META}" --metadata "${NF_DIR}/docs/meta_cohort_metadata.tsv" --out_dir "${FIG}" \
    --from_coloc "${DS}/coloc/full/coloc_all_results.tsv" --loci_dir "${DS}/coloc/loci_full" --n_loci 3
step fig2_proportion_phenos "${RSCRIPT}" "${MF}/figure2_celltype_proportion_phenotypes.R"
step fig2_validation_focus  "${RSCRIPT}" "${MF}/figure2_celltype_validation_focus.R"
# figure2_ctp_variation_merged.R (the "revised" Figure 2 with the LDSC h2 panel D) is
# deliberately NOT run: that figure is not used (2026-09-16).  To produce it later:
#   MF_H2_FILE=<ldsc_h2_summary.tsv> FIG2_OUT_DIR=<dir> Rscript figure2_ctp_variation_merged.R

echo ""; echo "Done  $(date -Is)"
echo "Figures ($(ls "${FIG}"/*.png 2>/dev/null | wc -l) png):"; ls "${FIG}"/*.png 2>/dev/null | sed 's/^/  /'
if [ -s "${DS}/.failed_steps" ]; then echo ""; echo "FAILED steps:"; sed 's/^/  /' "${DS}/.failed_steps"; exit 1; fi
