/*
 * HETEROGENEITY_BREAKDOWN
 *
 * Decomposes meta-analysis heterogeneity across four biological and technical
 * dimensions (cohorts, diagnosis groups, ancestry classes, analysis types):
 *
 *   1. heterogeneity_breakdown.py  — extracts per-cohort betas from REGENIE
 *      sumstats for all suggestive top hits, computes within-group Cochran's Q
 *      and I² for each dimension, writes three TSVs.
 *
 *   2. plot_heterogeneity_breakdown.R — generates forest plots (one PDF per top
 *      hit), dimension summary violin plots, grouped I² heatmaps, and a combined
 *      multi-panel figure.
 *
 * Inputs:
 *   meta_tbls           — list of METAL .tbl files (one per cell type)
 *   cohort_sumstat_dirs — list of per-cohort regenie_step2 directories
 *   cohort_metadata     — docs/meta_cohort_metadata.tsv (with diagnosis_group
 *                         and ancestry_class columns)
 *   py_script           — path to scripts/heterogeneity_breakdown.py
 *   r_script            — path to scripts/plot_heterogeneity_breakdown.R
 *   output_dir          — publishDir target
 *   cohorts             — space-separated ordered cohort list (matches METAL
 *                         Direction string order, i.e. alphabetical)
 *   results_root        — root results directory (default "results")
 *   gw_p_thresh         — genome-wide p threshold (default 5e-8)
 *   suggestive_p_thresh — suggestive p threshold (default 1e-5)
 *   max_forest_plots    — max number of forest plot PDFs to generate (default 50)
 */

process HETEROGENEITY_BREAKDOWN {
    label 'medium_memory'
    tag 'heterogeneity_breakdown'

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    input:
    tuple path(meta_tbls),
          path(cohort_metadata),
          path(py_script),
          path(r_script),
          val(output_dir),
          val(cohorts),
          val(results_root),
          val(sumstat_pattern),
          val(gw_p_thresh),
          val(suggestive_p_thresh),
          val(max_forest_plots)

    output:
    path "heterogeneity_breakdown/per_cohort_effects.tsv",   emit: per_cohort_effects
    path "heterogeneity_breakdown/group_i2_summary.tsv",     emit: group_i2_summary
    path "heterogeneity_breakdown/overall_group_i2.tsv",     emit: overall_group_i2
    path "heterogeneity_breakdown/forest_plots/*.pdf",        emit: forest_plots,        optional: true
    path "heterogeneity_breakdown/figures/*.pdf",             emit: figures,             optional: true

    script:
    def tbl_list      = meta_tbls instanceof List ? meta_tbls : [meta_tbls]
    def tbl_args      = tbl_list.collect { "\"${it}\"" }.join(' ')
    def cohort_args   = cohorts
    def pat           = sumstat_pattern ?: "${results_root}/{cohort}/regenie_step2/{cohort}_{cell_type}_step2.regenie.raw_p"
    def gw            = gw_p_thresh ?: "5e-8"
    def sug           = suggestive_p_thresh ?: "1e-5"
    def max_f         = max_forest_plots ?: 50
    """
    # Step 1: extract per-cohort effects and compute group-level I²
    python3 "${py_script}" \\
        --meta-tbls ${tbl_args} \\
        --cohort-sumstat-pattern "${pat}" \\
        --cohort-metadata "${cohort_metadata}" \\
        --output-dir heterogeneity_breakdown \\
        --cohorts ${cohort_args} \\
        --sumstat-results-dir "${results_root}" \\
        --gw-p-thresh "${gw}" \\
        --suggestive-p-thresh "${sug}"

    # Step 2: generate forest plots and summary figures
    Rscript "${r_script}" \\
        --per-cohort-effects heterogeneity_breakdown/per_cohort_effects.tsv \\
        --group-i2-summary   heterogeneity_breakdown/group_i2_summary.tsv \\
        --overall-group-i2   heterogeneity_breakdown/overall_group_i2.tsv \\
        --output-dir         heterogeneity_breakdown/figures \\
        --forest-dir         heterogeneity_breakdown/forest_plots \\
        --max-forest-plots   ${max_f}
    """
}
