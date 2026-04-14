/*
 * PLOT_META_RESULTS
 *
 * Generates three sets of visualizations from METAL meta-analysis output:
 *   1. Manhattan + QQ plots     (per cell type)
 *   2. I² heterogeneity histograms (per cell type + combined overview)
 *   3. Cross-cell-type P-value heatmap
 *
 * Input:
 *   results_dir  — directory containing *.tbl files (no "1.tbl" stale files)
 *   plot_script  — path to plot_meta_results.R
 *   output_dir   — where to write all plots
 *
 * Output:
 *   plot_dir     — directory path of generated plots (for downstream steps)
 */

process PLOT_META_RESULTS {
    label 'medium_memory'
    tag  "plot_meta"

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(results_dir),
          path(plot_script),
          val(output_dir),
          val(p_thresh),
          val(gw_thresh)

    output:
    path "plots/**",        emit: plot_files
    path "meta_analysis_summary.tsv", emit: summary_tsv, optional: true

    script:
    """
    mkdir -p plots

    Rscript ${plot_script} \
        --results_dir ${results_dir} \
        --output_dir  plots \
        --p_thresh    ${p_thresh} \
        --gw_thresh   ${gw_thresh}

    # Copy summary table to top level for publishDir pickup
    [ -f plots/meta_analysis_summary.tsv ] && \
        cp plots/meta_analysis_summary.tsv meta_analysis_summary.tsv || true
    """
}
