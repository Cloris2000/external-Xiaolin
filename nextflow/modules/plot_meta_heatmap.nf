/*
 * PLOT_META_HEATMAP  (aggregate step)
 *
 * Runs once after all per-CT PLOT_META_RESULTS jobs complete.
 * Reads the intermediate TSVs written by those jobs to produce:
 *   - Cross-cell-type P-value heatmap
 *   - Combined I² overview panel (all cell types)
 *   - meta_analysis_summary.tsv
 *
 * Does NOT re-read any .tbl files — only consumes tiny intermediate TSVs
 * (typically a few MB total), so memory requirement is small.
 */

process PLOT_META_HEATMAP {
    label 'small_memory'
    tag  "plot_heatmap"

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    input:
    tuple path(intermediate_files),     // collected *_top_hits / *_het_stats / *_summary_stats TSVs
          path(plot_script),
          val(output_dir),
          val(p_thresh),
          val(gw_thresh)

    output:
    path "plots/heatmap/**",                      emit: heatmap_files, optional: true
    path "plots/het/all_celltypes_*.png",          emit: het_overview,  optional: true
    path "plots/het/all_celltypes_het_summary.tsv",emit: het_summary,   optional: true
    path "meta_analysis_summary.tsv",              emit: summary_tsv,   optional: true

    script:
    """
    mkdir -p plots/heatmap plots/het intermediate

    # Stage all collected intermediate TSVs into a single flat directory
    for f in ${intermediate_files}; do
        cp "\$f" intermediate/
    done

    Rscript ${plot_script} \
        --mode             aggregate \
        --intermediate_dir intermediate \
        --output_dir       plots \
        --p_thresh         ${p_thresh} \
        --gw_thresh        ${gw_thresh}

    [ -f plots/meta_analysis_summary.tsv ] && \
        cp plots/meta_analysis_summary.tsv meta_analysis_summary.tsv || true
    """
}
