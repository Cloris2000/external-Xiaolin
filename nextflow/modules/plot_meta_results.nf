/*
 * PLOT_META_RESULTS  (per cell type)
 *
 * Runs for a single cell type: generates Manhattan + QQ + I² histogram plots
 * and writes three small intermediate TSV files used by PLOT_META_HEATMAP:
 *
 *   {cell_type}_top_hits.tsv      — variants passing p_thresh (for heatmap)
 *   {cell_type}_het_stats.tsv     — per-CT I² summary stats
 *   {cell_type}_summary_stats.tsv — n_variants / lambda / n_gw / n_sugg
 *
 * Runs in parallel: one SLURM job per cell type.
 * The aggregate cross-CT plots (heatmap, I² overview, summary table) are
 * produced by PLOT_META_HEATMAP after all per-CT jobs complete.
 */

process PLOT_META_RESULTS {
    label 'small_memory'
    tag  "${cell_type}_plot"

    publishDir "${output_dir}/manhattan", mode: 'copy', overwrite: true,
               saveAs: { name -> name.endsWith(".png") && name.contains("manhattan") ? name : null }
    publishDir "${output_dir}/qq",        mode: 'copy', overwrite: true,
               saveAs: { name -> name.endsWith(".png") && name.contains("qq") ? name : null }
    publishDir "${output_dir}/het",       mode: 'copy', overwrite: true,
               saveAs: { name -> name.endsWith(".png") && name.contains("het") ? name : null }

    input:
    tuple val(cell_type),
          path(tbl_file),
          path(plot_script),
          val(results_dir),
          val(output_dir),
          val(p_thresh),
          val(gw_thresh)

    output:
    tuple val(cell_type), path("plots/**/*.png"),                        emit: plot_files
    tuple val(cell_type), path("plots/${cell_type}_top_hits.tsv"),       emit: top_hits
    tuple val(cell_type), path("plots/${cell_type}_het_stats.tsv"),      emit: het_stats,     optional: true
    tuple val(cell_type), path("plots/${cell_type}_summary_stats.tsv"),  emit: summary_stats

    script:
    """
    mkdir -p plots/manhattan plots/qq plots/het plots/heatmap

    Rscript ${plot_script} \
        --mode        per_celltype \
        --results_dir ${results_dir} \
        --cell_type   ${cell_type} \
        --output_dir  plots \
        --p_thresh    ${p_thresh} \
        --gw_thresh   ${gw_thresh}
    """
}
