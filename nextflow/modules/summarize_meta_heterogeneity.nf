/*
 * SUMMARIZE_META_HETEROGENEITY
 *
 * Aggregates METAL heterogeneity fields across all cell types and writes:
 *   1) per-cell-type heterogeneity summary
 *   2) suggestive/high-heterogeneity variant table for MR-MEGA triage
 */

process SUMMARIZE_META_HETEROGENEITY {
    label 'low_memory'
    tag 'meta_heterogeneity_summary'

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    input:
    tuple path(meta_tbls),
          path(script_file),
          val(output_dir),
          val(gw_p_thresh),
          val(suggestive_p_thresh),
          val(i2_thresh),
          val(het_p_thresh)

    output:
    path "meta_heterogeneity_summary.tsv", emit: summary_tsv
    path "meta_heterogeneity_top_hits.tsv", emit: top_hits_tsv

    script:
    def file_list = meta_tbls instanceof List ? meta_tbls : [meta_tbls]
    def tbl_args = file_list.collect { "\"${it}\"" }.join(' ')
    """
    python3 "${script_file}" \\
        --meta-tbls ${tbl_args} \\
        --output-summary-tsv meta_heterogeneity_summary.tsv \\
        --output-top-hits-tsv meta_heterogeneity_top_hits.tsv \\
        --gw-p-thresh "${gw_p_thresh}" \\
        --suggestive-p-thresh "${suggestive_p_thresh}" \\
        --i2-thresh "${i2_thresh}" \\
        --het-p-thresh "${het_p_thresh}"
    """
}
