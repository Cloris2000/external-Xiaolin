/*
 * Module: Harmonize Meta Summary Statistics
 *
 * Runs pre-METAL harmonization for a single cell type across all cohorts:
 *   - normalise variant IDs to chr:pos:ref:alt
 *   - remove strand-ambiguous A/T and C/G SNPs
 *   - apply minimum MAF filter
 *   - require variant present in >= min_present_cohorts
 *   - flag (and optionally hard-filter) cross-cohort AF discrepancies
 *   - write a per-cell-type drop log and summary TSV
 */

process HARMONIZE_META_SUMSTATS {
    label 'medium_memory'
    tag "${cell_type}_harmonize"

    input:
    tuple val(cell_type),
          path(raw_p_files),
          val(cohort_names_json),
          path(harmonize_script),
          val(output_dir),
          val(drop_strand_ambiguous),
          val(min_meta_maf),
          val(min_present_cohorts),
          val(af_delta_report_threshold),
          val(af_delta_filter_threshold)

    output:
    tuple val(cell_type), path("harmonized/*.raw_p"),          emit: harmonized_files
    path "${cell_type}_drop_log.tsv",                          emit: drop_log
    path "${cell_type}_harmonization_summary.tsv",             emit: summary

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    script:
    def file_list = raw_p_files instanceof List ? raw_p_files : [raw_p_files]
    def files_arg = file_list.collect { it.toString() }.join(" ")
    def af_filter_arg = af_delta_filter_threshold ? "--af-delta-filter-threshold ${af_delta_filter_threshold}" : ""
    """
    mkdir -p harmonized

    python3 ${harmonize_script} \
        --input-files ${files_arg} \
        --cohort-names '${cohort_names_json}' \
        --cell-type ${cell_type} \
        --output-dir harmonized \
        --drop-log ${cell_type}_drop_log.tsv \
        --summary-tsv ${cell_type}_harmonization_summary.tsv \
        --drop-strand-ambiguous ${drop_strand_ambiguous} \
        --min-meta-maf ${min_meta_maf} \
        --min-present-cohorts ${min_present_cohorts} \
        --af-delta-report-threshold ${af_delta_report_threshold} \
        ${af_filter_arg}
    """
}
