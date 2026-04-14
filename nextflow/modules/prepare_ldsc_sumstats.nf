/*
 * PREPARE_LDSC_SUMSTATS
 *
 * Converts one METAL table plus the corresponding harmonized cohort files into
 * an LDSC-ready tabular input. Per-variant N is computed by summing the REGENIE
 * N column across the harmonized cohort files that contributed that variant.
 */

process PREPARE_LDSC_SUMSTATS {
    label 'low_memory'
    tag "${cell_type}_ldsc_prep"

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(cell_type),
          path(meta_tbl),
          path(harmonized_files),
          path(script_file),
          val(output_dir),
          val(variant_map)

    output:
    tuple val(cell_type), path("${cell_type}.ldsc_input.tsv"), emit: ldsc_input
    path "${cell_type}.ldsc_input.summary.tsv", emit: summary

    script:
    def file_list = harmonized_files instanceof List ? harmonized_files : [harmonized_files]
    def harmonized_args = file_list.collect { "\"${it}\"" }.join(' ')
    def variant_map_arg = variant_map ? "--variant-map \"${variant_map}\"" : ""
    """
    python3 "${script_file}" \\
        --meta-tsv "${meta_tbl}" \\
        --harmonized-files ${harmonized_args} \\
        --output-tsv "${cell_type}.ldsc_input.tsv" \\
        --summary-tsv "${cell_type}.ldsc_input.summary.tsv" \\
        ${variant_map_arg}
    """
}
