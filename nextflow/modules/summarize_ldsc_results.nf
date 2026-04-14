/*
 * SUMMARIZE_LDSC_RESULTS
 *
 * Parses LDSC log files into one summary TSV.
 */

process SUMMARIZE_LDSC_RESULTS {
    label 'low_memory'
    tag 'ldsc_summary'

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    input:
    tuple path(ldsc_logs),
          path(script_file),
          val(output_dir)

    output:
    path "ldsc_h2_summary.tsv", emit: summary_tsv

    script:
    def log_list = ldsc_logs instanceof List ? ldsc_logs : [ldsc_logs]
    def log_args = log_list.collect { "\"${it}\"" }.join(' ')
    """
    python3 "${script_file}" \\
        --ldsc-logs ${log_args} \\
        --output-tsv ldsc_h2_summary.tsv
    """
}
