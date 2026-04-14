/*
 * Module: Meta Cohort Audit
 * Builds a cohort-level QC matrix before multi-cohort meta-analysis.
 */

process META_COHORT_AUDIT {
    label 'low_memory'
    tag "meta_cohort_audit"

    input:
    tuple val(cohorts_json), path(metadata_file), path(audit_script), val(base_dir), val(output_dir)

    output:
    path "meta_cohort_qc_matrix.tsv", emit: audit_tsv
    path "meta_cohort_qc_summary.md", emit: audit_md

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    script:
    """
    python3 ${audit_script} \
        --metadata ${metadata_file} \
        --base-dir ${base_dir} \
        --cohorts-json '${cohorts_json}' \
        --output-tsv meta_cohort_qc_matrix.tsv \
        --output-md meta_cohort_qc_summary.md
    """
}
