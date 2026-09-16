/*
 * Module: Annotate METAL output with explicit REF/ALT and effect-allele columns
 *
 * Writes <cell_type>_meta_analysis_<suffix>.annotated.tsv next to the .tbl,
 * plus a small count summary.  The .tbl itself is left untouched so scripts
 * that read it by column position (coloc, LDSC prep) keep working.
 * See scripts/annotate_meta_alleles.py for column definitions.
 */

process ANNOTATE_META_ALLELES {
    label 'small_memory'
    tag "${cell_type}_annotate"

    input:
    tuple val(cell_type), path(meta_tbl), path(script), val(output_dir)

    output:
    tuple val(cell_type), path("${meta_tbl.baseName}.annotated.tsv"), emit: annotated
    path "${meta_tbl.baseName}.annotated.summary.txt",                 emit: summary

    publishDir "${output_dir}", mode: 'copy', overwrite: true

    script:
    """
    python3 ${script} \\
        --tbl ${meta_tbl} \\
        --out ${meta_tbl.baseName}.annotated.tsv \\
        --summary ${meta_tbl.baseName}.annotated.summary.txt
    """
}
