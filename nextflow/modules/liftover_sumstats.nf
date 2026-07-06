/*
 * Module: Liftover Summary Statistics (hg38 → hg19)
 *
 * Converts CHROM, GENPOS, and ID columns in a single cohort's
 * .regenie.raw_p file from hg38 to hg19 using a UCSC chain file
 * (via pyliftover).  Variants that cannot be mapped are excluded and
 * written to an unmapped log.
 *
 * This process runs once per (cell_type, cohort) pair for cohorts
 * whose genome_build is "hg38" in params.cohort_genome_builds.
 */

process LIFTOVER_SUMSTATS {
    label 'small_memory'
    tag "${cohort}_${cell_type}_liftover"

    input:
    tuple val(cell_type),
          val(cohort),
          path(raw_p_file),
          path(chain_file),
          path(liftover_script),
          val(output_dir)

    output:
    tuple val(cell_type),
          val(cohort),
          path("${cohort}_${cell_type}.hg19_lifted.regenie.raw_p"), emit: lifted
    path "${cohort}_${cell_type}_hg38_unmapped.tsv",               emit: unmapped_log

    publishDir "${output_dir}/liftover_logs", mode: 'copy', overwrite: true,
               saveAs: { name -> name.endsWith('.tsv') ? name : null }

    script:
    """
    python3 ${liftover_script} \
        --input-file  ${raw_p_file} \
        --cohort      ${cohort} \
        --chain-file  ${chain_file} \
        --output-file ${cohort}_${cell_type}.hg19_lifted.regenie.raw_p \
        --unmapped-log ${cohort}_${cell_type}_hg38_unmapped.tsv
    """
}
