/*
 * Module: Convert logp to p
 * Converts log10 p-values to p-values in Regenie output files for METAL meta-analysis
 */

process CONVERT_LOGP_TO_P {
    label 'low_memory'
    tag "${cohort}_${cell_type}_logp_to_p"
    
    input:
    tuple val(cohort), val(cell_type), path(regenie_file), val(output_dir)
    
    output:
    path "${cohort}_${cell_type}_step2.regenie.raw_p", emit: raw_p_file
    path "${cohort}_${cell_type}_logp_to_p.done", emit: done_file
    
    publishDir "${output_dir}", mode: 'copy', overwrite: true
    
    script:
    """
    awk 'NR==1 {print \$0, "P"} NR>1 {print \$0, 10^(-\$12)}' ${regenie_file} > ${cohort}_${cell_type}_step2.regenie.raw_p
    
    touch ${cohort}_${cell_type}_logp_to_p.done
    """
}

