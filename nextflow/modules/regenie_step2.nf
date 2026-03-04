/*
 * Module: Regenie Step 2
 * Performs Regenie step 2 (association testing) for a specific cohort and cell type
 */

process REGENIE_STEP2 {
    label 'high_memory'
    tag "${cohort}_${cell_type}_step2"
    
    input:
    tuple val(cohort), val(cell_type), path(pgen_files), path(pheno_file), path(covar_file), path(pred_list_file), val(regenie_path), val(threads), val(bsize), val(output_dir)
    
    output:
    path "${cohort}_${cell_type}_step2_${cell_type}.regenie", emit: regenie_file
    path "${cohort}_${cell_type}_step2.log", emit: log_file, optional: true
    path "${cohort}_${cell_type}_step2.done", emit: done_file
    
    publishDir "${output_dir}", mode: 'copy', overwrite: true
    
    script:
    // pgen_files is a list of [.pgen, .pvar, .psam] - get the prefix from the first file
    def pgen_file = pgen_files[0]
    def pgen_prefix = pgen_file.toString().replace('.pgen', '')
    """
    ${regenie_path} \\
        --step 2 \\
        --threads ${threads} \\
        --verbose \\
        --pgen ${pgen_prefix} \\
        --phenoFile ${pheno_file} \\
        --covarFile ${covar_file} \\
        --bsize ${bsize} \\
        --pred ${pred_list_file} \\
        --out ${cohort}_${cell_type}_step2
    
    touch ${cohort}_${cell_type}_step2.done
    """
}

