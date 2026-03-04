/*
 * Module: Regenie Step 1
 * Performs Regenie step 1 (whole-genome regression) for a specific cohort and cell type
 */

process REGENIE_STEP1 {
    label 'high_memory'
    tag "${cohort}_${cell_type}_step1"
    
    input:
    tuple val(cohort), val(cell_type), path(pgen_files), path(prune_in_file), path(pheno_file), path(covar_file), val(regenie_path), val(threads), val(bsize), val(output_dir)
    
    output:
    path "${cohort}_${cell_type}_step1_pred.list", emit: pred_list
    path "${cohort}_${cell_type}_step1.log", emit: log_file, optional: true
    path "${cohort}_${cell_type}_step1.done", emit: done_file
    
    publishDir "${output_dir}", mode: 'copy', overwrite: true
    
    script:
    // pgen_files is a list of [.pgen, .pvar, .psam] - get the prefix from the first file
    def pgen_file = pgen_files[0]
    def pgen_prefix = pgen_file.toString().replace('.pgen', '')
    def plink_path = params.get('plink_path', '')
    """
    # Pre-filter prune.in to SNPs with MAC>=1 in analysis samples (avoids REGENIE "low variance" error)
    # REGENIE fails on monomorphic SNPs in the phenotyped subset; PLINK filter prevents that.
    EXTRACT_FILE="${prune_in_file}"
    if [ -n "${plink_path}" ] && [ -x "${plink_path}" ]; then
        awk 'NR>1 && NF>=3 && \$3!="" && \$3!="NA" {print \$1,\$2}' ${pheno_file} > keep_regenie_samples.txt
        if [ -s keep_regenie_samples.txt ]; then
            # plink2 --write-snplist writes to <--out prefix>.snplist (no filename as arg)
            ${plink_path} --pfile ${pgen_prefix} --keep keep_regenie_samples.txt --extract ${prune_in_file} --mac 1 --out prune_in_mac1 --write-snplist 2>/dev/null || true
            if [ -s prune_in_mac1.snplist ]; then
                EXTRACT_FILE="prune_in_mac1.snplist"
            fi
        fi
    fi

    # Run REGENIE step1
    ${regenie_path} \\
        --step 1 \\
        --threads ${threads} \\
        --verbose \\
        --pgen ${pgen_prefix} \\
        --extract \${EXTRACT_FILE} \\
        --phenoFile ${pheno_file} \\
        --covarFile ${covar_file} \\
        --bsize ${bsize} \\
        --out ${cohort}_${cell_type}_step1
    
    touch ${cohort}_${cell_type}_step1.done
    """
}

