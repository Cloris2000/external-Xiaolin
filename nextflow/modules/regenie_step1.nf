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
    def snp_filter = params.regenie_step1_snp_filter
    if (!snp_filter && params.june9_replay) { snp_filter = 'mac1' }
    if (!snp_filter) { snp_filter = 'hom_alt' }
    """
    EXTRACT_FILE="${prune_in_file}"
    if [ -n "${plink_path}" ] && [ -x "${plink_path}" ]; then
        awk 'NR>1 && NF>=3 && \$3!="" && \$3!="NA" {print \$1,\$2}' ${pheno_file} > keep_regenie_samples.txt
        if [ -s keep_regenie_samples.txt ]; then
            if [ "${snp_filter}" = "mac1" ]; then
                # April 7 / June 9: drop only MAC=0 in the phenotyped subset.
                ${plink_path} --pfile ${pgen_prefix} --keep keep_regenie_samples.txt \\
                    --extract ${prune_in_file} --mac 1 --out prune_in_mac1 --write-snplist 2>/dev/null || true
                if [ -s prune_in_mac1.snplist ]; then
                    EXTRACT_FILE="prune_in_mac1.snplist"
                fi
            else
                # Current default: also require a hom-ALT (added 2026-07-06).
                ${plink_path} --pfile ${pgen_prefix} --keep keep_regenie_samples.txt \\
                    --extract ${prune_in_file} --geno-counts --out snp_geno_subset 2>/dev/null || true
                if [ -s snp_geno_subset.gcount ]; then
                    awk 'NR>1 && \$7+0 >= 1 && (\$6+0 + 2*\$7+0) >= 2 {print \$2}' snp_geno_subset.gcount > prune_in_filtered.snplist
                    if [ -s prune_in_filtered.snplist ]; then
                        EXTRACT_FILE="prune_in_filtered.snplist"
                    fi
                fi
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

