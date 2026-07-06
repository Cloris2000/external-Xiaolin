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
    # Pre-filter prune.in to remove zero-variance SNPs in the analysis sample subset.
    # REGENIE step1 crashes with "low variance" on SNPs that have no homozygous-ALT
    # genotypes in the phenotyped subset (e.g. HET-only variants), even when MAC>=2.
    # Root cause: REGENIE's dosage variance = Var(0/1/2 counts); a SNP with only 0s
    # and 1s (no hom-ALT) can have near-zero variance after covariate residualization.
    # Fix: use plink2 --geno-counts to identify and exclude SNPs where TWO_ALT_GENO_CTS==0
    # within the exact sample subset that REGENIE will analyse.
    EXTRACT_FILE="${prune_in_file}"
    if [ -n "${plink_path}" ] && [ -x "${plink_path}" ]; then
        # Build keep list: samples with a non-NA phenotype value (any column ≥3 non-empty)
        awk 'NR>1 && NF>=3 && \$3!="" && \$3!="NA" {print \$1,\$2}' ${pheno_file} > keep_regenie_samples.txt
        if [ -s keep_regenie_samples.txt ]; then
            # Compute per-SNP genotype counts restricted to phenotyped samples.
            # --geno-counts output: #CHROM ID REF ALT HOM_REF_CT HET_REF_ALT_CTS TWO_ALT_GENO_CTS ...
            ${plink_path} --pfile ${pgen_prefix} --keep keep_regenie_samples.txt \
                --extract ${prune_in_file} --geno-counts --out snp_geno_subset 2>/dev/null || true
            # Keep only SNPs where TWO_ALT_GENO_CTS >= 1 (hom-ALT exists) AND HET+2*HOM_ALT >= 2 (MAC>=2).
            # This ensures REGENIE sees genuine dosage variance (not just 0/1 heterozygotes only).
            if [ -s snp_geno_subset.gcount ]; then
                awk 'NR>1 && \$7+0 >= 1 && (\$6+0 + 2*\$7+0) >= 2 {print \$2}' snp_geno_subset.gcount > prune_in_filtered.snplist
                if [ -s prune_in_filtered.snplist ]; then
                    EXTRACT_FILE="prune_in_filtered.snplist"
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

