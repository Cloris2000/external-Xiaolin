/*
 * Module: Covariate Preparation (Step 4 - After PCA Calculation)
 * Merges PCA with clinical covariates
 * REQUIRES PCA FROM GENOTYPE QC
 */

process COV_PREP {
    label 'low_memory'
    
    publishDir "${params.outdir}/${params.study}", mode: 'copy', overwrite: true
    
    input:
    path pca_file
    path clinical_cov_file
    path samples_file
    path script_file
    val output_dir
    val study_name
    
    output:
    path "covariates.txt", emit: covariate_file
    
    script:
    def conda_init = '''
        # Set variable to prevent unbound variable error in conda deactivation scripts
        export xml_catalog_files_libxml2="${xml_catalog_files_libxml2:-}"
        # Initialize conda with error suppression
        set +u
        if [ -f "$HOME/.anaconda3/etc/profile.d/conda.sh" ]; then
            source "$HOME/.anaconda3/etc/profile.d/conda.sh" 2>/dev/null || true
        elif [ -f "/nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh" ]; then
            source "/nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh" 2>/dev/null || true
        elif [ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]; then
            source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
        fi
        set -u
    '''
    def conda_activate = params.conda_env ? """
        set +u
        if ! conda activate ${params.conda_env} 2>/dev/null; then
            echo "ERROR: Failed to activate conda environment '${params.conda_env}'" >&2
            exit 1
        fi
        set -u
        echo 'Activated conda env: ${params.conda_env}'
        if ! which R >/dev/null 2>&1; then
            echo "ERROR: R not found after conda activation" >&2
            exit 1
        fi
    """ : ""
    """
    ${conda_init}
    ${conda_activate}
    Rscript "${script_file}" \\
        --pca_file "${pca_file}" \\
        --clinical_cov_file "${clinical_cov_file}" \\
        --samples_file "${samples_file}" \\
        --output_file "covariates.txt" \\
        --output_dir "${output_dir}"
    """
}

