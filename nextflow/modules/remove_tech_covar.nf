/*
 * Module: Remove Batch Effects and Technical Covariates
 * Identifies top tech covariates and removes batch effects
 */

process REMOVE_TECH_COVAR {
    label 'high_memory'
    
    input:
    path zscore_data
    path metadata
    path script_file
    val top_n_tech_cov
    val output_dir
    
    output:
    path "corrected_data.RData", emit: corrected_data
    path "metadata_cleaned.csv", emit: metadata_cleaned
    path "top_tech_covariates.txt", emit: top_tech_cov
    path "eigencor_plot*.png", emit: plots, optional: true
    
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
        # Activate conda - suppress stderr warnings but check if activation succeeded
        if ! conda activate ${params.conda_env} 2>/dev/null; then
            echo "ERROR: Failed to activate conda environment '${params.conda_env}'" >&2
            exit 1
        fi
        set -u
        echo 'Activated conda env: ${params.conda_env}'
        # Verify R is available
        if ! which R >/dev/null 2>&1; then
            echo "ERROR: R not found after conda activation" >&2
            exit 1
        fi
    """ : ""
    """
    ${conda_init}
    ${conda_activate}
    Rscript "${script_file}" \\
        --zscore_data "${zscore_data}" \\
        --metadata "${metadata}" \\
        --top_n_tech_cov ${top_n_tech_cov} \\
        --tech_cov_mode "${params.tech_cov_mode}" \\
        ${params.tech_covariates_file ? "--tech_covariates_file \"${params.tech_covariates_file}\"" : ""} \\
        --batch_covariates "${params.batch_covariates}" \\
        ${params.disable_batch_correction ? "--disable_batch_correction" : ""} \\
        ${params.col_diagnosis ? "--col_diagnosis \"${params.col_diagnosis}\"" : ""} \\
        ${params.col_msex ? "--col_msex \"${params.col_msex}\"" : ""} \\
        ${params.batch_recode ? "--batch_recode \"${params.batch_recode}\"" : ""} \\
        --output_dir "${output_dir}" \\
        --corrected_output "corrected_data.RData" \\
        --metadata_output "metadata_cleaned.csv" \\
        --tech_cov_output "top_tech_covariates.txt"
    """
}