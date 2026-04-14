/*
 * Module: PCA and Technical Covariate Computation
 * Processes count matrix CSV file, computes PCA, and generates z-score normalized data
 */

process PCA_TECH_COV {
    label 'high_memory'
    
    input:
    path count_matrix_file
    path metadata_file
    path script_file
    val tissue_filter
    val output_dir
    val col_sample_id_for_matching
    
    output:
    path "zscore_data.RData", emit: zscore_data
    path "metadata_DLPFC.RData", emit: metadata
    path "combined_metrics.csv", emit: combined_metrics
    
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
        --count_matrix_file "${count_matrix_file}" \\
        --metadata_file "${metadata_file}" \\
        --tissue_filter "${tissue_filter}" \\
        --output_dir "${output_dir}" \\
        --zscore_output "zscore_data.RData" \\
        --metadata_output "metadata_DLPFC.RData" \\
        --metrics_output "combined_metrics.csv" \\
        ${col_sample_id_for_matching && col_sample_id_for_matching != '' ? "--col_sample_id_for_matching \"${col_sample_id_for_matching}\"" : ""} \\
        ${params.extra_metrics_file ? "--extra_metrics_file \"${params.extra_metrics_file}\"" : ""}
    """
}