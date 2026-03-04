/*
 * Module: Cell Type Proportion Estimation
 * Flexible module supporting multiple deconvolution tools and reference taxonomies
 */

process CELL_TYPE_DECONV {
    label 'medium_memory'
    
    input:
    path corrected_data
    path metadata
    path script_file
    val deconv_tool
    val reference_taxonomy
    path marker_file
    path hgnc_mapping_file
    val output_dir
    
    output:
    path "cell_proportions.csv", emit: cell_proportions
    path "cell_proportions_scaled.csv", emit: cell_proportions_scaled
    path "deconv_summary.txt", emit: summary, optional: true
    
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
        --corrected_data "${corrected_data}" \\
        --metadata "${metadata}" \\
        --deconv_tool "${deconv_tool}" \\
        --reference_taxonomy "${reference_taxonomy}" \\
        --marker_file "${marker_file}" \\
        --hgnc_mapping_file "${hgnc_mapping_file}" \\
        --output_dir "${output_dir}" \\
        --proportions_output "cell_proportions.csv" \\
        --proportions_scaled_output "cell_proportions_scaled.csv"
    """
}