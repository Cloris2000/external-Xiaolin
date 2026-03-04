/*
 * Module: Phenotype and Covariate File Preparation for GWAS
 * Prepares final phenotype and covariate files with RINT transformation
 */

process PHENO_COV_PREP {
    label 'low_memory'
    
    input:
    path cell_proportions
    path metadata
    path script_file
    path wgs_psam_file
    path pca_file
    path biospecimen_file  // Optional file - will be empty channel if not provided
    val output_dir
    val study_name  // Study/cohort name (e.g., "ROSMAP", "CMC_MSSM")
    val fid_method
    val fid_format
    val col_msex  // Optional column name for sex/gender
    val biospec_col_individual  // Optional column name for individual ID in biospecimen file
    val biospec_col_specimen  // Optional column name for specimen ID in biospecimen file
    
    output:
    path "phenotypes_RINT.txt", emit: phenotype_file
    path "covariates.txt", emit: covariate_file
    path "samples_with_phenotypes.txt", emit: samples_file
    
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
        --cell_proportions "${cell_proportions}" \\
        --metadata "${metadata}" \\
        --wgs_psam_file "${wgs_psam_file}" \\
        --pca_file "${pca_file}" \\
        ${biospecimen_file ? "--biospecimen_file \"${biospecimen_file}\"" : ""} \\
        --output_dir "${output_dir}" \\
        --study "${study_name}" \\
        --phenotype_output "phenotypes_RINT.txt" \\
        --covariate_output "covariates.txt" \\
        --samples_output "samples_with_phenotypes.txt" \\
        --fid_method "${fid_method}" \\
        ${fid_format && fid_format != '' ? "--fid_format \"${fid_format}\"" : ""} \\
        ${col_msex && col_msex != '' ? "--col_msex \"${col_msex}\"" : ""} \\
        ${biospec_col_individual && biospec_col_individual != '' ? "--biospec_col_individual \"${biospec_col_individual}\"" : ""} \\
        ${biospec_col_specimen && biospec_col_specimen != '' ? "--biospec_col_specimen \"${biospec_col_specimen}\"" : ""}
    """
}