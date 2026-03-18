/*
 * Module: Phenotype Preparation (Step 1 - Before Genotype QC)
 * Creates phenotypes, sample list, and clinical covariates
 * NO PCA OR GENOTYPE DATA REQUIRED
 */

process PHENO_PREP {
    label 'low_memory'
    
    publishDir "${params.outdir}/${params.study}", mode: 'copy', overwrite: true
    
    input:
    path cell_proportions
    path metadata
    path script_file
    path biospecimen_file  // Optional - will be empty if not provided
    val output_dir
    val study_name
    val fid_method
    val fid_format
    val col_msex
    val col_age
    val col_individualID
    val clinical_metadata_file
    val biospec_col_individual
    val biospec_col_specimen
    
    output:
    path "phenotypes_RINT.txt", emit: phenotype_file
    path "samples_with_phenotypes.txt", emit: samples_file
    path "clinical_covariates.txt", emit: clinical_cov_file
    
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
        --cell_proportions "${cell_proportions}" \\
        --metadata "${metadata}" \\
        ${biospecimen_file.name != '.empty_biospecimen' ? "--biospecimen_file \"${biospecimen_file}\"" : ""} \\
        --output_dir "${output_dir}" \\
        --study "${study_name}" \\
        --phenotype_output "phenotypes_RINT.txt" \\
        --samples_output "samples_with_phenotypes.txt" \\
        --clinical_cov_output "clinical_covariates.txt" \\
        --fid_method "${fid_method}" \\
        ${fid_format && fid_format != '' ? "--fid_format \"${fid_format}\"" : ""} \\
        ${col_msex && col_msex != '' ? "--col_msex \"${col_msex}\"" : ""} \\
        ${col_age && col_age != '' ? "--col_age \"${col_age}\"" : ""} \\
        ${col_individualID && col_individualID != '' ? "--col_individualID \"${col_individualID}\"" : ""} \\
        ${clinical_metadata_file && clinical_metadata_file != '' ? "--clinical_metadata \"${clinical_metadata_file}\"" : ""} \\
        ${biospec_col_individual && biospec_col_individual != '' ? "--biospec_col_individual \"${biospec_col_individual}\"" : ""} \\
        ${biospec_col_specimen && biospec_col_specimen != '' ? "--biospec_col_specimen \"${biospec_col_specimen}\"" : ""}
    """
}

