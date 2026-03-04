/*
 * Generic Phenotype Processing Pipeline for GWAS Analysis
 * 
 * This pipeline processes gene expression data through:
 * 1. PCA and technical covariate identification
 * 2. Batch effect removal and technical covariate regression
 * 3. Cell type proportion estimation (MGP or other tools)
 * 4. Phenotype and covariate file preparation for GWAS
 * 
 * Supports multiple cohorts: ROSMAP, Mayo, MSBB, CommonMind, NABEC, GTEx
 */

// Include modules
include { PCA_TECH_COV } from './modules/pca_tech_cov.nf'
include { REMOVE_TECH_COVAR } from './modules/remove_tech_covar.nf'
include { CELL_TYPE_DECONV } from './modules/cell_type_deconv.nf'
include { PHENO_COV_PREP } from './modules/pheno_cov_prep.nf'

workflow PHENOTYPE_PIPELINE {
    
    main:
    
    // Stage 1: PCA and technical covariate computation
    PCA_TECH_COV (
        params.count_matrix_file,
        params.metadata_file,
        file("${projectDir}/scripts/pca_tech_cov.R"),
        params.tissue_filter,
        params.output_dir,
        params.col_sample_id_for_matching
    )
    
    // Stage 2: Remove batch effects and top tech covariates
    REMOVE_TECH_COVAR (
        PCA_TECH_COV.out.zscore_data,
        PCA_TECH_COV.out.metadata,
        file("${projectDir}/scripts/remove_tech_covar.R"),
        params.top_n_tech_cov,
        params.output_dir
    )
    
    // Stage 3: Cell type proportion estimation
    CELL_TYPE_DECONV (
        REMOVE_TECH_COVAR.out.corrected_data,
        REMOVE_TECH_COVAR.out.metadata_cleaned,
        file("${projectDir}/scripts/cell_type_deconv.R"),
        params.deconv_tool,
        params.reference_taxonomy,
        params.marker_file,
        params.hgnc_mapping_file,
        params.output_dir
    )
    
    // Stage 4: Prepare phenotype and covariate files for GWAS
    // Handle optional biospecimen_file - create empty file if not provided
    biospecimen_ch = params.biospecimen_file && params.biospecimen_file.toString() != '' && params.biospecimen_file.toString() != 'null'
        ? file(params.biospecimen_file) 
        : file("${projectDir}/.empty_biospecimen")  // Empty placeholder file
    
    PHENO_COV_PREP (
        CELL_TYPE_DECONV.out.cell_proportions,
        PCA_TECH_COV.out.metadata,
        file("${projectDir}/scripts/pheno_cov_prep.R"),
        file(params.wgs_psam_file),
        file(params.pca_file),
        biospecimen_ch,
        params.output_dir,
        params.study,  // Add study name parameter
        params.fid_method,
        params.fid_format,
        params.col_msex ?: '',
        params.biospec_col_individual ?: '',
        params.biospec_col_specimen ?: ''
    )
    
    emit:
    zscore_data = PCA_TECH_COV.out.zscore_data
    corrected_data = REMOVE_TECH_COVAR.out.corrected_data
    cell_proportions = CELL_TYPE_DECONV.out.cell_proportions
    phenotype_file = PHENO_COV_PREP.out.phenotype_file
    covariate_file = PHENO_COV_PREP.out.covariate_file
    samples_file = PHENO_COV_PREP.out.samples_file
}

// Main workflow entry point
workflow {
    PHENOTYPE_PIPELINE()
}

