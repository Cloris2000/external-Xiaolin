/*
 * Combined Pipeline V2: Phenotype Prep -> Genotyping QC -> Covariate Prep -> GWAS
 * 
 * FIXED: Eliminates circular dependency by splitting phenotype/covariate preparation
 * 
 * This unified pipeline orchestrates in CORRECT SCIENTIFIC ORDER:
 * 1. PHENOTYPE_PREP: Create phenotypes + sample list (NO genotype data needed)
 * 2. GENOTYPING_QC_PIPELINE: QC genotyping data FILTERED to phenotyped samples + PCA
 * 3. COV_PREP: Merge PCA with clinical covariates
 * 4. GWAS_PIPELINE: Perform GWAS analysis
 * 
 * Execution Order (SEQUENTIAL - NOT PARALLEL):
 * 
 *   PHENOTYPE_PREP
 *           ↓
 *      (generates phenotypes, sample list, clinical covariates)
 *           ↓
 *   GENOTYPING_QC_PIPELINE (filtered to phenotyped samples)
 *           ↓
 *      (generates PCA on QC'd genotypes)
 *           ↓
 *   COV_PREP
 *           ↓
 *      (merges PCA + clinical covariates)
 *           ↓
 *   GWAS_PIPELINE
 *           ↓
 *   Meta-analysis results
 * 
 * KEY DIFFERENCE FROM V1:
 * - V1: Genotype QC first → causes circular dependency
 * - V2: Phenotype prep first → scientifically correct, no circularity
 * 
 * SCIENTIFIC JUSTIFICATION:
 * - QC metrics (MAF, HWE, heterozygosity) should reflect the ANALYSIS cohort
 * - Filtering genotypes to phenotyped samples BEFORE QC ensures proper population metrics
 * - Matches standard GWAS best practices (Anderson et al. 2010, Nature Protocols)
 * 
 * Supports multiple cohorts: ROSMAP, Mayo, MSBB, GTEx, NABEC, CMC (MSSM, PITT, PENN)
 */

// Include workflows and modules
include { GENOTYPING_QC_PIPELINE } from './genotyping_qc_pipeline.nf'
include { GWAS_PIPELINE } from './gwas_pipeline.nf'

// Include the phenotype pipeline components (RNA-seq processing)
include { PCA_TECH_COV } from './modules/pca_tech_cov.nf'
include { REMOVE_TECH_COVAR } from './modules/remove_tech_covar.nf'
include { CELL_TYPE_DECONV } from './modules/cell_type_deconv.nf'

// Include the NEW split modules
include { PHENO_PREP } from './modules/pheno_prep.nf'
include { COV_PREP } from './modules/cov_prep.nf'

workflow COMBINED_PIPELINE {
    
    main:
    
    println """
    
    ╔════════════════════════════════════════════════════════════════╗
    ║     COMBINED PIPELINE V2 - SCIENTIFICALLY CORRECT ORDER        ║
    ╠════════════════════════════════════════════════════════════════╣
    ║  Stage 1: RNA-seq Processing (PCA, batch correction, deconv)  ║
    ║  Stage 2: Phenotype Prep (NO genotype data)                    ║
    ║  Stage 3: Genotype QC (filtered to phenotyped samples) + PCA   ║
    ║  Stage 4: Covariate Prep (merge PCA + clinical)                ║
    ║  Stage 5: GWAS Analysis                                        ║
    ╠════════════════════════════════════════════════════════════════╣
    ║  KEY FIX: Phenotypes created BEFORE genotype QC                ║
    ║           → No circular dependency                             ║
    ║           → QC metrics reflect analysis cohort                 ║
    ║           → Matches published best practices                   ║
    ╚════════════════════════════════════════════════════════════════╝
    
    """
    
    // ============================================================================
    // STAGE 1: RNA-seq Processing Pipeline
    // ============================================================================
    
    println """
    ========================================
    Stage 1: RNA-seq Processing
    ========================================
    Processing expression data...
    - PCA and technical covariate identification
    - Batch effect removal
    - Cell type deconvolution
    """
    
    // Stage 1a: PCA and technical covariate computation
    PCA_TECH_COV (
        params.count_matrix_file,
        params.metadata_file,
        file("${projectDir}/scripts/pca_tech_cov.R"),
        params.tissue_filter,
        params.output_dir,
        params.col_sample_id_for_matching
    )
    
    // Stage 1b: Remove batch effects and top tech covariates
    REMOVE_TECH_COVAR (
        PCA_TECH_COV.out.zscore_data,
        PCA_TECH_COV.out.metadata,
        file("${projectDir}/scripts/remove_tech_covar.R"),
        params.top_n_tech_cov,
        params.output_dir
    )
    
    // Stage 1c: Cell type proportion estimation
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
    
    // ============================================================================
    // STAGE 2: Phenotype Preparation (NO GENOTYPE DATA REQUIRED)
    // ============================================================================
    
    println """
    
    ========================================
    Stage 2: Phenotype Preparation
    ========================================
    Creating phenotypes and sample list...
    - RINT-transformed cell proportions
    - Sample list for genotype filtering
    - Clinical covariates (sex, age)
    
    NOTE: This step does NOT require genotype data!
    """
    
    // Handle optional biospecimen file
    biospecimen_ch = params.biospecimen_file && params.biospecimen_file.toString() != '' && params.biospecimen_file.toString() != 'null'
        ? file(params.biospecimen_file) 
        : file("${projectDir}/.empty_biospecimen")
    
    PHENO_PREP (
        CELL_TYPE_DECONV.out.cell_proportions,
        PCA_TECH_COV.out.metadata,
        file("${projectDir}/scripts/pheno_prep.R"),
        biospecimen_ch,
        params.output_dir,
        params.study,
        params.fid_method,
        params.fid_format,
        params.col_msex ?: '',
        params.col_age ?: '',
        params.col_individualID ?: '',
        params.clinical_metadata_file ?: '',
        params.biospec_col_individual ?: '',
        params.biospec_col_specimen ?: ''
    )
    
    // ============================================================================
    // STAGE 3: Genotyping QC Pipeline (FILTERED to phenotyped samples)
    // ============================================================================
    
    println """
    
    ========================================
    Stage 3: Genotyping QC Pipeline
    ========================================
    QC'ing genotypes FILTERED to phenotyped samples...
    - Sample filtering (using sample list from Stage 2)
    - Variant QC (MAF, HWE, missingness)
    - LD pruning
    - Heterozygosity filtering
    - PCA calculation
    
    IMPORTANT: QC metrics reflect the ANALYSIS cohort!
    """
    
    // Pass samples_with_phenotypes.txt to genotyping QC for filtering
    GENOTYPING_QC_PIPELINE(PHENO_PREP.out.samples_file)
    
    // ============================================================================
    // STAGE 4: Covariate Preparation (Merge PCA + Clinical)
    // ============================================================================
    
    println """
    
    ========================================
    Stage 4: Covariate Preparation
    ========================================
    Merging PCA with clinical covariates...
    - PCA from genotype QC (Stage 3)
    - Clinical covariates from phenotype prep (Stage 2)
    """
    
    COV_PREP (
        GENOTYPING_QC_PIPELINE.out.pca_file,
        PHENO_PREP.out.clinical_cov_file,
        PHENO_PREP.out.samples_file,
        file("${projectDir}/scripts/cov_prep.R"),
        params.output_dir,
        params.study
    )
    
    // ============================================================================
    // STAGE 5: GWAS Pipeline
    // ============================================================================
    
    println """
    
    ========================================
    Stage 5: GWAS Analysis
    ========================================
    Running GWAS for all cell types...
    - Using phenotypes from Stage 2
    - Using covariates from Stage 4
    - Using QC'd genotypes from Stage 3
    """
    
    // Convert file paths to file channels
    def pgen_file_ch = GENOTYPING_QC_PIPELINE.out.final_pgen.map { file(it) }
    def prune_in_file_ch = GENOTYPING_QC_PIPELINE.out.prune_in_file.map { file(it) }
    
    GWAS_PIPELINE(
        PHENO_PREP.out.phenotype_file,
        pgen_file_ch,
        prune_in_file_ch,
        COV_PREP.out.covariate_file
    )
    
    emit:
    // Genotyping QC pipeline outputs
    final_pgen = GENOTYPING_QC_PIPELINE.out.final_pgen
    final_psam = GENOTYPING_QC_PIPELINE.out.final_psam
    final_pvar = GENOTYPING_QC_PIPELINE.out.final_pvar
    pca_file = GENOTYPING_QC_PIPELINE.out.pca_file
    
    // Phenotype outputs
    phenotype_file = PHENO_PREP.out.phenotype_file
    samples_file = PHENO_PREP.out.samples_file
    clinical_cov_file = PHENO_PREP.out.clinical_cov_file
    
    // Covariate outputs
    covariate_file = COV_PREP.out.covariate_file
    
    // RNA-seq processing outputs
    cell_proportions = CELL_TYPE_DECONV.out.cell_proportions
    corrected_data = REMOVE_TECH_COVAR.out.corrected_data
    zscore_data = PCA_TECH_COV.out.zscore_data
    
    // GWAS pipeline outputs
    regenie_step1_pred = GWAS_PIPELINE.out.regenie_step1_pred
    regenie_step2_results = GWAS_PIPELINE.out.regenie_step2_results
    raw_p_files = GWAS_PIPELINE.out.raw_p_files
    meta_analysis_results = GWAS_PIPELINE.out.meta_analysis_results
}

// Main workflow entry point
workflow {
    COMBINED_PIPELINE()
}

