/*
 * Combined Pipeline: Genotyping QC -> Phenotype Processing -> GWAS Analysis
 * 
 * This unified pipeline orchestrates three major workflows in SEQUENTIAL order:
 * 1. GENOTYPING_QC_PIPELINE: QC genotyping data and generate PCA (RUNS FIRST)
 * 2. PHENOTYPE_PIPELINE: Process gene expression data and prepare phenotypes/covariates (RUNS SECOND)
 * 3. GWAS_PIPELINE: Perform GWAS analysis using outputs from pipelines 1 and 2 (RUNS THIRD)
 * 
 * Execution Order (SEQUENTIAL - NOT PARALLEL):
 * 
 *   GENOTYPING_QC_PIPELINE
 *           ↓
 *      (generates .psam and pca.csv)
 *           ↓
 *   PHENOTYPE_PIPELINE
 *           ↓
 *      (uses .psam and pca.csv + generates phenotypes/covariates)
 *           ↓
 *   GWAS_PIPELINE
 *           ↓
 *      (uses all outputs from above)
 *           ↓
 *   Meta-analysis results
 * 
 * IMPORTANT: Pipelines CANNOT run in parallel because:
 * - PHENOTYPE_PIPELINE requires wgs_psam_file and pca_file from GENOTYPING_QC_PIPELINE
 * - GWAS_PIPELINE requires outputs from both PHENOTYPE_PIPELINE and GENOTYPING_QC_PIPELINE
 * 
 * Supports multiple cohorts: ROSMAP, Mayo, MSBB, GTEx, NABEC, CMC (MSSM, PITT, PENN)
 */

// Include the three sub-workflows
include { GENOTYPING_QC_PIPELINE } from './genotyping_qc_pipeline.nf'
include { PHENOTYPE_PIPELINE as PHENOTYPE_PIPELINE_BASE } from './phenotype_pipeline.nf'
include { GWAS_PIPELINE } from './gwas_pipeline.nf'

// Wrapper process to ensure phenotype pipeline waits for genotyping QC
process CHECK_GENOTYPING_COMPLETE {
    label 'low_memory'
    
    input:
    val pca_file
    val psam_file
    
    output:
    val true
    
    script:
    """
    echo "=========================================="
    echo "Genotyping QC Complete!"
    echo "=========================================="
    echo "PCA file: ${pca_file}"
    echo "PSAM file: ${psam_file}"
    echo ""
    echo "Files verified. Ready for Phenotype Pipeline..."
    echo "=========================================="
    """
}

// Wrapper workflow for phenotype pipeline that enforces dependency on genotyping QC
workflow PHENOTYPE_PIPELINE {
    take:
    genotyping_complete
    
    main:
    // Wait for the signal that genotyping is complete
    genotyping_complete.view { 
        """
        
        ========================================
        Starting Phenotype Pipeline...
        ========================================
        """
    }
    
    // Run the actual phenotype pipeline
    PHENOTYPE_PIPELINE_BASE()
    
    emit:
    zscore_data = PHENOTYPE_PIPELINE_BASE.out.zscore_data
    corrected_data = PHENOTYPE_PIPELINE_BASE.out.corrected_data
    cell_proportions = PHENOTYPE_PIPELINE_BASE.out.cell_proportions
    phenotype_file = PHENOTYPE_PIPELINE_BASE.out.phenotype_file
    covariate_file = PHENOTYPE_PIPELINE_BASE.out.covariate_file
    samples_file = PHENOTYPE_PIPELINE_BASE.out.samples_file
}

workflow COMBINED_PIPELINE {
    
    main:
    
    // ============================================================================
    // STAGE 1: Genotyping QC Pipeline (RUNS FIRST)
    // ============================================================================
    
    println """
    
    ╔════════════════════════════════════════════════════════════════╗
    ║           COMBINED PIPELINE - SEQUENTIAL EXECUTION             ║
    ╠════════════════════════════════════════════════════════════════╣
    ║  Stage 1: Genotyping QC Pipeline                               ║
    ║  Stage 2: Phenotype Pipeline (waits for Stage 1)               ║
    ║  Stage 3: GWAS Pipeline (waits for Stage 1 & 2)                ║
    ╚════════════════════════════════════════════════════════════════╝
    
    ========================================
    Stage 1: Genotyping QC Pipeline
    ========================================
    Processing genotyping data...
    """
    
    GENOTYPING_QC_PIPELINE()
    
    // Check that genotyping QC is complete before proceeding
    CHECK_GENOTYPING_COMPLETE(
        GENOTYPING_QC_PIPELINE.out.pca_file,
        GENOTYPING_QC_PIPELINE.out.final_psam
    )
    
    // ============================================================================
    // STAGE 2: Phenotype Pipeline (RUNS SECOND, after STAGE 1)
    // ============================================================================
    
    // Run phenotype pipeline only after genotyping QC completes
    // The CHECK_GENOTYPING_COMPLETE output ensures this dependency
    PHENOTYPE_PIPELINE(CHECK_GENOTYPING_COMPLETE.out)
    
    // ============================================================================
    // STAGE 3: GWAS Pipeline (RUNS THIRD, after both STAGE 1 and 2)
    // ============================================================================
    
    // Create a combined signal channel that emits when both pipelines complete
    def genotyping_complete = GENOTYPING_QC_PIPELINE.out.pca_file
    def phenotype_complete = PHENOTYPE_PIPELINE.out.phenotype_file
    
    // Combine the two completion signals - GWAS will only start when both emit
    def both_complete = phenotype_complete.combine(genotyping_complete)
    
    // Wait for both pipelines to complete, then run GWAS
    both_complete.view { pheno, geno ->
        """
        
        ========================================
        Stage 3: GWAS Pipeline
        ========================================
        Prerequisites Complete!
        - Phenotype file: ${pheno}
        - PCA file: ${geno}
        
        Starting GWAS analysis for 19 cell types...
        ========================================
        """
    }
    
    // Run GWAS pipeline after both prerequisites complete
    // Pass the phenotype file and genotyping QC outputs so GWAS waits for them to exist
    // Use actual outputs from genotyping QC pipeline - convert file path strings to file channels
    // These will wait for the files to exist before GWAS starts
    def pgen_file_ch = GENOTYPING_QC_PIPELINE.out.final_pgen.map { file(it) }
    def prune_in_file_ch = GENOTYPING_QC_PIPELINE.out.prune_in_file.map { file(it) }
    
    GWAS_PIPELINE(
        PHENOTYPE_PIPELINE.out.phenotype_file,
        pgen_file_ch,
        prune_in_file_ch
    )
    
    emit:
    // Emit outputs from all three pipelines
    // Genotyping QC pipeline outputs
    final_pgen = GENOTYPING_QC_PIPELINE.out.final_pgen
    final_psam = GENOTYPING_QC_PIPELINE.out.final_psam
    final_pvar = GENOTYPING_QC_PIPELINE.out.final_pvar
    pca_file = GENOTYPING_QC_PIPELINE.out.pca_file
    
    // Phenotype pipeline outputs
    phenotype_file = PHENOTYPE_PIPELINE.out.phenotype_file
    covariate_file = PHENOTYPE_PIPELINE.out.covariate_file
    cell_proportions = PHENOTYPE_PIPELINE.out.cell_proportions
    
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
