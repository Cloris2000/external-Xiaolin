/*
 * Genotyping QC Pipeline for GWAS Analysis
 * 
 * This pipeline processes VCF files (WGS or SNP array) through:
 * 1. Create pgen files - per chromosome (uses pre-normalized VCFs from vcf_pattern or normalized_vcf_dir)
 * 2. Merge pgen files across chromosomes
 * 3. Standard GWAS QC (MAF, HWE, missingness, imputation quality)
 * 4. LD pruning
 * 5. Heterozygosity calculation and filtering
 * 6. Generate final QC'd pgen/pvar/psam files
 * 7. Calculate PCA
 * 
 * Note: VCF normalization step removed - use pre-normalized VCFs.
 * Supports multiple cohorts: ROSMAP, Mayo, MSBB, GTEx, NABEC, CMC (MSSM, PITT, PENN)
 */

// Include module with different names for each step
include { GENOTYPING_QC_STEP as CREATE_PGEN } from './modules/genotyping_qc_step.nf'
include { GENOTYPING_QC_STEP as CREATE_PGEN_CHR } from './modules/genotyping_qc_step.nf'
include { GENOTYPING_QC_STEP as MERGE_PGEN } from './modules/genotyping_qc_step.nf'
include { GENOTYPING_QC_STEP as STANDARD_QC } from './modules/genotyping_qc_step.nf'
include { GENOTYPING_QC_STEP as LD_PRUNING } from './modules/genotyping_qc_step.nf'
include { GENOTYPING_QC_STEP as HETEROZYGOSITY } from './modules/genotyping_qc_step.nf'
include { GENOTYPING_QC_STEP as FINAL_QC } from './modules/genotyping_qc_step.nf'
include { GENOTYPING_QC_STEP as CALCULATE_PCA } from './modules/genotyping_qc_step.nf'

workflow GENOTYPING_QC_PIPELINE {
    
    take:
    samples_file_ch  // Optional: Channel containing samples_with_phenotypes.txt file
    
    main:
    
    // Helper: Provide defaults for optional parameters
    def vcf_dir = params.vcf_dir ?: ''
    def vcf_pattern = params.vcf_pattern ?: ''
    def normalized_vcf_dir = params.normalized_vcf_dir ?: ''
    def mach_r2_filter = params.mach_r2_filter ?: ''
    def vcf_min_gq = params.vcf_min_gq ?: ''
    def vcf_min_dp = params.vcf_min_dp ?: ''
    def vcf_min_qual = params.vcf_min_qual ?: ''
    def log_dir = params.log_dir ?: ''
    def output_dir = params.output_dir ?: ''
    
    // Prepare samples_to_keep: Convert file channel to path string for process calls
    def samples_to_keep_val = samples_file_ch.map { it.toString() }
    
    // Log sample filtering
    samples_to_keep_val.view { path ->
        """
        ========================================
        GENOTYPING QC: Filtering to phenotyped samples
        ========================================
        Sample filter file: ${path}
        This ensures QC metrics reflect the ANALYSIS cohort!
        ========================================
        """
    }
    
    // Stage 1: Create pgen files (per chromosome) - PARALLEL FOR ALL COHORTS
    // Uses pre-normalized VCFs from vcf_pattern or normalized_vcf_dir
    def start_marker = file("${projectDir}/.pipeline_start")
    
    def chrom_channel = params.include_x ? Channel.from(1..22, 'X') : Channel.from(1..22)
    
    // Combine chromosome with samples file path so each task gets the actual path (not a channel ref).
    // Fixes "DataflowVariable(value=null)" / syntax error in --samples_to_keep.
    def pgen_params = chrom_channel.combine(samples_to_keep_val).map { chrom, samples_path ->
        tuple(
            'create_pgen_single_chr',
            params.study,
            file("${projectDir}/scripts/genotyping_qc.py"),
            start_marker,  // dependency marker
            params.qc_version,
            params.work_dir,
            chrom.toString(),  // chromosome - varies per task
            vcf_dir,
            vcf_pattern,
            normalized_vcf_dir,  // Use pre-normalized VCFs when available
            params.plink_path,
            params.bcftools_path,
            params.maf_threshold,
            params.hwe_threshold,
            params.geno_threshold,
            params.mind_threshold,
            mach_r2_filter,
            vcf_min_gq,
            vcf_min_dp,
            vcf_min_qual,
            params.prune_window_size,
            params.prune_step_size,
            params.prune_r2_threshold,
            params.include_x,
            params.sort_vars,
            params.autosome_only,
            samples_path,  // resolved path string, not channel
            false,  // normalize_vcf (deprecated - always use pre-normalized)
            params.use_slurm,
            params.slurm_time,
            log_dir,
            output_dir
        )
    }
    
    // Call the process with each parameter tuple - Nextflow parallelizes automatically
    def create_pgen_output = CREATE_PGEN_CHR(
        pgen_params.map { it[0] },   // step
        pgen_params.map { it[1] },   // study
        pgen_params.map { it[2] },   // script_file
        pgen_params.map { it[3] },   // prev_output
        pgen_params.map { it[4] },   // qc_version
        pgen_params.map { it[5] },   // work_dir
        pgen_params.map { it[6] },   // chromosome
        pgen_params.map { it[7] },   // vcf_dir
        pgen_params.map { it[8] },   // vcf_pattern
        pgen_params.map { it[9] },   // normalized_vcf_dir
        pgen_params.map { it[10] },  // plink_path
        pgen_params.map { it[11] },  // bcftools_path
        pgen_params.map { it[12] },  // maf_threshold
        pgen_params.map { it[13] },  // hwe_threshold
        pgen_params.map { it[14] },  // geno_threshold
        pgen_params.map { it[15] },  // mind_threshold
        pgen_params.map { it[16] },  // mach_r2_filter
        pgen_params.map { it[17] },  // vcf_min_gq
        pgen_params.map { it[18] },  // vcf_min_dp
        pgen_params.map { it[19] },  // vcf_min_qual
        pgen_params.map { it[20] },  // prune_window_size
        pgen_params.map { it[21] },  // prune_step_size
        pgen_params.map { it[22] },  // prune_r2_threshold
        pgen_params.map { it[23] },  // include_x
        pgen_params.map { it[24] },  // sort_vars
        pgen_params.map { it[25] },  // autosome_only
        pgen_params.map { it[26] },  // samples_to_keep
        pgen_params.map { it[27] },  // normalize_vcf
        pgen_params.map { it[28] },  // use_slurm
        pgen_params.map { it[29] },  // slurm_time
        pgen_params.map { it[30] },  // log_dir
        pgen_params.map { it[31] }   // output_dir
    ).output_files.collect()  // Collect all chromosome completion markers
    
    // Stage 3: Merge pgen files across chromosomes (wait for create_pgen)
    def merge_output = MERGE_PGEN (
        'merge',
        params.study,
        file("${projectDir}/scripts/genotyping_qc.py"),
        create_pgen_output,  // Pass previous output to create dependency
        params.qc_version,  // Parameter version to invalidate cache when QC params change
        params.work_dir,
        '',  // chromosome (not used for merge)
        vcf_dir,
        vcf_pattern,
        normalized_vcf_dir,
        params.plink_path,
        params.bcftools_path,
        params.maf_threshold,
        params.hwe_threshold,
        params.geno_threshold,
        params.mind_threshold,
        mach_r2_filter,
        vcf_min_gq,
        vcf_min_dp,
        vcf_min_qual,
        params.prune_window_size,
        params.prune_step_size,
        params.prune_r2_threshold,
        params.include_x,
        params.sort_vars,
        params.autosome_only,
        samples_to_keep_val,
        params.normalize_vcf,
        params.use_slurm,
        params.slurm_time,
        log_dir,
        output_dir
    ).output_files
    
    // Stage 4: Standard GWAS QC (wait for merge)
    def standard_qc_output = STANDARD_QC (
        'standard_qc',
        params.study,
        file("${projectDir}/scripts/genotyping_qc.py"),
        merge_output,  // Pass previous output to create dependency
        params.qc_version,  // Parameter version to invalidate cache when QC params change
        params.work_dir,
        '',  // chromosome (not used)
        vcf_dir,
        vcf_pattern,
        normalized_vcf_dir,
        params.plink_path,
        params.bcftools_path,
        params.maf_threshold,
        params.hwe_threshold,
        params.geno_threshold,
        params.mind_threshold,
        mach_r2_filter,
        vcf_min_gq,
        vcf_min_dp,
        vcf_min_qual,
        params.prune_window_size,
        params.prune_step_size,
        params.prune_r2_threshold,
        params.include_x,
        params.sort_vars,
        params.autosome_only,
        samples_to_keep_val,
        params.normalize_vcf,
        params.use_slurm,
        params.slurm_time,
        log_dir,
        output_dir
    ).output_files
    
    // Stage 5: LD pruning (wait for standard_qc)
    def pruning_output = LD_PRUNING (
        'ld_pruning',
        params.study,
        file("${projectDir}/scripts/genotyping_qc.py"),
        standard_qc_output,  // Pass previous output to create dependency
        params.qc_version,  // Parameter version to invalidate cache when QC params change
        params.work_dir,
        '',  // chromosome (not used)
        vcf_dir,
        vcf_pattern,
        normalized_vcf_dir,
        params.plink_path,
        params.bcftools_path,
        params.maf_threshold,
        params.hwe_threshold,
        params.geno_threshold,
        params.mind_threshold,
        mach_r2_filter,
        vcf_min_gq,
        vcf_min_dp,
        vcf_min_qual,
        params.prune_window_size,
        params.prune_step_size,
        params.prune_r2_threshold,
        params.include_x,
        params.sort_vars,
        params.autosome_only,
        samples_to_keep_val,
        params.normalize_vcf,
        params.use_slurm,
        params.slurm_time,
        log_dir,
        output_dir
    ).output_files
    
    // Stage 6: Calculate heterozygosity (wait for pruning)
    def heterozygosity_output = HETEROZYGOSITY (
        'heterozygosity',
        params.study,
        file("${projectDir}/scripts/genotyping_qc.py"),
        pruning_output,  // Pass previous output to create dependency
        params.qc_version,  // Parameter version to invalidate cache when QC params change
        params.work_dir,
        '',  // chromosome (not used)
        vcf_dir,
        vcf_pattern,
        normalized_vcf_dir,
        params.plink_path,
        params.bcftools_path,
        params.maf_threshold,
        params.hwe_threshold,
        params.geno_threshold,
        params.mind_threshold,
        mach_r2_filter,
        vcf_min_gq,
        vcf_min_dp,
        vcf_min_qual,
        params.prune_window_size,
        params.prune_step_size,
        params.prune_r2_threshold,
        params.include_x,
        params.sort_vars,
        params.autosome_only,
        samples_to_keep_val,
        params.normalize_vcf,
        params.use_slurm,
        params.slurm_time,
        log_dir,
        output_dir
    ).output_files
    
    // Stage 7: Generate final QC'd files (wait for heterozygosity)
    def final_qc_output = FINAL_QC (
        'final_qc',
        params.study,
        file("${projectDir}/scripts/genotyping_qc.py"),
        heterozygosity_output,  // Pass previous output to create dependency
        params.qc_version,  // Parameter version to invalidate cache when QC params change
        params.work_dir,
        '',  // chromosome (not used)
        vcf_dir,
        vcf_pattern,
        normalized_vcf_dir,
        params.plink_path,
        params.bcftools_path,
        params.maf_threshold,
        params.hwe_threshold,
        params.geno_threshold,
        params.mind_threshold,
        mach_r2_filter,
        vcf_min_gq,
        vcf_min_dp,
        vcf_min_qual,
        params.prune_window_size,
        params.prune_step_size,
        params.prune_r2_threshold,
        params.include_x,
        params.sort_vars,
        params.autosome_only,
        samples_to_keep_val,
        params.normalize_vcf,
        params.use_slurm,
        params.slurm_time,
        log_dir,
        output_dir
    ).output_files
    
    // Stage 8: Calculate PCA (wait for final_qc)
    def pca_output = CALCULATE_PCA (
        'pca',
        params.study,
        file("${projectDir}/scripts/genotyping_qc.py"),
        final_qc_output,  // Pass previous output to create dependency
        params.qc_version,  // Parameter version to invalidate cache when QC params change
        params.work_dir,
        '',  // chromosome (not used)
        vcf_dir,
        vcf_pattern,
        normalized_vcf_dir,
        params.plink_path,
        params.bcftools_path,
        params.maf_threshold,
        params.hwe_threshold,
        params.geno_threshold,
        params.mind_threshold,
        mach_r2_filter,
        vcf_min_gq,
        vcf_min_dp,
        vcf_min_qual,
        params.prune_window_size,
        params.prune_step_size,
        params.prune_r2_threshold,
        params.include_x,
        params.sort_vars,
        params.autosome_only,
        samples_to_keep_val,
        params.normalize_vcf,
        params.use_slurm,
        params.slurm_time,
        log_dir,
        output_dir
    ).output_files
    
    // Output summary - emit actual file channels that wait for process completion
    emit:
    final_pgen = CALCULATE_PCA.out.output_files.map { file("${params.output_dir}/${params.study}.QC.final.pgen") }
    final_psam = CALCULATE_PCA.out.output_files.map { file("${params.output_dir}/${params.study}.QC.final.psam") }
    final_pvar = CALCULATE_PCA.out.output_files.map { file("${params.output_dir}/${params.study}.QC.final.pvar") }
    prune_in_file = LD_PRUNING.out.output_files.map { file("${params.output_dir}/${params.study}.QC.prune.in") }
    pca_file = CALCULATE_PCA.out.output_files.map { file("${params.output_dir}/pca.csv") }
}

// Default workflow entry point
workflow {
    GENOTYPING_QC_PIPELINE()
}
