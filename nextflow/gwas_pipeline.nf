/*
 * GWAS Pipeline for Cell Type Proportions
 * 
 * This pipeline performs:
 * 1. Regenie Step 1 (whole-genome regression) for each cohort and cell type
 * 2. Regenie Step 2 (association testing) for each cohort and cell type
 * 3. Convert log10 p-values to p-values for METAL compatibility
 * 4. METAL meta-analysis across cohorts for each cell type
 * 
 * Supports multiple cohorts: ROSMAP, Mayo, MSBB, GTEx, NABEC, CMC (MSSM, PITT, PENN)
 */

// Include modules
include { EXTRACT_PHENOTYPE_COLUMN } from './modules/extract_phenotype_column.nf'
include { REGENIE_STEP1 } from './modules/regenie_step1.nf'
include { REGENIE_STEP2 } from './modules/regenie_step2.nf'
include { CONVERT_LOGP_TO_P } from './modules/convert_logp_to_p.nf'
include { METAL_META_ANALYSIS } from './modules/metal_meta_analysis.nf'

workflow GWAS_PIPELINE {
    
    take:
    pheno_file_ch  // Channel containing the phenotype file (waits for phenotype pipeline)
    pgen_file_ch   // Channel containing the pgen file (waits for genotyping QC pipeline)
    prune_in_file_ch  // Channel containing the prune.in file (waits for genotyping QC pipeline)
    
    main:
    
    // Create sorted cohorts list once to avoid concurrent modification
    def sorted_cohorts = params.cohorts.sort { -it.length() }
    def cohort_suffix = params.cohorts.sort().join('_')
    
    // Define cell types
    def cell_types = [
        "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.ET", "L5.6.IT.Car3", 
        "L5.6.NP", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", 
        "Oligodendrocyte", "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC"
    ]
    
    // Cell type name mapping for different cohort versions
    def cell_type_mapping = [
        "Astrocyte": "Astrocyte",
        "Endothelial": "Endothelial",
        "IT": "IT",
        "L4.IT": "L4_IT",
        "L5.ET": "L5_ET",
        "L5.6.IT.Car3": "L5.6_IT_Car3",
        "L5.6.NP": "L5.6_NP",
        "L6.CT": "L6_CT",
        "L6b": "L6b",
        "LAMP5": "LAMP5",
        "Microglia": "Microglia",
        "OPC": "OPC",
        "Oligodendrocyte": "Oligodendrocyte",
        "PAX6": "PAX6",
        "PVALB": "PVALB",
        "Pericyte": "Pericyte",
        "SST": "SST",
        "VIP": "VIP",
        "VLMC": "VLMC"
    ]
    
    // Create channels for cohort-cell type combinations
    def cohort_ch = Channel.fromList(params.cohorts)
    def celltype_ch = Channel.fromList(cell_types)
    def cohort_celltype_ch = cohort_ch.combine(celltype_ch)
    
    // Stage 0: Extract phenotype columns for each cohort and cell type
    // Use the actual phenotype file from the input channel (waits for it to exist)
    def extract_pheno_input = cohort_celltype_ch
        .combine(pheno_file_ch)
        .map { cohort, cell_type, pheno_file ->
            // Use the actual file from the phenotype pipeline, not just the path
            tuple(cohort, cell_type, pheno_file)
        }
    
    // Extract phenotype columns - call process directly on tuple channel
    EXTRACT_PHENOTYPE_COLUMN(extract_pheno_input)
    
    // Create a keyed channel for extracted phenotype files
    def extracted_pheno_keyed = EXTRACT_PHENOTYPE_COLUMN.out.pheno_file
        .map { pheno_file ->
            // Extract cohort and cell_type from filename: {cohort}_phenotypes_{cell_type}.txt
            def file_name = pheno_file.toString()
            def basename = file_name.split('/').last()
            def parts = basename.replace('_phenotypes_', '|').replace('.txt', '').split('\\|')
            def cohort = parts[0]
            def cell_type = parts[1]
            [[cohort, cell_type], pheno_file]
        }
    
    // Stage 1: Regenie Step 1 for each cohort and cell type
    // Wait for genotyping QC outputs to be available
    // Combine phenotype channel with genotyping QC outputs
    // First, wait for pgen and prune_in files to be ready
    def genotyping_ready = pgen_file_ch.combine(prune_in_file_ch)
    
    def regenie_step1_input = extracted_pheno_keyed
        .combine(genotyping_ready)
        .map { key, pheno_file, pgen_file, prune_in_file ->
            def (cohort, cell_type) = key
            // Use the actual files from channels, not from params
            def pgen_prefix = pgen_file.toString().replace('.pgen', '')
            def pgen_files = [
                file("${pgen_prefix}.pgen"),
                file("${pgen_prefix}.pvar"),
                file("${pgen_prefix}.psam")
            ]
            // Get cohort-specific parameters for other files
            def cohort_params = params.cohort_configs[cohort]
            def covar_file = file("${cohort_params.covar_file}")
            // Create cohort-specific output directory
            def step1_dir = "${projectDir}/results/${cohort}/regenie_step1"
            
            tuple(cohort, cell_type, pgen_files, prune_in_file, pheno_file, covar_file, params.regenie_path, params.regenie_threads, params.regenie_bsize, step1_dir)
        }
    
    // Call Regenie Step 1 directly on tuple channel
    REGENIE_STEP1(regenie_step1_input)
    
    // Stage 2: Regenie Step 2 for each cohort and cell type
    // Create a channel with pred_list files keyed by cohort and cell_type
    // Include both newly generated and existing pred_list files
    def pred_list_keyed = REGENIE_STEP1.out.pred_list
        .map { pred_list_file ->
            def file_name = pred_list_file.toString()
            def basename = file_name.split('/').last()
            // Extract cohort and cell_type from filename: {cohort}_{cell_type}_step1_pred.list
            def parts = basename.replace('_step1_pred.list', '')
            // Split at first occurrence of cell type pattern (comes after cohort name)
            // Match cohort name from params.cohorts list (check longest names first)
            def cohort = null
            def cell_type = null
            for (c in sorted_cohorts) {
                if (parts.startsWith(c + '_')) {
                    cohort = c
                    cell_type = parts.substring(c.length() + 1)
                    break
                }
            }
            [[cohort, cell_type], pred_list_file]
        }
    
    // Join pred_list with extracted phenotype files
    def regenie_step2_input = pred_list_keyed
        .join(extracted_pheno_keyed, by: 0)
        .map { key, pred_list_file, extracted_pheno_file ->
            def (cohort, cell_type) = key
            def cohort_params = params.cohort_configs[cohort]
            def covar_file = file("${cohort_params.covar_file}")
            // Get all three pgen companion files
            def pgen_prefix = cohort_params.pgen_file.toString().replace('.pgen', '')
            def pgen_files = [
                file("${pgen_prefix}.pgen"),
                file("${pgen_prefix}.pvar"),
                file("${pgen_prefix}.psam")
            ]
            // Create cohort-specific output directory
            def step2_dir = "${projectDir}/results/${cohort}/regenie_step2"
            
            tuple(cohort, cell_type, pgen_files, extracted_pheno_file, covar_file, pred_list_file, params.regenie_path, params.regenie_threads, params.regenie_bsize, step2_dir)
        }
    
    // Call Regenie Step 2 directly on tuple channel
    REGENIE_STEP2(regenie_step2_input)
    
    // Stage 3: Convert logp to p for METAL compatibility
    def convert_input = REGENIE_STEP2.out.regenie_file
        .map { regenie_file ->
            def file_name = regenie_file.toString()
            // Extract cohort and cell_type from filename: {cohort}_{cell_type}_step2_{cell_type}.regenie
            def basename = file_name.split('/').last()
            def parts = basename.replace('_step2_', '|').replace('.regenie', '').split('\\|')
            def cohort_celltype = parts[0]
            // Match cohort name from params.cohorts list (check longest names first)
            def cohort = null
            def cell_type = null
            for (c in sorted_cohorts) {
                if (cohort_celltype.startsWith(c + '_')) {
                    cohort = c
                    cell_type = cohort_celltype.substring(c.length() + 1)
                    break
                }
            }
            // Create output directory for converted files (same as step2 output directory)
            // Safety check: ensure cohort and cell_type are not null
            if (cohort == null || cell_type == null) {
                throw new Exception("Failed to parse cohort and cell_type from filename: ${basename}")
            }
            def output_dir = "${projectDir}/results/${cohort}/regenie_step2"
            tuple(cohort, cell_type, regenie_file, output_dir)
        }
    
    // Call conversion process directly on tuple channel
    CONVERT_LOGP_TO_P(convert_input)
    
    // Stage 4: METAL meta-analysis for each cell type
    // Group raw_p files by cell type
    def raw_p_by_celltype = CONVERT_LOGP_TO_P.out.raw_p_file
        .map { raw_p_file ->
            def file_name = raw_p_file.toString()
            def basename = file_name.split('/').last()
            // Extract cohort and cell_type: {cohort}_{cell_type}_step2.regenie.raw_p
            def parts = basename.replace('_step2.regenie.raw_p', '')
            
            // Match cohort name from cohorts list (check longest names first to handle multi-part names like CMC_MSSM)
            def cohort = null
            def cell_type = null
            for (c in sorted_cohorts) {
                if (parts.startsWith(c + '_')) {
                    cohort = c
                    cell_type = parts.substring(c.length() + 1)
                    break
                }
            }
            
            // Map cell type name back to standard name (reverse mapping for CMC)
            def standard_cell_type = cell_type_mapping.find { k, v -> v == cell_type }?.key ?: cell_type
            [standard_cell_type, raw_p_file]
        }
        .groupTuple(by: 0)
        .map { cell_type, raw_p_files ->
            // Convert to list if needed
            def file_list = raw_p_files instanceof List ? raw_p_files : [raw_p_files]
            // Meta-analysis directory is shared across cohorts
            def meta_dir = "${projectDir}/results/meta_analysis"
            tuple(cell_type, file_list, params.metal_path, meta_dir, cohort_suffix)
        }
    
    // Call METAL meta-analysis
    METAL_META_ANALYSIS(raw_p_by_celltype)
    
    emit:
    regenie_step1_pred = REGENIE_STEP1.out.pred_list
    regenie_step2_results = REGENIE_STEP2.out.regenie_file
    raw_p_files = CONVERT_LOGP_TO_P.out.raw_p_file
    meta_analysis_results = METAL_META_ANALYSIS.out.meta_result
}

// Main workflow entry point
workflow {
    GWAS_PIPELINE()
}

