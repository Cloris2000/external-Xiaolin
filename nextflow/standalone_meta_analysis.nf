#!/usr/bin/env nextflow

/*
 * Standalone Meta-Analysis Pipeline
 * 
 * This pipeline takes existing regenie output files from multiple cohorts,
 * converts LOG10P to raw P-values, and performs meta-analysis using METAL.
 */

nextflow.enable.dsl=2

// Include modules
include { METAL_META_ANALYSIS } from './modules/metal_meta_analysis.nf'

workflow {
    
    // Define cell types
    def cell_types = [
        "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.ET", "L5.6.IT.Car3", 
        "L5.6.NP", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", 
        "Oligodendrocyte", "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC"
    ]
    
    // Create channel of existing raw_p files for each cohort and cell type
    // Use only .regenie.raw_p pattern to avoid duplicates (GTEx has both patterns)
    raw_p_files = Channel.fromPath("${params.base_dir}/results/*/regenie_step2/*.regenie.raw_p")
        .map { file ->
            def cohort = file.parent.parent.name
            def basename = file.name.replaceAll(/\.regenie\.raw_p$/, "")
            def cell_type = basename.replaceAll(/${cohort}_/, "").replaceAll(/_step2$/, "")
            [cell_type, file]
        }
    
    // Group by cell type for meta-analysis
    raw_p_by_celltype = raw_p_files
        .groupTuple()
        .map { cell_type, files ->
            def cohort_suffix = params.cohorts.sort().join('_')
            [cell_type, files, params.metal_path, params.output_dir, cohort_suffix]
        }
    
    // Run meta-analysis
    METAL_META_ANALYSIS(raw_p_by_celltype)
}

workflow.onComplete {
    println ""
    println "=========================================="
    println "Meta-Analysis Complete!"
    println "=========================================="
    println "Results are in: ${params.output_dir}"
    println "=========================================="
}
