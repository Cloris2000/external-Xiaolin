#!/usr/bin/env nextflow

/*
 * Standalone Multi-Ancestry Meta-Analysis Pipeline
 *
 * Steps:
 *   1. Cohort QC audit  — produce a per-cohort audit matrix
 *   2. Harmonization    — remove strand-ambiguous SNPs, apply MAF and
 *                         min-presence filters, flag AF discrepancies
 *   3. METAL            — fixed-effect inverse-variance meta-analysis
 *                         with ANALYZE HETEROGENEITY (Q, I²)
 *
 * Primary analysis  : METAL SCHEME STDERR fixed-effect
 * Secondary/future  : MR-MEGA (run_mr_mega = true in config; disabled by default)
 *
 * Config entrypoints:
 *   nextflow.config.standalone_meta              — original 8-cohort meta
 *   nextflow.config.standalone_meta_all_cohorts  — all 12 cohorts including
 *                                                  HBCC platforms and GVEX
 */

nextflow.enable.dsl=2

include { META_COHORT_AUDIT }       from './modules/meta_cohort_audit.nf'
include { HARMONIZE_META_SUMSTATS } from './modules/harmonize_meta_sumstats.nf'
include { METAL_META_ANALYSIS }     from './modules/metal_meta_analysis.nf'
include { PLOT_META_RESULTS }       from './modules/plot_meta_results.nf'

workflow {

    // ------------------------------------------------------------------ //
    // Stage 0: Cohort QC audit
    // ------------------------------------------------------------------ //
    def audit_dir = "${params.audit_output_dir}"
    audit_ch = Channel.of(
        tuple(
            groovy.json.JsonOutput.toJson(params.cohorts),
            file("${params.cohort_metadata_file}"),
            file("${projectDir}/scripts/audit_meta_cohorts.py"),
            params.base_dir,
            audit_dir
        )
    )
    META_COHORT_AUDIT(audit_ch)

    // ------------------------------------------------------------------ //
    // Stage 1: Collect raw_p files, filter by cohort_include
    // ------------------------------------------------------------------ //
    raw_p_files = Channel
        .fromPath("${params.base_dir}/results/*/regenie_step2/*.regenie.raw_p")
        .map { f ->
            def cohort    = f.parent.parent.name
            def basename  = f.name.replaceAll(/\.regenie\.raw_p$/, "")
            def cell_type = basename.replaceAll(/${cohort}_/, "").replaceAll(/_step2$/, "")
            [cell_type, cohort, f]
        }
        .filter { cell_type, cohort, f ->
            (params.cohort_include == null || params.cohort_include.size() == 0) ||
            (cohort in params.cohort_include)
        }

    // Group by cell type: [cell_type, [files], [cohorts]]
    by_celltype = raw_p_files
        .groupTuple(by: 0)
        .map { cell_type, cohorts, files ->
            [cell_type, files instanceof List ? files : [files],
             cohorts instanceof List ? cohorts : [cohorts]]
        }

    // ------------------------------------------------------------------ //
    // Stage 2: Harmonization before METAL
    // ------------------------------------------------------------------ //
    harmonize_input = by_celltype.map { cell_type, files, cohorts ->
        tuple(
            cell_type,
            files,
            groovy.json.JsonOutput.toJson(cohorts),
            file("${projectDir}/scripts/harmonize_meta_sumstats.py"),
            "${params.harmonized_output_dir}/${cell_type}",
            params.drop_strand_ambiguous.toString(),
            params.min_meta_maf.toString(),
            params.min_present_cohorts.toString(),
            params.af_delta_report_threshold.toString(),
            params.af_delta_filter_threshold ? params.af_delta_filter_threshold.toString() : ""
        )
    }
    HARMONIZE_META_SUMSTATS(harmonize_input)

    // ------------------------------------------------------------------ //
    // Stage 3: METAL fixed-effect meta with heterogeneity output
    // ------------------------------------------------------------------ //
    def cohort_suffix = params.cohorts.sort().join('_')
    metal_input = HARMONIZE_META_SUMSTATS.out.harmonized_files
        .map { cell_type, harm_files ->
            def file_list = harm_files instanceof List ? harm_files : [harm_files]
            tuple(cell_type, file_list, params.metal_path, params.output_dir, cohort_suffix)
        }
    METAL_META_ANALYSIS(metal_input)

    // ------------------------------------------------------------------ //
    // Stage 4: Visualization — Manhattan, QQ, I² histograms, heatmap
    // Runs after ALL METAL jobs complete (collectFile triggers on last output)
    // ------------------------------------------------------------------ //
    // Collect a sentinel from all METAL jobs so plotting waits for all of them
    METAL_META_ANALYSIS.out.meta_result
        .collect()
        .map { _tbl_files ->
            tuple(
                params.output_dir,
                file("${projectDir}/scripts/plot_meta_results.R"),
                "${params.output_dir}/plots",
                params.plot_p_thresh  ?: "1e-5",
                params.plot_gw_thresh ?: "5e-8"
            )
        } | PLOT_META_RESULTS

    // ------------------------------------------------------------------ //
    // MR-MEGA placeholder (disabled by default; set run_mr_mega = true in config)
    // When enabled in a future step, this is where the MR-MEGA process
    // would consume the same harmonized files for a sensitivity run.
    // ------------------------------------------------------------------ //
}

workflow.onComplete {
    println ""
    println "=========================================="
    println "Multi-Ancestry Meta-Analysis Complete!"
    println "=========================================="
    println "Primary method    : METAL fixed-effect (SCHEME STDERR) + ANALYZE HETEROGENEITY"
    println "Secondary (future): MR-MEGA (run_mr_mega = ${params.run_mr_mega})"
    println "Results           : ${params.output_dir}"
    println "Harmonized inputs : ${params.harmonized_output_dir}"
    println "Cohort audit      : ${params.audit_output_dir}"
    println "=========================================="
}
