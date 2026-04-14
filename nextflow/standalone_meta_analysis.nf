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
 *   nextflow.config.standalone_meta_all_cohorts  — all 13 cohorts including
 *                                                  HBCC platforms, GVEX, and ROSMAP_array
 */

nextflow.enable.dsl=2

include { META_COHORT_AUDIT }       from './modules/meta_cohort_audit.nf'
include { HARMONIZE_META_SUMSTATS } from './modules/harmonize_meta_sumstats.nf'
include { METAL_META_ANALYSIS }     from './modules/metal_meta_analysis.nf'
include { PLOT_META_RESULTS }       from './modules/plot_meta_results.nf'
include { SUMMARIZE_META_HETEROGENEITY } from './modules/summarize_meta_heterogeneity.nf'
include { PREPARE_LDSC_SUMSTATS }   from './modules/prepare_ldsc_sumstats.nf'
include { RUN_LDSC_H2 }             from './modules/run_ldsc_h2.nf'
include { SUMMARIZE_LDSC_RESULTS }  from './modules/summarize_ldsc_results.nf'

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
        .fromPath("${params.meta_input_results_dir ?: "${params.base_dir}/results"}/*/regenie_step2/*.regenie.raw_p")
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
    // Stage 3b: Post-meta heterogeneity summary for MR-MEGA triage
    // ------------------------------------------------------------------ //
    if (params.run_heterogeneity_summary == null || params.run_heterogeneity_summary) {
        METAL_META_ANALYSIS.out.meta_result
            .collect()
            .map { tbl_files ->
                tuple(
                    tbl_files,
                    file("${projectDir}/scripts/summarize_meta_heterogeneity.py"),
                    params.heterogeneity_output_dir,
                    params.heterogeneity_gw_p_thresh ?: "5e-8",
                    params.heterogeneity_suggestive_p_thresh ?: "1e-5",
                    params.heterogeneity_i2_thresh ?: "50",
                    params.heterogeneity_het_p_thresh ?: "0.05"
                )
            } | SUMMARIZE_META_HETEROGENEITY
    } else {
        println "INFO: Skipping Stage 3b heterogeneity summary (run_heterogeneity_summary=false)"
    }

    if (params.run_ldsc) {
        def mergeAllelesPath = params.ldsc_merge_alleles ?: params.ldsc_hm3_snplist ?: ""
        if (!params.ldsc_ref_ld_chr || !params.ldsc_w_ld_chr) {
            throw new IllegalArgumentException(
                "run_ldsc=true requires ldsc_ref_ld_chr and ldsc_w_ld_chr"
            )
        }

        def ldsc_prep_input = METAL_META_ANALYSIS.out.meta_result_keyed
            .join(HARMONIZE_META_SUMSTATS.out.harmonized_files, by: 0)
            .map { cell_type, meta_tbl, harmonized_files ->
                tuple(
                    cell_type,
                    meta_tbl,
                    harmonized_files,
                    file("${projectDir}/scripts/prepare_ldsc_sumstats.py"),
                    params.ldsc_sumstats_outdir,
                    params.ldsc_variant_map ?: "",
                    params.ldsc_max_missing_variant_map_frac ?: ""
                )
            }
        PREPARE_LDSC_SUMSTATS(ldsc_prep_input)

        def ldsc_run_input = PREPARE_LDSC_SUMSTATS.out.ldsc_input
            .map { cell_type, ldsc_input_file ->
                tuple(
                    cell_type,
                    ldsc_input_file,
                    params.ldsc_conda_env ?: "",
                    mergeAllelesPath,
                    params.ldsc_ref_ld_chr,
                    params.ldsc_w_ld_chr,
                    params.ldsc_sumstats_outdir,
                    params.ldsc_results_outdir
                )
            }
        RUN_LDSC_H2(ldsc_run_input)

        RUN_LDSC_H2.out.h2_log
            .collect()
            .map { log_files ->
                tuple(
                    log_files,
                    file("${projectDir}/scripts/summarize_ldsc_results.py"),
                    params.ldsc_summary_outdir
                )
            } | SUMMARIZE_LDSC_RESULTS
    }

    // ------------------------------------------------------------------ //
    // Stage 4: Visualization — Manhattan, QQ, I² histograms, heatmap
    // Runs after ALL METAL jobs complete (collectFile triggers on last output)
    // ------------------------------------------------------------------ //
    if (params.run_plot_meta == null || params.run_plot_meta) {
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
    } else {
        println "INFO: Skipping Stage 4 plotting (run_plot_meta=false)"
    }

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
    println "LDSC post-meta    : ${params.run_ldsc ?: false}"
    println "Plotting enabled  : ${params.run_plot_meta == null || params.run_plot_meta}"
    println "Results           : ${params.output_dir}"
    println "Harmonized inputs : ${params.harmonized_output_dir}"
    println "Cohort audit      : ${params.audit_output_dir}"
    println "=========================================="
}
