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

include { META_COHORT_AUDIT }            from './modules/meta_cohort_audit.nf'
include { LIFTOVER_SUMSTATS }            from './modules/liftover_sumstats.nf'
include { HARMONIZE_META_SUMSTATS }      from './modules/harmonize_meta_sumstats.nf'
include { METAL_META_ANALYSIS }          from './modules/metal_meta_analysis.nf'
include { PLOT_META_RESULTS }            from './modules/plot_meta_results.nf'
include { PLOT_META_HEATMAP }            from './modules/plot_meta_heatmap.nf'
include { SUMMARIZE_META_HETEROGENEITY } from './modules/summarize_meta_heterogeneity.nf'
include { HETEROGENEITY_BREAKDOWN }      from './modules/heterogeneity_breakdown.nf'
include { PREPARE_LDSC_SUMSTATS }        from './modules/prepare_ldsc_sumstats.nf'
include { RUN_LDSC_H2 }                  from './modules/run_ldsc_h2.nf'
include { SUMMARIZE_LDSC_RESULTS }       from './modules/summarize_ldsc_results.nf'

workflow {

    // ------------------------------------------------------------------ //
    // Stage 0: Cohort QC audit (optional; disabled for per-cell-type array workers)
    // ------------------------------------------------------------------ //
    if (params.run_cohort_audit == null || params.run_cohort_audit) {
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
    } else {
        println "INFO: Skipping Stage 0 cohort audit (run_cohort_audit=false)"
    }

    // ------------------------------------------------------------------ //
    // Stages 1–2: Collect inputs → optional liftover → harmonization
    //   OR reuse pre-harmonized per-cohort files from a prior meta run.
    // ------------------------------------------------------------------ //
    def cohort_suffix = params.cohorts.sort().join('_')
    def harmonized_files_ch

    if (params.pre_harmonized_input_dir) {
        println "INFO: Using pre-harmonized inputs from ${params.pre_harmonized_input_dir} (skipping liftover + harmonization)"

        harmonized_files_ch = Channel
            .fromPath("${params.pre_harmonized_input_dir}/*/harmonized/*_harmonized.raw_p")
            .map { f ->
                def cell_type = f.parent.parent.name
                def basename  = f.name.replaceAll(/_harmonized\.raw_p$/, "")
                def cohort    = basename.replaceAll(/_${cell_type}$/, "")
                [cell_type, cohort, f]
            }
            .filter { cell_type, cohort, f ->
                def ctypes = params.cell_type_include ?: []
                if (ctypes instanceof String) {
                    ctypes = [ctypes]
                }
                def cohort_ok = (params.cohort_include == null || params.cohort_include.size() == 0) ||
                                (cohort in params.cohort_include)
                def ctype_ok  = (ctypes.size() == 0) || (cell_type in ctypes)
                cohort_ok && ctype_ok
            }
            .groupTuple(by: 0)
            .map { cell_type, cohorts, files ->
                [cell_type, files instanceof List ? files : [files]]
            }
    } else {
        raw_p_files = Channel
            .fromPath("${params.meta_input_results_dir ?: "${params.base_dir}/results"}/*/regenie_step2/*.regenie.raw_p")
            .map { f ->
                def cohort    = f.parent.parent.name
                def basename  = f.name.replaceAll(/\.regenie\.raw_p$/, "")
                def cell_type = basename.replaceAll(/${cohort}_/, "").replaceAll(/_step2$/, "")
                [cell_type, cohort, f]
            }
            .filter { cell_type, cohort, f ->
                def ctypes = params.cell_type_include ?: []
                if (ctypes instanceof String) {
                    ctypes = [ctypes]
                }
                def cohort_ok = (params.cohort_include == null || params.cohort_include.size() == 0) ||
                                (cohort in params.cohort_include)
                def ctype_ok  = (ctypes.size() == 0) || (cell_type in ctypes)
                cohort_ok && ctype_ok
            }

        def builds = params.cohort_genome_builds ?: [:]

        needs_liftover_ch = raw_p_files.filter { cell_type, cohort, f ->
            builds.get(cohort) == "hg38"
        }
        already_hg19_ch = raw_p_files.filter { cell_type, cohort, f ->
            builds.get(cohort) != "hg38"
        }

        if (builds.any { k, v -> v == "hg38" }) {
            if (!params.liftover_chain_hg38_to_hg19) {
                error "cohort_genome_builds contains hg38 cohorts but liftover_chain_hg38_to_hg19 is not set"
            }
            chain_file = file(params.liftover_chain_hg38_to_hg19)
            liftover_input_ch = needs_liftover_ch.map { cell_type, cohort, f ->
                tuple(
                    cell_type,
                    cohort,
                    f,
                    chain_file,
                    file("${projectDir}/scripts/liftover_sumstats.py"),
                    "${params.harmonized_output_dir}"
                )
            }
            LIFTOVER_SUMSTATS(liftover_input_ch)
            all_raw_p_files = already_hg19_ch.mix(LIFTOVER_SUMSTATS.out.lifted)
        } else {
            all_raw_p_files = already_hg19_ch
        }

        by_celltype = all_raw_p_files
            .groupTuple(by: 0)
            .map { cell_type, cohorts, files ->
                [cell_type, files instanceof List ? files : [files],
                 cohorts instanceof List ? cohorts : [cohorts]]
            }

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
        harmonized_files_ch = HARMONIZE_META_SUMSTATS.out.harmonized_files
    }

    // ------------------------------------------------------------------ //
    // Stage 3: METAL fixed-effect meta with heterogeneity output
    // ------------------------------------------------------------------ //
    metal_input = harmonized_files_ch
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

    // ------------------------------------------------------------------ //
    // Stage 3c: Heterogeneity breakdown across dimensions (optional)
    //
    // Decomposes per-variant I²/Q across cohort groups defined by:
    //   - diagnosis_group  (AD_neurological / psychiatric_mixed / neurologically_normal)
    //   - ancestry_class   (EUR_homogeneous / AFR_enriched / mixed_ancestry)
    //   - analysis_type    (DLPFC_array_imputed / DLPFC_WGS / non_DLPFC_brain / control_brain)
    //
    // Outputs forest plots per top hit + I² violin/heatmap summary figures.
    // Disabled by default (run_heterogeneity_breakdown = false).
    // ------------------------------------------------------------------ //
    if (params.run_heterogeneity_breakdown) {
        def het_bd_results_root  = params.meta_input_results_dir ?: "${params.base_dir}/results"
        def het_bd_output_dir    = params.heterogeneity_breakdown_dir ?: "${params.output_dir}/heterogeneity_breakdown"
        def het_bd_cohorts       = params.cohorts.sort().join(' ')
        def het_bd_pattern       = params.heterogeneity_breakdown_sumstat_pattern ?:
            "${het_bd_results_root}/{cohort}/regenie_step2/{cohort}_{cell_type}_step2.regenie.raw_p"

        METAL_META_ANALYSIS.out.meta_result
            .collect()
            .map { tbl_files ->
                tuple(
                    tbl_files,
                    file("${params.cohort_metadata_file}"),
                    file("${projectDir}/scripts/heterogeneity_breakdown.py"),
                    file("${projectDir}/scripts/plot_heterogeneity_breakdown.R"),
                    het_bd_output_dir,
                    het_bd_cohorts,
                    het_bd_results_root,
                    het_bd_pattern,
                    params.heterogeneity_gw_p_thresh ?: "5e-8",
                    params.heterogeneity_suggestive_p_thresh ?: "1e-5",
                    params.heterogeneity_breakdown_max_forest_plots ?: 50
                )
            } | HETEROGENEITY_BREAKDOWN
    } else {
        println "INFO: Skipping Stage 3c heterogeneity breakdown (run_heterogeneity_breakdown=false)"
    }

    if (params.run_ldsc) {
        def mergeAllelesPath = params.ldsc_merge_alleles ?: params.ldsc_hm3_snplist ?: ""
        if (!params.ldsc_ref_ld_chr || !params.ldsc_w_ld_chr) {
            throw new IllegalArgumentException(
                "run_ldsc=true requires ldsc_ref_ld_chr and ldsc_w_ld_chr"
            )
        }

        def ldsc_prep_input = METAL_META_ANALYSIS.out.meta_result_keyed
            .join(harmonized_files_ch, by: 0)
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
    // Stage 4: Visualization — 19 parallel per-CT jobs + 1 aggregate job
    //
    //   Stage 4a (PLOT_META_RESULTS × 19):
    //     One SLURM job per cell type — reads one .tbl file, produces
    //     Manhattan + QQ + I² histogram + 3 intermediate TSVs.
    //
    //   Stage 4b (PLOT_META_HEATMAP × 1):
    //     Runs after all 19 per-CT jobs complete — reads the intermediate
    //     TSVs (not the full .tbl files) to produce the cross-cell-type
    //     heatmap, combined I² overview, and meta_analysis_summary.tsv.
    // ------------------------------------------------------------------ //
    if (params.run_plot_meta == null || params.run_plot_meta) {
        def plot_script  = file("${projectDir}/scripts/plot_meta_results.R")
        def p_thresh     = params.plot_p_thresh  ?: "1e-5"
        def gw_thresh    = params.plot_gw_thresh ?: "5e-8"
        def plots_dir    = "${params.output_dir}/plots"

        // Stage 4a: one job per cell type (uses keyed output from METAL)
        per_ct_plot_input = METAL_META_ANALYSIS.out.meta_result_keyed
            .map { cell_type, tbl_file ->
                tuple(
                    cell_type,
                    tbl_file,
                    plot_script,
                    params.output_dir,   // results_dir: R script discovers .tbl from here
                    plots_dir,
                    p_thresh,
                    gw_thresh
                )
            }
        PLOT_META_RESULTS(per_ct_plot_input)

        // Stage 4b: aggregate — collect all intermediate TSVs from Stage 4a
        PLOT_META_RESULTS.out.top_hits
            .map { ct, f -> f }
            .mix(
                PLOT_META_RESULTS.out.het_stats.map    { ct, f -> f },
                PLOT_META_RESULTS.out.summary_stats.map { ct, f -> f }
            )
            .collect()
            .map { all_tsvs ->
                tuple(
                    all_tsvs,
                    plot_script,
                    plots_dir,
                    p_thresh,
                    gw_thresh
                )
            } | PLOT_META_HEATMAP
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
    println "Het breakdown     : ${params.run_heterogeneity_breakdown ?: false}"
    println "LDSC post-meta    : ${params.run_ldsc ?: false}"
    println "Plotting enabled  : ${params.run_plot_meta == null || params.run_plot_meta}"
    println "Results           : ${params.output_dir}"
    println "Harmonized inputs : ${params.harmonized_output_dir}"
    println "Cohort audit      : ${params.audit_output_dir}"
    println "=========================================="
}
