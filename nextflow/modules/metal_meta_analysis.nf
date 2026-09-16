/*
 * Module: METAL Meta-Analysis
 * Performs meta-analysis using METAL for a specific cell type across multiple cohorts
 */

process METAL_META_ANALYSIS {
    label 'medium_memory'
    tag "${cell_type}_meta"
    
    input:
    tuple val(cell_type), path(raw_p_files), val(metal_path), val(output_dir), val(cohort_suffix)
    
    output:
    tuple val(cell_type), path("${cell_type}_meta_analysis_${cohort_suffix}.tbl"), emit: meta_result_keyed
    path "${cell_type}_meta_analysis_${cohort_suffix}.tbl", emit: meta_result
    path "${cell_type}_meta_analysis_${cohort_suffix}.tbl.info", emit: meta_info, optional: true
    path "${cell_type}_meta.done", emit: done_file
    
    publishDir "${output_dir}", mode: 'copy', overwrite: true
    
    script:
    // METAL allele/frequency orientation -- do not "restore" the old lines.
    //
    // REGENIE's effect allele is ALLELE1 (ALT): BETA is per-copy of ALLELE1 and
    // A1FREQ is ALLELE1's frequency.  METAL's directive is
    //     ALLELE <allele the EFFECT/FREQ columns describe> <other allele>
    // so ALLELE1 must come first.
    //
    // This block previously read "ALLELE ALLELE0 ALLELE1" plus "FLIP OFF".  Two
    // bugs, verified against this METAL build (2011-03-25) on toy and real data:
    //
    //   1. The allele order was reversed.
    //   2. This build has no OFF argument for FLIP; it matches the FLIP keyword
    //      and prints "## All effects will be flipped".  "FLIP OFF" turned
    //      flipping ON.
    //
    // The two cancelled for EFFECT (reversed order, then flipped back) but not
    // for FREQ, which FLIP does not touch.  So past runs carried correct effect
    // sizes, p-values and heterogeneity, but Freq1/MinFreq/MaxFreq were reported
    // as 1-AF.  Single-cohort identity test on real chr1 data under the settings
    // below: 162,728/162,728 variants exact on both frequency and effect.
    // Re-running a 3-cohort meta with this fix leaves Allele1, P-value and Effect
    // byte-identical for all 230,464 variants and complements Freq1 for 230,437.
    //
    // Note that METAL also applies an internal canonical allele sort before
    // writing (statgen/METAL issue #9): output pairs are only ever A/C, A/G, A/T,
    // C/G, T/C, T/G, and Allele1 is NOT necessarily the allele named first here.
    // When it reorders, it flips Effect and Freq1 to match, so the table stays
    // self-consistent -- but any downstream code must read the effect allele from
    // the Allele1 column rather than assuming it is the ALT of MarkerName.
    //
    // Create METAL script - raw_p_files is a tuple, need to convert to list
    def file_list = raw_p_files instanceof List ? raw_p_files : [raw_p_files]
    def process_commands = file_list.collect { file ->
        "PROCESS ${file}"
    }.join('\n')
    
    """
    # Add METAL to PATH
    export PATH="${metal_path}:\$PATH"
    
    cat > ${cell_type}_metal_script.txt << 'EOF'
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
MARKER ID
ALLELE ALLELE1 ALLELE0
FREQ A1FREQ
EFFECT BETA
STDERR SE
PVAL P
${process_commands}
OUTFILE metal_out .tbl
ANALYZE HETEROGENEITY
QUIT
EOF

    metal ${cell_type}_metal_script.txt

    # METAL must write inside the task work directory, never straight into the
    # publishDir: OUTFILE + ANALYZE HETEROGENEITY emits <prefix>1.tbl, so pointing
    # it at the shared results dir left both <name>1.tbl and the published
    # <name>.tbl there, and the old `ls -S` pick then resolved to whichever was
    # biggest across ALL previous runs rather than the one just written. That is
    # how the pre-liftover May 21 tables kept being re-published over later runs.
    # Here the glob can only ever match this task's own output.
    meta_output=\$(ls metal_out*.tbl 2>/dev/null | head -1)

    if [ -z "\$meta_output" ] || [ ! -s "\$meta_output" ]; then
        echo "ERROR: METAL produced no non-empty output for ${cell_type}" >&2
        ls -la metal_out* 2>/dev/null || echo "No metal_out* files written"
        exit 1
    fi
    echo "Using METAL output: \$meta_output (\$(wc -c < "\$meta_output") bytes)"

    mv "\$meta_output" ${cell_type}_meta_analysis_${cohort_suffix}.tbl
    if [ -s "\$meta_output.info" ]; then
        mv "\$meta_output.info" ${cell_type}_meta_analysis_${cohort_suffix}.tbl.info
    fi

    touch ${cell_type}_meta.done
    """
}

