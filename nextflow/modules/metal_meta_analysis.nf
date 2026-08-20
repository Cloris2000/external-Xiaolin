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
    // Create METAL script - raw_p_files is a tuple, need to convert to list
    def file_list = raw_p_files instanceof List ? raw_p_files : [raw_p_files]
    def process_commands = file_list.collect { file ->
        "PROCESS ${file}"
    }.join('\n')
    
    """
    mkdir -p ${output_dir}
    
    # Add METAL to PATH
    export PATH="${metal_path}:\$PATH"
    
    cat > ${cell_type}_metal_script.txt << 'EOF'
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
FLIP OFF
MARKER ID
ALLELE ALLELE0 ALLELE1
FREQ A1FREQ
EFFECT BETA
STDERR SE
PVAL P
${process_commands}
OUTFILE ${output_dir}/${cell_type}_meta_analysis_${cohort_suffix} .tbl
ANALYZE HETEROGENEITY
QUIT
EOF

    metal ${cell_type}_metal_script.txt
    
    # METAL OUTFILE + ANALYZE HETEROGENEITY writes prefix1.tbl (real data).
    # A bare prefix.tbl can be empty; never pick size-0 or lexicographic head -1.
    meta_output=\$(ls -S ${output_dir}/${cell_type}_meta_analysis_${cohort_suffix}*1.tbl ${output_dir}/${cell_type}_meta_analysis_${cohort_suffix}.tbl 2>/dev/null | while read f; do [ -s "\$f" ] && echo "\$f" && break; done)
    meta_info=\$(ls ${output_dir}/${cell_type}_meta_analysis_${cohort_suffix}*1.tbl.info ${output_dir}/${cell_type}_meta_analysis_${cohort_suffix}.tbl.info 2>/dev/null | head -1)
    
    if [ -z "\$meta_output" ]; then
        echo "ERROR: METAL output file not found (or all empty) in ${output_dir}/" >&2
        ls -la ${output_dir}/${cell_type}_meta_analysis* 2>/dev/null || echo "No files found matching pattern"
        exit 1
    fi
    echo "Using METAL output: \$meta_output (\$(wc -c < "\$meta_output") bytes)"
    
    # Copy output files to work directory for Nextflow
    cp "\$meta_output" ${cell_type}_meta_analysis_${cohort_suffix}.tbl
    if [ -n "\$meta_info" ]; then
        cp "\$meta_info" ${cell_type}_meta_analysis_${cohort_suffix}.tbl.info
    fi
    
    touch ${cell_type}_meta.done
    """
}

