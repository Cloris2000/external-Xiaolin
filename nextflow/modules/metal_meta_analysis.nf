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
ANALYZE
QUIT
EOF

    metal ${cell_type}_metal_script.txt
    
    # METAL appends a number (1, 2, etc.) when there's only one input file
    # Find the actual output file (could be ${cell_type}_meta_analysis_${cohort_suffix}.tbl or ${cell_type}_meta_analysis_${cohort_suffix}1.tbl, etc.)
    meta_output=\$(ls ${output_dir}/${cell_type}_meta_analysis_${cohort_suffix}*.tbl 2>/dev/null | head -1)
    meta_info=\$(ls ${output_dir}/${cell_type}_meta_analysis_${cohort_suffix}*.tbl.info 2>/dev/null | head -1)
    
    if [ -z "\$meta_output" ]; then
        echo "ERROR: METAL output file not found in ${output_dir}/" >&2
        ls -la ${output_dir}/${cell_type}_meta_analysis* 2>/dev/null || echo "No files found matching pattern"
        exit 1
    fi
    
    # Copy output files to work directory for Nextflow
    cp "\$meta_output" ${cell_type}_meta_analysis_${cohort_suffix}.tbl
    if [ -n "\$meta_info" ]; then
        cp "\$meta_info" ${cell_type}_meta_analysis_${cohort_suffix}.tbl.info
    fi
    
    touch ${cell_type}_meta.done
    """
}

