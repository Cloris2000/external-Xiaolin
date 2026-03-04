/*
 * Module: Extract All Phenotype Columns
 * Extracts all cell type columns from a multi-column phenotype file in one pass
 * This is much faster than extracting columns one at a time
 */

process EXTRACT_ALL_PHENOTYPE_COLUMNS {
    label 'low_memory'
    tag "${cohort}_extract_all"
    
    input:
    val cohort
    val cell_types_str  // Comma-separated string of cell types
    path pheno_file  // Multi-column phenotype file
    
    output:
    path "${cohort}_phenotypes_*.txt", emit: pheno_files
    
    script:
    """
    # Extract all columns in one pass - much faster than 19 separate processes
    # Convert comma-separated string to space-separated for awk
    cell_types_space=\$(echo "${cell_types_str}" | tr ',' ' ')
    # Use awk to extract all columns at once
    awk -F'\t' -v OFS='\t' -v cell_types="\${cell_types_space}" '
    BEGIN {
        split(cell_types, ct_array, " ")
        # Read header to find column indices
        getline
        for (i=1; i<=NF; i++) {
            col_idx[\$i] = i
        }
        # Create output files for each cell type
        for (j=1; j<=length(ct_array); j++) {
            ct = ct_array[j]
            if (ct in col_idx) {
                idx = col_idx[ct]
                outfile = "${cohort}_phenotypes_" ct ".txt"
                print \$1, \$2, \$idx > outfile
            }
        }
    }
    {
        # Process data rows
        for (j=1; j<=length(ct_array); j++) {
            ct = ct_array[j]
            if (ct in col_idx) {
                idx = col_idx[ct]
                outfile = "${cohort}_phenotypes_" ct ".txt"
                print \$1, \$2, \$idx >> outfile
            }
        }
    }
    ' ${pheno_file}
    """
}

