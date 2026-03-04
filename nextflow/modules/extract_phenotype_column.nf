/*
 * Module: Extract Phenotype Column
 * Extracts a specific cell type column from a multi-column phenotype file
 */

process EXTRACT_PHENOTYPE_COLUMN {
    label 'low_memory'
    tag "${cohort}_${cell_type}_extract"
    
    input:
    tuple val(cohort), val(cell_type), path(pheno_file)
    
    output:
    path "${cohort}_phenotypes_${cell_type}.txt", emit: pheno_file
    
    script:
    """
    # Extract FID, IID, and the specified cell type column
    awk -v col="${cell_type}" '
    BEGIN {
        FS="\\t"
        OFS="\\t"
    }
    NR==1 {
        # Find the column index for the cell type
        for (i=1; i<=NF; i++) {
            if (\$i == col) {
                col_idx = i
                break
            }
        }
        # Print header: FID, IID, and cell type name
        print \$1, \$2, col
    }
    NR>1 {
        # Print FID, IID, and the cell type column value
        print \$1, \$2, \$col_idx
    }' ${pheno_file} > ${cohort}_phenotypes_${cell_type}.txt
    """
}

