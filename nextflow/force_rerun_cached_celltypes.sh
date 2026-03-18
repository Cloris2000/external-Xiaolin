#!/bin/bash
# Force rerun of the 9 cached cell types from December run
# This script removes the cached results to force Nextflow to regenerate them

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

echo "================================================"
echo "Force Rerun of 9 Cached Cell Types"
echo "================================================"
echo "Started at: $(date)"
echo ""

# List of cell types that were cached (from Dec 17)
CACHED_CELLTYPES=(
    "IT"
    "L4.IT"
    "L5.6.IT.Car3"
    "L5.ET"
    "L6b"
    "LAMP5"
    "Oligodendrocyte"
    "PAX6"
    "VLMC"
)

echo "Cell types to be rerun:"
for celltype in "${CACHED_CELLTYPES[@]}"; do
    echo "  - $celltype"
done
echo ""

# Backup the old results (optional safety measure)
BACKUP_DIR="results/meta_analysis_backup_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"

echo "Creating backup of cached results in: $BACKUP_DIR"
echo ""

# Remove the cached results to force rerun
echo "Removing cached results and work directories..."
for celltype in "${CACHED_CELLTYPES[@]}"; do
    echo "Processing: $celltype"
    
    # Backup the old results
    if [ -f "results/meta_analysis/${celltype}_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP1.tbl" ]; then
        cp "results/meta_analysis/${celltype}_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP1.tbl" "$BACKUP_DIR/"
        cp "results/meta_analysis/${celltype}_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP1.tbl.info" "$BACKUP_DIR/" 2>/dev/null
    fi
    
    # Remove the final output file to force rerun
    rm -f "results/meta_analysis/${celltype}_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP1.tbl"
    rm -f "results/meta_analysis/${celltype}_meta_analysis_CMC_MSSM_CMC_PENN_CMC_PITT_GTEx_MSBB_Mayo_NABEC_ROSMAP1.tbl.info"
    
    # Remove the .done marker if it exists
    rm -f "results/meta_analysis/${celltype}_meta.done"
    
    echo "  ✓ Removed cached results for $celltype"
done

echo ""
echo "================================================"
echo "Cleanup Complete!"
echo "================================================"
echo ""
echo "Now rerunning meta-analysis for these 9 cell types..."
echo ""

# Rerun the pipeline - it will now process the 9 cell types that had their results removed
nohup nextflow run standalone_meta_analysis.nf \
    -c nextflow.config.standalone_meta \
    -resume \
    > pipeline_rerun_cached_celltypes.log 2>&1 &

PIPELINE_PID=$!

echo "Pipeline restarted with PID: $PIPELINE_PID"
echo "Log file: pipeline_rerun_cached_celltypes.log"
echo ""
echo "To monitor progress:"
echo "  tail -f pipeline_rerun_cached_celltypes.log"
echo ""
echo "To check running jobs:"
echo "  squeue -u \$USER"
echo ""
echo "Backup of old results: $BACKUP_DIR"
echo "================================================"
