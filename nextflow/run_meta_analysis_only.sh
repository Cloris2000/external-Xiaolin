#!/bin/bash
# Run Meta-Analysis ONLY for all 8 cohorts
# This script uses existing GWAS results and runs only the meta-analysis step

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

echo "================================================"
echo "Running Meta-Analysis for All 8 Cohorts"
echo "================================================"
echo "Started at: $(date)"
echo "Cohorts: ROSMAP, Mayo, MSBB, CMC_MSSM, CMC_PENN, CMC_PITT, GTEx, NABEC"
echo "Cell types: 19"
echo "================================================"

# Check that all cohorts have GWAS results
echo "Checking existing GWAS results..."
all_good=true
for cohort in ROSMAP Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx NABEC; do
    count=$(ls results/$cohort/regenie_step2/*.regenie 2>/dev/null | wc -l)
    if [ "$count" -eq 19 ]; then
        echo "  ✓ $cohort: $count/19 regenie files"
    else
        echo "  ✗ $cohort: $count/19 regenie files (INCOMPLETE!)"
        all_good=false
    fi
done

if [ "$all_good" = false ]; then
    echo ""
    echo "ERROR: Some cohorts have incomplete GWAS results!"
    echo "Please complete GWAS for all cohorts first."
    exit 1
fi

echo ""
echo "All GWAS results present! Starting meta-analysis..."
echo "================================================"

# Run standalone meta-analysis pipeline - directly processes existing regenie files
nohup nextflow run standalone_meta_analysis.nf \
    -c nextflow.config.standalone_meta \
    -resume \
    > pipeline_run_meta_analysis_all_cohorts.log 2>&1 &

META_PID=$!

echo "Meta-analysis started with PID: $META_PID"
echo ""
echo "To monitor progress:"
echo "  tail -f pipeline_run_meta_analysis_all_cohorts.log"
echo ""
echo "To check SLURM jobs:"
echo "  squeue -u xzhou"
echo ""
echo "Results will be in:"
echo "  results/meta_analysis/"
echo ""
echo "================================================"
