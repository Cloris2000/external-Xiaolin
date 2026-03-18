#!/bin/bash
# Setup GTEx to use exact method from original script
# This creates a merged metadata file with specific ages

set -e

echo "=========================================="
echo "GTEx Setup - Using Exact Original Method"
echo "=========================================="
echo ""

echo "Step 1: Creating merged metadata file with specific ages..."
Rscript /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/scripts/prepare_gtex_metadata_with_ages.R

if [ ! -f "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/data/GTEx_FC_sample_metadata_with_ages.csv" ]; then
    echo "❌ ERROR: Failed to create merged metadata file"
    exit 1
fi

echo ""
echo "Step 2: Updating GTEx config to use merged metadata..."
echo "   - metadata_file → GTEx_FC_sample_metadata_with_ages.csv"
echo "   - col_msex → SEX.y (from merged file)"
echo "   - col_age → AGE.y (specific numeric ages from merged file)"
echo ""

echo "✅ Setup complete!"
echo ""
echo "Next steps:"
echo "  1. Update nextflow.config.combined.gtex:"
echo "     - metadata_file = '.../GTEx_FC_sample_metadata_with_ages.csv'"
echo "     - col_msex = 'SEX.y'"
echo "     - col_age = 'AGE.y'"
echo ""
echo "  2. Run pipeline:"
echo "     bash run_8_cohorts_truly_parallel.sh"
echo ""
