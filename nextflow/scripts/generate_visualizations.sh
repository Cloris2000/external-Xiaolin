#!/bin/bash
# Quick script to generate GWAS visualizations

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS_DIR="${PROJECT_DIR}/results"
OUTPUT_DIR="${PROJECT_DIR}/visualizations"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "Error: Rscript not found. Please install R."
    exit 1
fi

# Check if required R packages are installed
echo "Checking R packages..."
Rscript -e "if (!require('qqman')) install.packages('qqman', repos='https://cloud.r-project.org')" 2>&1 | grep -v "Loading required"

# Generate visualizations for ROSMAP
if [ -d "${RESULTS_DIR}/ROSMAP/regenie_step2" ]; then
    echo "Generating visualizations for ROSMAP..."
    Rscript "${SCRIPT_DIR}/generate_gwas_report.R" \
        --results_dir "${RESULTS_DIR}/ROSMAP/regenie_step2" \
        --meta_dir "${RESULTS_DIR}/meta_analysis" \
        --output_dir "${OUTPUT_DIR}/ROSMAP" \
        --cohort "ROSMAP"
else
    echo "Warning: ROSMAP results directory not found"
fi

# Generate visualizations for GTEx (when available)
if [ -d "${RESULTS_DIR}/GTEx/regenie_step2" ]; then
    echo "Generating visualizations for GTEx..."
    Rscript "${SCRIPT_DIR}/generate_gwas_report.R" \
        --results_dir "${RESULTS_DIR}/GTEx/regenie_step2" \
        --meta_dir "${RESULTS_DIR}/meta_analysis" \
        --output_dir "${OUTPUT_DIR}/GTEx" \
        --cohort "GTEx"
fi

echo "Visualization generation complete!"
echo "Results saved to: ${OUTPUT_DIR}"

