#!/bin/bash
#SBATCH --job-name=plot_meta_15cohorts
#SBATCH --partition=medium
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/meta_15cohorts/plot_meta_15cohorts_%j.out
#SBATCH --error=logs/meta_15cohorts/plot_meta_15cohorts_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
RESULTS_DIR=${NF_DIR}/results/meta_analysis_15cohorts
PLOTS_DIR=${RESULTS_DIR}/plots

echo "Checking meta-analysis results..."
tbl_count=$(ls ${RESULTS_DIR}/*.tbl 2>/dev/null | wc -l)
echo "  Found ${tbl_count}/19 .tbl files"

if [ "$tbl_count" -lt 19 ]; then
    echo "WARNING: Not all 19 cell types have completed. Plotting available cell types only."
fi

mkdir -p ${PLOTS_DIR}

echo ""
echo "Generating Manhattan + QQ plots for all available cell types..."
Rscript ${NF_DIR}/scripts/plot_meta_results.R \
    --results_dir ${RESULTS_DIR} \
    --output_dir  ${PLOTS_DIR} \
    --p_thresh    1e-5 \
    --gw_thresh   5e-8

echo ""
echo "Done. Plots written to: ${PLOTS_DIR}/"
ls ${PLOTS_DIR}/manhattan/*.png 2>/dev/null | wc -l | xargs -I{} echo "  {} Manhattan plots generated"
ls ${PLOTS_DIR}/qq/*.png 2>/dev/null | wc -l | xargs -I{} echo "  {} QQ plots generated"
