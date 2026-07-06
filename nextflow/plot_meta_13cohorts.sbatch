#!/bin/bash
#SBATCH --job-name=plot_meta_13cohorts
#SBATCH --partition=medium
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/plot_meta_13cohorts_%j.out
#SBATCH --error=logs/plot_meta_13cohorts_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

echo "================================================"
echo "Generating Plots for 13-Cohort Meta-Analysis"
echo "================================================"
echo "Results dir: ${NF_DIR}/results/meta_analysis_13cohorts"
echo "Output dir: ${NF_DIR}/results/meta_analysis_13cohorts/plots"
echo "================================================"
echo ""

/nethome/kcni/xzhou/.anaconda3/envs/test/bin/Rscript \
    ${NF_DIR}/scripts/plot_meta_results.R \
    --results_dir ${NF_DIR}/results/meta_analysis_13cohorts \
    --output_dir ${NF_DIR}/results/meta_analysis_13cohorts/plots \
    --p_thresh 1e-5 \
    --gw_thresh 5e-8

echo ""
echo "================================================"
echo "Plotting complete!"
echo "Results saved to: ${NF_DIR}/results/meta_analysis_13cohorts/plots/"
echo "================================================"
