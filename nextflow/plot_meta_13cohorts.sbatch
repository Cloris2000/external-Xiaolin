#!/bin/bash
#SBATCH --job-name=plot_meta_13cohorts
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/scratch/zhoux156/logs/plot_meta_13cohorts_%j.out
#SBATCH --error=/scratch/zhoux156/logs/plot_meta_13cohorts_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

echo "================================================"
echo "Generating Plots for 13-Cohort Meta-Analysis"
echo "================================================"
echo "Results dir: ${NF_DIR}/results/meta_analysis_13cohorts"
echo "Output dir: ${NF_DIR}/results/meta_analysis_13cohorts/plots"
echo "================================================"
echo ""

/home/zhoux156/miniforge3/envs/test/bin/Rscript \
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
