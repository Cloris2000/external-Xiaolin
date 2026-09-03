#!/bin/bash
#SBATCH --job-name=plot_meta
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --output=/scratch/zhoux156/logs/plot_meta_%j.out
#SBATCH --error=/scratch/zhoux156/logs/plot_meta_%j.err
#SBATCH --account=rrg-shreejoy

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
RESULTS_DIR="${NF_DIR}/results/meta_analysis_all_cohorts"
PLOTS_DIR="${NF_DIR}/results/meta_analysis_all_cohorts/plots"
SCRIPT="${NF_DIR}/scripts/plot_meta_results.R"

mkdir -p "${NF_DIR}/logs"
mkdir -p "${PLOTS_DIR}"

echo "============================================"
echo "Meta-Analysis Visualization"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Results    : ${RESULTS_DIR}"
echo "Output     : ${PLOTS_DIR}"
echo "Started    : $(date)"
echo "============================================"

# Verify .tbl files exist
N_TBL=$(ls "${RESULTS_DIR}"/*.tbl 2>/dev/null | grep -v "1\.tbl$" | wc -l)
if [[ ${N_TBL} -eq 0 ]]; then
    echo "ERROR: No .tbl files found in ${RESULTS_DIR}"
    exit 1
fi
echo "Found ${N_TBL} meta-analysis result files."

# Build list of cell types that still need plotting.
# "Complete" means all three per-cell outputs exist:
#   manhattan/<ct>_manhattan.png
#   manhattan/<ct>_qq.png
#   het/<ct>_het_histogram.png
mkdir -p "${PLOTS_DIR}/manhattan" "${PLOTS_DIR}/het" "${PLOTS_DIR}/heatmap"
missing_cts=()
while IFS= read -r tbl; do
    bn=$(basename "${tbl}")
    ct="${bn%%_meta_analysis_*}"
    man="${PLOTS_DIR}/manhattan/${ct}_manhattan.png"
    qq="${PLOTS_DIR}/manhattan/${ct}_qq.png"
    het="${PLOTS_DIR}/het/${ct}_het_histogram.png"
    if [[ -s "${man}" && -s "${qq}" && -s "${het}" ]]; then
        echo "[skip] ${ct} already completed"
    else
        missing_cts+=("${ct}")
    fi
done < <(ls "${RESULTS_DIR}"/*.tbl 2>/dev/null | grep -v "1\.tbl$")

echo "Cell types needing (re)run: ${#missing_cts[@]}"
if [[ ${#missing_cts[@]} -eq 0 ]]; then
    echo "All per-cell Manhattan/QQ/heterogeneity plots already complete."
else
    for ct in "${missing_cts[@]}"; do
        echo "Running plots for ${ct} ..."
        Rscript "${SCRIPT}" \
            --results_dir "${RESULTS_DIR}" \
            --output_dir  "${PLOTS_DIR}" \
            --cell_type   "${ct}" \
            --p_thresh    1e-5 \
            --gw_thresh   5e-8
    done
fi

echo ""
echo "============================================"
echo "Plotting complete: $(date)"
echo "Output directory : ${PLOTS_DIR}"
echo ""
echo "Contents:"
echo "  ${PLOTS_DIR}/manhattan/   — per-cell-type Manhattan + QQ plots"
echo "  ${PLOTS_DIR}/het/         — I² histograms + combined overview"
echo "  ${PLOTS_DIR}/heatmap/     — cross-cell-type P-value heatmap"
echo "  ${PLOTS_DIR}/meta_analysis_summary.tsv — λ and hit counts per cell type"
echo "============================================"
