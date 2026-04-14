#!/bin/bash
#SBATCH --job-name=plot_meta_13cohorts
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=logs/plot_meta_13cohorts_%j.out
#SBATCH --error=logs/plot_meta_13cohorts_%j.err

set -euo pipefail

NF_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
RESULTS_DIR="${NF_DIR}/results/meta_analysis_13cohorts_design_matrix"
PLOTS_DIR="${NF_DIR}/results/meta_analysis_13cohorts_design_matrix/plots"
SCRIPT="${NF_DIR}/scripts/plot_meta_results.R"

mkdir -p "${NF_DIR}/logs"
mkdir -p "${PLOTS_DIR}/manhattan" "${PLOTS_DIR}/het" "${PLOTS_DIR}/heatmap"

echo "============================================"
echo "13-Cohort Meta-Analysis Visualization"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Results    : ${RESULTS_DIR}"
echo "Output     : ${PLOTS_DIR}"
echo "Started    : $(date)"
echo "============================================"

# Verify .tbl files exist
N_TBL=$(ls "${RESULTS_DIR}"/*.tbl 2>/dev/null | grep -v "1\.tbl$" | wc -l)
if [[ ${N_TBL} -eq 0 ]]; then
    echo "ERROR: No .tbl files found in ${RESULTS_DIR}"
    echo "Run run_meta_13cohorts.sbatch first."
    exit 1
fi
echo "Found ${N_TBL} meta-analysis result files."

# Collect cell types that still need (re)plotting
missing_cts=()
while IFS= read -r tbl; do
    bn=$(basename "${tbl}")
    ct="${bn%%_meta_analysis_*}"
    man="${PLOTS_DIR}/manhattan/${ct}_manhattan.png"
    qq="${PLOTS_DIR}/manhattan/${ct}_qq.png"
    het="${PLOTS_DIR}/het/${ct}_het_histogram.png"
    if [[ -s "${man}" && -s "${qq}" && -s "${het}" ]]; then
        echo "[skip] ${ct} — plots already complete"
    else
        missing_cts+=("${ct}")
    fi
done < <(ls "${RESULTS_DIR}"/*.tbl 2>/dev/null | grep -v "1\.tbl$")

echo "Cell types to plot: ${#missing_cts[@]}"
if [[ ${#missing_cts[@]} -eq 0 ]]; then
    echo "All plots already complete."
    exit 0
fi

# ── Parallel execution: one Rscript per cell type, 8 at a time ──
# 96 GB / 8 jobs ≈ 12 GB per process — well within limits.
plot_one() {
    local ct="$1"
    local NF_DIR="$2"
    local SCRIPT="$3"
    local RESULTS_DIR="$4"
    local PLOTS_DIR="$5"
    local SLURM_JOB_ID="$6"
    echo "[$(date +%H:%M:%S)] START ${ct}"
    Rscript "${SCRIPT}" \
        --results_dir "${RESULTS_DIR}" \
        --output_dir  "${PLOTS_DIR}" \
        --cell_type   "${ct}" \
        --p_thresh    1e-5 \
        --gw_thresh   5e-8 \
        2>&1 | tee "${NF_DIR}/logs/plot_${ct}_${SLURM_JOB_ID}.log"
    echo "[$(date +%H:%M:%S)] DONE  ${ct}"
}
export -f plot_one

printf '%s\n' "${missing_cts[@]}" | \
    parallel --jobs 8 --halt never \
        plot_one {} "${NF_DIR}" "${SCRIPT}" "${RESULTS_DIR}" "${PLOTS_DIR}" "${SLURM_JOB_ID}"

echo ""
echo "Per-cell plots complete: $(date)"

# ── Cross-cell-type heatmap: single run after all per-cell plots done ──
echo "Generating cross-cell-type heatmap and summary..."
Rscript "${SCRIPT}" \
    --results_dir "${RESULTS_DIR}" \
    --output_dir  "${PLOTS_DIR}" \
    --p_thresh    1e-5 \
    --gw_thresh   5e-8

echo ""
echo "============================================"
echo "All plots complete: $(date)"
echo "  ${PLOTS_DIR}/manhattan/  — Manhattan + QQ plots (19 cell types)"
echo "  ${PLOTS_DIR}/het/        — I² heterogeneity histograms"
echo "  ${PLOTS_DIR}/heatmap/    — cross-cell-type P-value heatmap"
echo "  ${PLOTS_DIR}/meta_analysis_summary.tsv"
echo "============================================"
