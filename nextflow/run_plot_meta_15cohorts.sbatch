#!/bin/bash
#SBATCH --job-name=plot_meta_15cohorts
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=4:00:00
#SBATCH --output=logs/meta_15cohorts/plot_meta_15cohorts_%j.out
#SBATCH --error=logs/meta_15cohorts/plot_meta_15cohorts_%j.err

# Standalone re-plotting of an existing 15-cohort meta-analysis.
# Normally plots are produced in-pipeline by PLOT_META_RESULTS / PLOT_META_HEATMAP;
# use this only to regenerate them without re-running METAL.
#
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail

NF_DIR=/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow
RESULTS_DIR=${NF_DIR}/results/meta_analysis_15cohorts
PLOTS_DIR=${RESULTS_DIR}/plots
INTERMEDIATE_DIR=${PLOTS_DIR}/intermediate

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate test

echo "Checking meta-analysis results..."
tbl_count=$(ls ${RESULTS_DIR}/*.tbl 2>/dev/null | wc -l)
echo "  Found ${tbl_count}/19 .tbl files"

if [ "$tbl_count" -lt 19 ]; then
    echo "WARNING: Not all 19 cell types have completed. Plotting available cell types only."
fi

mkdir -p ${PLOTS_DIR} ${INTERMEDIATE_DIR}

# plot_meta_results.R requires --cell_type in its default per_celltype mode, so
# drive it once per cell type, then run the aggregate pass over the intermediates.
echo ""
echo "Generating Manhattan + QQ plots per cell type..."
while read -r cell_type; do
    [ -z "${cell_type}" ] && continue
    echo "  --> ${cell_type}"
    Rscript ${NF_DIR}/scripts/plot_meta_results.R \
        --mode        per_celltype \
        --results_dir ${RESULTS_DIR} \
        --cell_type   "${cell_type}" \
        --output_dir  ${PLOTS_DIR} \
        --p_thresh    1e-5 \
        --gw_thresh   5e-8
    for suffix in top_hits het_stats summary_stats; do
        [ -f "${PLOTS_DIR}/${cell_type}_${suffix}.tsv" ] && \
            cp "${PLOTS_DIR}/${cell_type}_${suffix}.tsv" "${INTERMEDIATE_DIR}/" || true
    done
done < ${NF_DIR}/bulk_rosmap_cell_types.txt

echo ""
echo "Generating cross-cell-type heatmap + I2 overview..."
Rscript ${NF_DIR}/scripts/plot_meta_results.R \
    --mode             aggregate \
    --intermediate_dir ${INTERMEDIATE_DIR} \
    --output_dir       ${PLOTS_DIR} \
    --p_thresh         1e-5 \
    --gw_thresh        5e-8

echo ""
echo "Done. Plots written to: ${PLOTS_DIR}/"
ls ${PLOTS_DIR}/manhattan/*.png 2>/dev/null | wc -l | xargs -I{} echo "  {} Manhattan plots generated"
ls ${PLOTS_DIR}/qq/*.png 2>/dev/null | wc -l | xargs -I{} echo "  {} QQ plots generated"
