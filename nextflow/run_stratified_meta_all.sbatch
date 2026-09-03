#!/bin/bash
# Submit all five stratified 15-cohort meta-analyses in parallel.
#
# Usage:
#   bash run_stratified_meta_all.sbatch
#
# After all jobs finish, compare against pooled meta:
#   sbatch run_compare_stratified_meta.sbatch

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd "${NF_DIR}"

for script in \
    run_stratified_meta_neurologically_normal.sbatch \
    run_stratified_meta_ad_neurological.sbatch \
    run_stratified_meta_psychiatric_mixed.sbatch \
    run_stratified_meta_eur_homogeneous.sbatch \
    run_stratified_meta_dlpfc_only.sbatch
do
    echo "Submitting ${script}..."
    sbatch "${NF_DIR}/${script}"
done

echo ""
echo "All stratified meta jobs submitted."
echo "Monitor: squeue -u \$USER | grep meta_strat"
echo "Compare when done: sbatch run_compare_stratified_meta.sbatch"
